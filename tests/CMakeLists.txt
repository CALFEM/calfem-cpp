find_package(Catch2 3 CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)

file(GLOB CALFEM_TEST_SOURCES CONFIGURE_DEPENDS test_*.cpp)

# reference/ holds a frozen copy of the calfem code ObjectiveFrame was built against
# (namespace calfem_ref), used to check that results are unchanged.

add_executable(calfem_tests ${CALFEM_TEST_SOURCES} reference/calfem_ref.cpp)
target_include_directories(calfem_tests PRIVATE reference)
target_link_libraries(calfem_tests PRIVATE calfem::calfem Catch2::Catch2WithMain nlohmann_json::nlohmann_json)

# Results computed with calfem-python (see reference_data/generate_calfem_python.py)
target_compile_definitions(calfem_tests PRIVATE
    CALFEM_PYTHON_REFERENCE="${CMAKE_CURRENT_SOURCE_DIR}/reference_data/calfem_python_reference.json"
    CALFEM_PYTHON_MESH_REFERENCE="${CMAKE_CURRENT_SOURCE_DIR}/reference_data/calfem_python_mesh_reference.json"
    CALFEM_PYTHON_VIS_REFERENCE="${CMAKE_CURRENT_SOURCE_DIR}/reference_data/calfem_python_vis_reference.json"
    CALFEM_PYTHON_VIS_FIELDS_REFERENCE="${CMAKE_CURRENT_SOURCE_DIR}/reference_data/calfem_python_vis_fields_reference.json")

include(Catch)
catch_discover_tests(calfem_tests)

# Header self-containment check: every public header must compile when included on its own.

file(GLOB_RECURSE CALFEM_PUBLIC_HEADERS CONFIGURE_DEPENDS RELATIVE ${PROJECT_SOURCE_DIR}/include
     ${PROJECT_SOURCE_DIR}/include/calfem/*.h)

set(CALFEM_HEADER_CHECK_SOURCES "")

foreach(header ${CALFEM_PUBLIC_HEADERS})
    string(MAKE_C_IDENTIFIER ${header} name)
    set(check_source ${CMAKE_CURRENT_BINARY_DIR}/header_check/${name}.cpp)
    file(CONFIGURE OUTPUT ${check_source} CONTENT "#include <${header}>\n")
    list(APPEND CALFEM_HEADER_CHECK_SOURCES ${check_source})
endforeach()

add_library(calfem_header_check OBJECT ${CALFEM_HEADER_CHECK_SOURCES})
target_link_libraries(calfem_header_check PRIVATE calfem::calfem)

# The Gmsh shared library must be next to the executables
if(CALFEM_WITH_GMSH AND Gmsh_RUNTIME)
    foreach(target calfem_tests)
        add_custom_command(TARGET ${target} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different ${Gmsh_RUNTIME} $<TARGET_FILE_DIR:${target}>)
    endforeach()
endif()

# Qt backend tests (separate executable, minimal platform: no windows appear; vcpkg's qtbase has no
# offscreen plugin on Windows)
if(CALFEM_WITH_QT)
    add_executable(calfem_qt_tests qt/test_vis_qt.cpp)
    target_include_directories(calfem_qt_tests PRIVATE ${PROJECT_SOURCE_DIR}/src/calfem_qt)
    target_link_libraries(calfem_qt_tests PRIVATE calfem::vis_qt Catch2::Catch2WithMain)
    calfem_deploy_qt_plugins(calfem_qt_tests)
    catch_discover_tests(calfem_qt_tests PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=minimal")
endif()
