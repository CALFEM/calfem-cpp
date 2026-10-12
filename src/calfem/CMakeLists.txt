file(GLOB_RECURSE CALFEM_HEADERS CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/include/calfem/*.h)
file(GLOB_RECURSE CALFEM_SOURCES CONFIGURE_DEPENDS *.cpp)
file(GLOB_RECURSE CALFEM_PRIVATE_HEADERS CONFIGURE_DEPENDS detail/*.h mesh/*.h)

# Geometry and mesh generation are optional (CALFEM_WITH_MESH)
if(NOT CALFEM_WITH_MESH)
    list(FILTER CALFEM_SOURCES EXCLUDE REGEX "/mesh/|geometry\\.cpp$")
    list(FILTER CALFEM_PRIVATE_HEADERS EXCLUDE REGEX "/mesh/")
endif()

add_library(calfem ${CALFEM_HEADERS} ${CALFEM_PRIVATE_HEADERS} ${CALFEM_SOURCES})
add_library(calfem::calfem ALIAS calfem)

target_include_directories(calfem
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(calfem PUBLIC Eigen3::Eigen)

if(CALFEM_WITH_MESH)
    target_link_libraries(calfem PRIVATE CDT::CDT)
    target_compile_definitions(calfem PUBLIC CALFEM_WITH_MESH)
endif()

if(CALFEM_WITH_GMSH)
    target_link_libraries(calfem PRIVATE Gmsh::Gmsh)
    target_compile_definitions(calfem PUBLIC CALFEM_WITH_GMSH)
    if(Gmsh_HAS_CWRAP)
        target_compile_definitions(calfem PRIVATE CALFEM_GMSH_CWRAP)
    endif()
endif()
target_compile_features(calfem PUBLIC cxx_std_17)

set_target_properties(calfem PROPERTIES
    VERSION ${PROJECT_VERSION}
    DEBUG_POSTFIX d
    WINDOWS_EXPORT_ALL_SYMBOLS ON)

if(MSVC)
    target_compile_options(calfem PRIVATE /W4 /permissive- /bigobj)
else()
    target_compile_options(calfem PRIVATE -Wall -Wextra)
endif()

source_group(TREE ${PROJECT_SOURCE_DIR}/include/calfem PREFIX "Header Files" FILES ${CALFEM_HEADERS})
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} PREFIX "Source Files" FILES ${CALFEM_SOURCES} ${CALFEM_PRIVATE_HEADERS})
