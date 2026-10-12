cmake_minimum_required(VERSION 3.24)

if(POLICY CMP0144)
    cmake_policy(SET CMP0144 NEW) # GMSH_ROOT is used by find_package(Gmsh)
endif()

project(calfem VERSION 0.1.0 DESCRIPTION "CALFEM finite element library for C++" LANGUAGES CXX)

# Tests and install rules default to ON only when calfem is the top-level project, so that
# consumers using add_subdirectory()/FetchContent (e.g. ObjectiveFrame) only get the library.
option(CALFEM_BUILD_TESTS "Build the calfem unit tests" ${PROJECT_IS_TOP_LEVEL})
option(CALFEM_BUILD_EXAMPLES "Build the calfem examples" ${PROJECT_IS_TOP_LEVEL})
option(CALFEM_INSTALL "Generate install and package export rules" ${PROJECT_IS_TOP_LEVEL})

# Geometry and mesh generation (requires the CDT library, vcpkg feature "mesh")
option(CALFEM_WITH_MESH "Build geometry and mesh generation" ${PROJECT_IS_TOP_LEVEL})

# Gmsh backend for mesh generation. Requires the Gmsh SDK (set GMSH_ROOT). Gmsh is GPL licensed:
# programs linked with calfem built with this option must comply with the GPL when distributed.
option(CALFEM_WITH_GMSH "Use the Gmsh library for mesh generation (GPL)" OFF)

# Interactive figure windows for calfem::vis (target calfem::vis_qt, vcpkg feature "qt")
option(CALFEM_WITH_QT "Build the Qt backend for calfem::vis (calfem::vis_qt)" OFF)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

if(PROJECT_IS_TOP_LEVEL)
    set(CMAKE_DEBUG_POSTFIX d)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

    # Executables (examples, tests) go to bin/<config>, e.g. bin/debug and bin/release.
    # The generator expression also keeps multi-config generators from adding their own subdirectory.
    if(NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/$<LOWER_CASE:$<CONFIG>>)
    endif()
endif()

find_package(Eigen3 CONFIG REQUIRED)

if(Eigen3_VERSION VERSION_LESS 3.4)
    message(FATAL_ERROR "calfem requires Eigen 3.4 or later (found ${Eigen3_VERSION})")
endif()

if(CALFEM_WITH_MESH)
    find_package(CDT CONFIG REQUIRED)
endif()

if(CALFEM_WITH_GMSH)
    if(NOT CALFEM_WITH_MESH)
        message(FATAL_ERROR "CALFEM_WITH_GMSH requires CALFEM_WITH_MESH")
    endif()
    find_package(Gmsh REQUIRED)
    message(STATUS "Gmsh: ${Gmsh_LIBRARY} (C API wrapper: ${Gmsh_HAS_CWRAP})")
endif()

if(CALFEM_WITH_QT)
    find_package(Qt6 REQUIRED COMPONENTS Widgets)
    message(STATUS "Qt: ${Qt6_VERSION}")
endif()

add_subdirectory(src/calfem)

if(CALFEM_WITH_QT)
    add_subdirectory(src/calfem_qt)
endif()

if(CALFEM_BUILD_EXAMPLES)
    add_subdirectory(src/examples)
endif()

if(CALFEM_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(CALFEM_INSTALL)
    include(GNUInstallDirs)
    include(CMakePackageConfigHelpers)

    set(CALFEM_INSTALL_TARGETS calfem)
    if(CALFEM_WITH_QT)
        list(APPEND CALFEM_INSTALL_TARGETS calfem_vis_qt)
    endif()

    install(TARGETS ${CALFEM_INSTALL_TARGETS} EXPORT calfemTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

    install(DIRECTORY include/calfem DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

    install(EXPORT calfemTargets
        FILE calfemTargets.cmake
        NAMESPACE calfem::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/calfem)

    configure_package_config_file(cmake/calfemConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/calfemConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/calfem)

    write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/calfemConfigVersion.cmake
        COMPATIBILITY SameMinorVersion)

    install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/calfemConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/calfemConfigVersion.cmake
        cmake/FindGmsh.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/calfem)
endif()
