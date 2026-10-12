# Finds the Gmsh SDK (C/C++ API).
#
# Set GMSH_ROOT (CMake or environment variable) to the SDK directory, i.e. the directory with
# include/gmshc.h and the library (lib/gmsh.lib, lib/gmsh.dll.lib, lib/libgmsh.so, ...). The
# official Gmsh SDK, Linux packages (libgmsh-dev) and the "gmsh" Python wheel all work.
#
# Defines
#   Gmsh_FOUND
#   Gmsh::Gmsh          imported target
#   Gmsh_HAS_CWRAP      TRUE if gmsh.h_cwrap (C++ API on top of the C API) is available
#   Gmsh_RUNTIME        shared library to deploy next to executables (Windows)

find_path(Gmsh_INCLUDE_DIR
    NAMES gmshc.h
    HINTS ${GMSH_ROOT} ENV GMSH_ROOT
    PATH_SUFFIXES include)

find_library(Gmsh_LIBRARY
    NAMES gmsh gmsh.dll libgmsh
    HINTS ${GMSH_ROOT} ENV GMSH_ROOT
    PATH_SUFFIXES lib Lib lib64)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Gmsh REQUIRED_VARS Gmsh_LIBRARY Gmsh_INCLUDE_DIR)

if(Gmsh_FOUND)
    set(Gmsh_HAS_CWRAP FALSE)
    if(EXISTS "${Gmsh_INCLUDE_DIR}/gmsh.h_cwrap")
        set(Gmsh_HAS_CWRAP TRUE)
    endif()

    if(WIN32)
        get_filename_component(_gmsh_libdir "${Gmsh_LIBRARY}" DIRECTORY)
        get_filename_component(_gmsh_root "${_gmsh_libdir}" DIRECTORY)
        file(GLOB _gmsh_dlls "${_gmsh_libdir}/gmsh*.dll" "${_gmsh_root}/bin/gmsh*.dll")
        list(GET _gmsh_dlls 0 Gmsh_RUNTIME)
    endif()

    if(NOT TARGET Gmsh::Gmsh)
        if(WIN32 AND Gmsh_RUNTIME)
            add_library(Gmsh::Gmsh SHARED IMPORTED)
            set_target_properties(Gmsh::Gmsh PROPERTIES
                IMPORTED_IMPLIB "${Gmsh_LIBRARY}"
                IMPORTED_LOCATION "${Gmsh_RUNTIME}"
                INTERFACE_COMPILE_DEFINITIONS "GMSH_DLL")
        else()
            add_library(Gmsh::Gmsh UNKNOWN IMPORTED)
            set_target_properties(Gmsh::Gmsh PROPERTIES IMPORTED_LOCATION "${Gmsh_LIBRARY}")
        endif()
        set_target_properties(Gmsh::Gmsh PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${Gmsh_INCLUDE_DIR}")
    endif()
endif()

mark_as_advanced(Gmsh_INCLUDE_DIR Gmsh_LIBRARY)
