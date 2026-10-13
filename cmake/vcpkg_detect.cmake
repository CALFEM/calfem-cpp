# Locates vcpkg and sets CMAKE_TOOLCHAIN_FILE to its toolchain. Must be included before project().
# Does nothing when a toolchain file is already given (-DCMAKE_TOOLCHAIN_FILE=..., a preset, or a parent project)
# or with -DCALFEM_VCPKG_DETECT=OFF. Search order: the VCPKG_ROOT environment variable, the vcpkg executable on
# PATH, VCPKG_INSTALLATION_ROOT (GitHub runners) and common installation directories.

option(CALFEM_VCPKG_DETECT "Locate vcpkg and use its toolchain file when none is given" ON)

function(_calfem_detect_vcpkg)
    set(candidates)
    if(DEFINED ENV{VCPKG_ROOT})
        list(APPEND candidates "$ENV{VCPKG_ROOT}")
    endif()
    find_program(_calfem_vcpkg_exe vcpkg NO_CACHE)
    if(_calfem_vcpkg_exe)
        # Resolve symlinks (e.g. /usr/local/bin/vcpkg) to reach the vcpkg checkout
        file(REAL_PATH "${_calfem_vcpkg_exe}" vcpkg_exe)
        get_filename_component(vcpkg_dir "${vcpkg_exe}" DIRECTORY)
        list(APPEND candidates "${vcpkg_dir}")
    endif()
    if(DEFINED ENV{VCPKG_INSTALLATION_ROOT})
        list(APPEND candidates "$ENV{VCPKG_INSTALLATION_ROOT}")
    endif()
    if(WIN32)
        list(APPEND candidates "$ENV{USERPROFILE}/vcpkg" "C:/vcpkg" "C:/src/vcpkg" "C:/dev/vcpkg")
    else()
        list(APPEND candidates "$ENV{HOME}/vcpkg" "/opt/vcpkg" "/usr/local/vcpkg")
    endif()

    foreach(root IN LISTS candidates)
        set(toolchain "${root}/scripts/buildsystems/vcpkg.cmake")
        if(EXISTS "${toolchain}")
            file(TO_CMAKE_PATH "${toolchain}" toolchain)
            message(STATUS "Using vcpkg toolchain: ${toolchain}")
            set(CMAKE_TOOLCHAIN_FILE "${toolchain}" CACHE FILEPATH "CMake toolchain file" FORCE)
            return()
        endif()
    endforeach()
    message(WARNING "vcpkg not found: set VCPKG_ROOT, put vcpkg on PATH or pass -DCMAKE_TOOLCHAIN_FILE. "
                    "Dependencies must then be found without vcpkg.")
endfunction()

# A vcpkg toolchain path that does not exist is searched for again; it is typically left in the cache by a
# configure with an unset VCPKG_ROOT ("$env{VCPKG_ROOT}/scripts/..." becomes "/scripts/...").
if(CMAKE_TOOLCHAIN_FILE MATCHES "scripts/buildsystems/vcpkg\\.cmake$" AND NOT EXISTS "${CMAKE_TOOLCHAIN_FILE}")
    message(STATUS "vcpkg toolchain ${CMAKE_TOOLCHAIN_FILE} not found, searching for vcpkg")
    set(_calfem_vcpkg_search ON)
elseif(NOT CMAKE_TOOLCHAIN_FILE)
    set(_calfem_vcpkg_search ON)
else()
    set(_calfem_vcpkg_search OFF)
endif()

if(CALFEM_VCPKG_DETECT AND _calfem_vcpkg_search AND CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    _calfem_detect_vcpkg()
endif()
