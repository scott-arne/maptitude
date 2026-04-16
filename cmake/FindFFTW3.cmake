# cmake/FindFFTW3.cmake
# Find FFTW3 library
#
# This module finds the FFTW3 double-precision library and creates an imported target.
#
# User can set these environment variables or CMake variables:
#   FFTW3_ROOT - Root directory of FFTW3 installation
#
# The following imported targets are created:
#   FFTW3::fftw3 - FFTW3 double-precision library
#
# The following variables are set:
#   FFTW3_FOUND        - TRUE if FFTW3 was found
#   FFTW3_INCLUDE_DIRS - Include directories
#   FFTW3_LIBRARIES    - Library files

# Common install locations
set(_FFTW3_SEARCH_PATHS
    ${FFTW3_ROOT}
    $ENV{FFTW3_ROOT}
    /usr/local
    /opt/homebrew
    /usr
)

find_path(FFTW3_INCLUDE_DIR
    NAMES fftw3.h
    PATHS ${_FFTW3_SEARCH_PATHS}
    PATH_SUFFIXES include
)

# When FFTW3_USE_STATIC is set, restrict search to static libraries only.
# This is required for universal2 builds where the dynamic library may not
# contain all requested architectures.
if(FFTW3_USE_STATIC)
    set(_FFTW3_ORIG_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
    if(WIN32)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a)
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    endif()
endif()

find_library(FFTW3_LIBRARY
    NAMES fftw3
    PATHS ${_FFTW3_SEARCH_PATHS}
    PATH_SUFFIXES lib lib64
)

if(FFTW3_USE_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${_FFTW3_ORIG_SUFFIXES})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW3
    REQUIRED_VARS FFTW3_LIBRARY FFTW3_INCLUDE_DIR
)

if(FFTW3_FOUND AND NOT TARGET FFTW3::fftw3)
    if(FFTW3_USE_STATIC)
        add_library(FFTW3::fftw3 STATIC IMPORTED)
    else()
        add_library(FFTW3::fftw3 UNKNOWN IMPORTED)
    endif()
    set_target_properties(FFTW3::fftw3 PROPERTIES
        IMPORTED_LOCATION "${FFTW3_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR}"
    )
    set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})
    set(FFTW3_LIBRARIES ${FFTW3_LIBRARY})
endif()

mark_as_advanced(FFTW3_INCLUDE_DIR FFTW3_LIBRARY)
