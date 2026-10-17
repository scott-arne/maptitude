# CMakeLists.txt
cmake_minimum_required(VERSION 3.16)

if(POLICY CMP0144)
    cmake_policy(SET CMP0144 NEW)
endif()

project(maptitude VERSION 0.6.2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# OpenEye's Windows SDK ships static libraries built against the release C
# runtime (/MD) only, and under _DEBUG Python's pyconfig.h links
# python3XX_d.lib, which the python.org installer does not ship. A build on the
# debug runtime (/MDd), which is what a Debug configuration defaults to, can
# link neither, so every configuration uses the release runtime; Debug keeps /Od
# and its debug information. A cache default, so it can still be overridden.
if(MSVC)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDLL" CACHE STRING
        "MSVC C runtime for every configuration")
endif()

# Options
option(MAPTITUDE_BUILD_TESTS "Build tests" ON)
option(MAPTITUDE_BUILD_PYTHON "Build Python bindings" ON)
option(MAPTITUDE_UNIVERSAL2 "Build universal2 binary for macOS" OFF)
option(MAPTITUDE_USE_STABLE_ABI "Use Python stable ABI" ON)

if(APPLE AND MAPTITUDE_UNIVERSAL2)
    set(CMAKE_OSX_ARCHITECTURES "x86_64;arm64" CACHE STRING "Build architectures" FORCE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "11.0" CACHE STRING "Minimum macOS version" FORCE)
endif()

# Fetch header-only dependencies
include(FetchContent)

FetchContent_Declare(
    nanoflann
    GIT_REPOSITORY https://github.com/jlblancoc/nanoflann.git
    GIT_TAG v1.5.5
    GIT_SHALLOW TRUE
)

FetchContent_Declare(
    cmake_openeye
    GIT_REPOSITORY https://github.com/scott-arne/cmake-openeye.git
    GIT_TAG v1.2.1
    GIT_SHALLOW TRUE
)

# nanoflann's examples and unit tests are not this project's to build or run.
# Its tests also register with this project's ctest, and their build step names
# no configuration, so under a multi-config generator a Release ctest run
# reports nanoflann's unit_tests_run as Not Run. nanoflann declares CMake 3.5,
# under which option() ignores a plain variable unless CMP0077 defaults to NEW.
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(NANOFLANN_BUILD_EXAMPLES OFF)
set(NANOFLANN_BUILD_TESTS OFF)

FetchContent_MakeAvailable(nanoflann cmake_openeye)
list(APPEND CMAKE_MODULE_PATH "${cmake_openeye_SOURCE_DIR}")

# Dependencies
if(MAPTITUDE_BUILD_PYTHON)
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
    # SWIG 4.4.0 picks Python 3.13-only C API calls from the header version
    # instead of Py_LIMITED_API, so its stable-ABI wrapper does not compile
    # against Python 3.13+ headers. 4.4.1 fixed it (swig/swig#3283).
    if(MAPTITUDE_USE_STABLE_ABI)
        find_package(SWIG 4.4.1)
    else()
        find_package(SWIG 4.0)
    endif()
    # The bindings below are skipped when SWIG is missing, and wheel builds
    # force MAPTITUDE_BUILD_PYTHON ON, so a too-old SWIG would otherwise produce a
    # wheel with no extension. Fail here instead, as oemmpa does.
    if(NOT SWIG_FOUND)
        message(FATAL_ERROR
            "MAPTITUDE_BUILD_PYTHON is ON but no suitable SWIG was found "
            "(the stable-ABI build needs 4.4.1 or newer; SWIG 4.4.0, which "
            "Ubuntu 26.04 ships, is too old). Install a newer SWIG (pip install "
            "swig==4.5.1) and pass -DSWIG_EXECUTABLE=<path>, deleting "
            "CMakeCache.txt first if an older SWIG is cached, or configure "
            "with -DMAPTITUDE_BUILD_PYTHON=OFF for a C++-only build.")
    endif()
endif()

# cmake-openeye forces OPENEYE_USE_SHARED OFF on Windows, where neither the SDK
# nor the openeye-toolkits wheel ships MSVC import libraries, but it does so
# inside find_package(OpenEye) below. Checked here first, the preset's ON would
# make a Windows configure require an openeye-toolkits the build never uses.
if(OPENEYE_USE_SHARED AND NOT WIN32)
    find_package(OpenEyePython REQUIRED)
endif()

find_package(OpenEye REQUIRED)
find_package(OpenMP)

# Library sources
set(MAPTITUDE_SOURCES
    src/Residue.cpp
    src/UnitCell.cpp
    src/SymOp.cpp
    src/ScatteringFactors.cpp
    src/Grid.cpp
    src/DensityCalculator.cpp
    src/Metric.cpp
    src/DensityScoreResult.cpp
    src/GridOps.cpp
    src/MapIO.cpp
    src/SpatialIndex.cpp
)

add_library(maptitude ${MAPTITUDE_SOURCES})

set_target_properties(maptitude PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(maptitude
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# PocketFFT is a vendored single header (third_party/pocketfft, BSD-3-Clause)
# rather than a fetched or system dependency: it has no releases to pin and the
# wheel builds run with FetchContent disconnected. It is PRIVATE because no
# public header names it.
#
# Multithreading is off because DensityCalculator only ever asks for one
# thread. The threaded build would add a process-wide thread pool, a
# pthread_atfork handler and a hardware_concurrency() call at load time --
# exactly the kind of global FFT state that made the previous FFTW dependency
# collide with the copy OpenEye statically embeds in liboegrid.
target_include_directories(maptitude
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third_party/pocketfft)
target_compile_definitions(maptitude PRIVATE POCKETFFT_NO_MULTITHREADING)

target_link_libraries(maptitude
    PUBLIC
        OpenEye::OEChem
        OpenEye::OEGrid
    PRIVATE
        nanoflann::nanoflann
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(maptitude PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(maptitude PRIVATE MAPTITUDE_USE_OPENMP)
endif()

# GCC 8's libstdc++ keeps std::filesystem in a separate archive; GCC 9 folded
# it into libstdc++ proper. A shared object may leave symbols unresolved that
# an executable may not, so without this the x86_64 wheel built in the Rocky 8
# container linked and then failed at import with an undefined
# std::filesystem::path::parent_path. PUBLIC because the library is static
# and the reference surfaces on whatever links it.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
    target_link_libraries(maptitude PUBLIC stdc++fs)
endif()

# Tests
if(MAPTITUDE_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Python bindings
if(MAPTITUDE_BUILD_PYTHON AND SWIG_FOUND)
    add_subdirectory(swig)
endif()

# Install
install(TARGETS maptitude
    EXPORT maptitude-targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(DIRECTORY include/maptitude DESTINATION include)

install(EXPORT maptitude-targets
    FILE MaptitudeTargets.cmake
    NAMESPACE Maptitude::
    DESTINATION lib/cmake/maptitude
)
