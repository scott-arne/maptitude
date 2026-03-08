# CMakeLists.txt
cmake_minimum_required(VERSION 3.16)

if(POLICY CMP0144)
    cmake_policy(SET CMP0144 NEW)
endif()

project(maptitude VERSION 0.1.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
option(MAPTITUDE_BUILD_TESTS "Build tests" ON)
option(MAPTITUDE_BUILD_PYTHON "Build Python bindings" ON)
option(MAPTITUDE_UNIVERSAL2 "Build universal2 binary for macOS" OFF)

if(APPLE AND MAPTITUDE_UNIVERSAL2)
    set(CMAKE_OSX_ARCHITECTURES "x86_64;arm64" CACHE STRING "Build architectures" FORCE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "11.0" CACHE STRING "Minimum macOS version" FORCE)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Dependencies
find_package(OpenEye REQUIRED)
find_package(FFTW3 REQUIRED)
find_package(OpenMP)

# Fetch header-only dependencies
include(FetchContent)

FetchContent_Declare(
    nanoflann
    GIT_REPOSITORY https://github.com/jlblancoc/nanoflann.git
    GIT_TAG v1.5.5
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(nanoflann)

if(MAPTITUDE_BUILD_PYTHON)
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
    find_package(SWIG 4.0)
endif()

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
    src/SpatialIndex.cpp
)

add_library(maptitude ${MAPTITUDE_SOURCES})

set_target_properties(maptitude PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(maptitude
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(maptitude
    PUBLIC
        OpenEye::OEChem
        OpenEye::OEGrid
    PRIVATE
        FFTW3::fftw3
        nanoflann::nanoflann
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(maptitude PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(maptitude PRIVATE MAPTITUDE_USE_OPENMP)
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
