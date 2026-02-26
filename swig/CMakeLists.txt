# swig/CMakeLists.txt
# CMake configuration for SWIG Python bindings

cmake_minimum_required(VERSION 3.16)

# Option for Python stable ABI (abi3) - builds wheel compatible with Python 3.10+
option(MAPTITUDE_USE_STABLE_ABI "Build with Python stable ABI (abi3) for cross-version compatibility" OFF)

# Find required packages
find_package(SWIG 4.0 REQUIRED)
if(MAPTITUDE_USE_STABLE_ABI)
    if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.26")
        find_package(Python3 COMPONENTS Interpreter Development.SABIModule REQUIRED)
    else()
        find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
    endif()
else()
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
endif()

include(${SWIG_USE_FILE})

# Set SWIG properties
set_source_files_properties(maptitude.i PROPERTIES
    CPLUSPLUS ON
    SWIG_MODULE_NAME _maptitude
)

# Set SWIG flags
set(CMAKE_SWIG_FLAGS "-doxygen")

# Create the SWIG library
# Use MODULE type instead of SHARED - on macOS this creates a bundle which is
# required for Python extension modules with Python 3.13+
swig_add_library(maptitude_python
    TYPE MODULE
    LANGUAGE python
    OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}
    SOURCES maptitude.i
)

# Configure stable ABI (abi3) if requested
if(MAPTITUDE_USE_STABLE_ABI)
    target_compile_definitions(maptitude_python PRIVATE
        Py_LIMITED_API=0x030A0000
    )
    message(STATUS "Building with Python stable ABI (abi3) for Python 3.10+")
endif()

# Set include directories
target_include_directories(maptitude_python PRIVATE
    ${Python3_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/include
)

# Platform-specific linking
if(APPLE)
    # On macOS, use -undefined dynamic_lookup to allow Python symbols to be
    # resolved at runtime by the Python interpreter that loads the module.
    # This is critical for Python 3.13+ multi-phase initialization.
    # Do NOT link against Python3::Python directly on macOS.
    target_link_libraries(maptitude_python
        PRIVATE
            maptitude
    )
    target_link_options(maptitude_python PRIVATE
        -undefined dynamic_lookup
    )
else()
    # On other platforms, link against Python normally
    target_link_libraries(maptitude_python
        PRIVATE
            maptitude
            Python3::Python
    )
endif()

# Set output name and properties
set_target_properties(maptitude_python PROPERTIES
    OUTPUT_NAME "_maptitude"
    PREFIX ""
)

# Platform-specific suffix
if(APPLE)
    set_target_properties(maptitude_python PROPERTIES SUFFIX ".so")
elseif(WIN32)
    set_target_properties(maptitude_python PROPERTIES SUFFIX ".pyd")
endif()

# Add RPATH for finding libraries
if(APPLE)
    if(OpenEye_LIBRARY_TYPE STREQUAL "SHARED")
        execute_process(
            COMMAND ${Python3_EXECUTABLE} -c "from openeye import libs; import os; print(os.path.basename(libs.FindOpenEyeDLLSDirectory()))"
            OUTPUT_VARIABLE OPENEYE_PLATFORM
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        set_target_properties(maptitude_python PROPERTIES
            INSTALL_RPATH "@loader_path;@loader_path/../openeye/libs/${OPENEYE_PLATFORM}"
            BUILD_WITH_INSTALL_RPATH TRUE
        )
        message(STATUS "Setting rpath for OpenEye Python package: @loader_path/../openeye/libs/${OPENEYE_PLATFORM}")
    else()
        set_target_properties(maptitude_python PROPERTIES
            INSTALL_RPATH "@loader_path"
            BUILD_WITH_INSTALL_RPATH TRUE
        )
    endif()
elseif(UNIX)
    if(OpenEye_LIBRARY_TYPE STREQUAL "SHARED")
        execute_process(
            COMMAND ${Python3_EXECUTABLE} -c "from openeye import libs; import os; print(os.path.basename(libs.FindOpenEyeDLLSDirectory()))"
            OUTPUT_VARIABLE OPENEYE_PLATFORM
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        set_target_properties(maptitude_python PROPERTIES
            INSTALL_RPATH "$ORIGIN:$ORIGIN/../openeye/libs/${OPENEYE_PLATFORM}"
            BUILD_WITH_INSTALL_RPATH TRUE
        )
    else()
        set_target_properties(maptitude_python PROPERTIES
            INSTALL_RPATH "$ORIGIN"
            BUILD_WITH_INSTALL_RPATH TRUE
        )
    endif()
endif()

# Define output directory for Python module
set(PYTHON_OUTPUT_DIR ${CMAKE_SOURCE_DIR}/python/maptitude)

# Custom command to copy the generated Python wrapper and shared library
add_custom_command(TARGET maptitude_python POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${PYTHON_OUTPUT_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/_maptitude.py
        ${PYTHON_OUTPUT_DIR}/maptitude.py
    COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_FILE:maptitude_python>
        ${PYTHON_OUTPUT_DIR}/$<TARGET_FILE_NAME:maptitude_python>
    COMMENT "Copying Python module to ${PYTHON_OUTPUT_DIR}"
)

# Also copy libmaptitude for standalone package
add_custom_command(TARGET maptitude_python POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_FILE:maptitude>
        ${PYTHON_OUTPUT_DIR}/$<TARGET_FILE_NAME:maptitude>
    COMMENT "Copying libmaptitude to ${PYTHON_OUTPUT_DIR}"
)

# Install targets for wheel building
install(TARGETS maptitude_python
    LIBRARY DESTINATION maptitude
    RUNTIME DESTINATION maptitude
)

# Install the SWIG-generated Python file
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/_maptitude.py
    DESTINATION maptitude
    RENAME maptitude.py
)

# Install the Python package __init__.py
install(FILES ${CMAKE_SOURCE_DIR}/python/maptitude/__init__.py
    DESTINATION maptitude
)

# Generate build info file with OpenEye version for runtime checking
if(DEFINED OPENEYE_TOOLKITS_VERSION AND OPENEYE_TOOLKITS_VERSION)
    set(_BUILD_INFO_VERSION "${OPENEYE_TOOLKITS_VERSION}")
else()
    set(_BUILD_INFO_VERSION "${OpenEye_VERSION}")
endif()
set(BUILD_INFO_CONTENT "# Auto-generated build information - do not edit
OPENEYE_BUILD_VERSION = \"${_BUILD_INFO_VERSION}\"
OPENEYE_LIBRARY_TYPE = \"${OpenEye_LIBRARY_TYPE}\"
")
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/_build_info.py "${BUILD_INFO_CONTENT}")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/_build_info.py
    DESTINATION maptitude
)
