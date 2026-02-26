# cmake/FindOpenEye.cmake
# Find OpenEye Toolkits installation
#
# This module finds the OpenEye C++ toolkits and creates imported targets.
#
# User can set these environment variables or CMake variables:
#   OPENEYE_ROOT or OE_DIR - Root directory of OpenEye installation (for headers)
#   OPENEYE_LIB_DIR        - Override library directory (e.g., from openeye-toolkits Python package)
#
# Options:
#   OPENEYE_USE_SHARED - Prefer shared libraries over static (default: OFF)
#                        Set to ON for wheels that depend on openeye-toolkits
#
# The following imported targets are created:
#   OpenEye::OEChem     - OEChem library (main chemistry library)
#   OpenEye::OESystem   - OESystem library
#   OpenEye::OEPlatform - OEPlatform library
#   OpenEye::OEMath     - OEMath library
#   OpenEye::OEGraphSim - OEGraphSim library (if available)
#   OpenEye::zstd       - Bundled zstd library (if available)
#
# The following variables are set:
#   OpenEye_FOUND       - TRUE if OpenEye was found
#   OpenEye_GraphSim_FOUND - TRUE if OEGraphSim was found
#   OpenEye_VERSION     - Version string (e.g., "2025.2.1")

option(OPENEYE_USE_SHARED "Prefer shared OpenEye libraries for dynamic linking" OFF)
set(OPENEYE_LIB_DIR "" CACHE PATH "Override OpenEye library directory (e.g., from openeye-toolkits Python package)")

# Look for the include directory
find_path(OPENEYE_INCLUDE_DIR
    NAMES openeye.h
    PATHS
        ${OPENEYE_ROOT}/include
        $ENV{OPENEYE_ROOT}/include
        $ENV{OE_DIR}/include
        /opt/openeye/include
        /usr/local/openeye/include
    PATH_SUFFIXES openeye
)

# Get the library directory - use override if provided, otherwise derive from include path
if(OPENEYE_LIB_DIR)
    message(STATUS "OpenEye: Using library directory override: ${OPENEYE_LIB_DIR}")
    set(_OPENEYE_LIB_SEARCH_PATHS ${OPENEYE_LIB_DIR})
elseif(OPENEYE_INCLUDE_DIR)
    get_filename_component(_DEFAULT_LIB_DIR "${OPENEYE_INCLUDE_DIR}/../lib" ABSOLUTE)
    set(_OPENEYE_LIB_SEARCH_PATHS ${_DEFAULT_LIB_DIR})
endif()

# Set library search order based on preference (save/restore to not affect other finds)
set(_SAVED_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
if(OPENEYE_USE_SHARED)
    # For shared linking, look for .dylib/.so first
    if(APPLE)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .dylib .a)
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES .so .a)
    endif()
    message(STATUS "OpenEye: Preferring shared libraries for dynamic linking")
endif()

# Helper macro to find OpenEye library, handling versioned names (e.g., liboechem-4.3.0.1.dylib)
macro(find_openeye_library VAR_NAME LIB_NAME)
    # First try to find versioned shared library in the override directory (openeye-toolkits Python package)
    if(OPENEYE_LIB_DIR AND OPENEYE_USE_SHARED)
        if(APPLE)
            file(GLOB _VERSIONED_LIB "${OPENEYE_LIB_DIR}/lib${LIB_NAME}-*.dylib")
        else()
            file(GLOB _VERSIONED_LIB "${OPENEYE_LIB_DIR}/lib${LIB_NAME}-*.so")
        endif()
        if(_VERSIONED_LIB)
            # Get the first match (should only be one)
            list(GET _VERSIONED_LIB 0 ${VAR_NAME})
            message(STATUS "OpenEye: Found versioned ${LIB_NAME}: ${${VAR_NAME}}")
        endif()
    endif()

    # Fall back to standard find_library if versioned library not found
    if(NOT ${VAR_NAME})
        find_library(${VAR_NAME}
            NAMES ${LIB_NAME}
            PATHS
                ${_OPENEYE_LIB_SEARCH_PATHS}
                ${OPENEYE_ROOT}/lib
                $ENV{OPENEYE_ROOT}/lib
                $ENV{OE_DIR}/lib
                /opt/openeye/lib
                /usr/local/openeye/lib
            NO_DEFAULT_PATH
        )
    endif()
endmacro()

# Find required libraries
find_openeye_library(OECHEM_LIBRARY oechem)
find_openeye_library(OESYSTEM_LIBRARY oesystem)
find_openeye_library(OEPLATFORM_LIBRARY oeplatform)
find_openeye_library(OEMATH_LIBRARY oemath)

# Find grid library (needed for OEScalarGrid implementation)
find_openeye_library(OEGRID_LIBRARY oegrid)

# Find fizzchem library (transitive dependency of oegrid)
find_openeye_library(OEFIZZCHEM_LIBRARY oefizzchem)

# Find optional libraries
find_openeye_library(OEGRAPHSIM_LIBRARY oegraphsim)

# Find bundled zstd library (OpenEye bundles this) - uses different naming
if(OPENEYE_LIB_DIR AND OPENEYE_USE_SHARED)
    file(GLOB _ZSTD_LIB "${OPENEYE_LIB_DIR}/libzstd*.dylib" "${OPENEYE_LIB_DIR}/libzstd*.so")
    if(_ZSTD_LIB)
        list(GET _ZSTD_LIB 0 OEZSTD_LIBRARY)
        message(STATUS "OpenEye: Found zstd: ${OEZSTD_LIBRARY}")
    endif()
endif()
if(NOT OEZSTD_LIBRARY)
    find_library(OEZSTD_LIBRARY
        NAMES zstd_static zstd
        PATHS
            ${_OPENEYE_LIB_SEARCH_PATHS}
            ${OPENEYE_ROOT}/lib
            $ENV{OPENEYE_ROOT}/lib
            $ENV{OE_DIR}/lib
            /opt/openeye/lib
            /usr/local/openeye/lib
        NO_DEFAULT_PATH
    )
endif()

# Restore CMAKE_FIND_LIBRARY_SUFFIXES before finding system libraries
set(CMAKE_FIND_LIBRARY_SUFFIXES ${_SAVED_CMAKE_FIND_LIBRARY_SUFFIXES})

# Find system zlib
find_package(ZLIB REQUIRED)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenEye
    REQUIRED_VARS
        OPENEYE_INCLUDE_DIR
        OECHEM_LIBRARY
        OESYSTEM_LIBRARY
        OEPLATFORM_LIBRARY
        OEMATH_LIBRARY
)

# Determine library type based on file extension
if(OpenEye_FOUND)
    # Check if library name contains .dylib or .so (handles versioned names like liboechem-4.3.0.1.dylib)
    get_filename_component(OECHEM_NAME "${OECHEM_LIBRARY}" NAME)
    if(OECHEM_NAME MATCHES "\\.dylib$" OR OECHEM_NAME MATCHES "\\.so$" OR OECHEM_NAME MATCHES "\\.so\\.")
        set(OPENEYE_LIBRARY_TYPE SHARED)
        message(STATUS "OpenEye: Using shared libraries (dynamic linking)")
    else()
        set(OPENEYE_LIBRARY_TYPE STATIC)
        message(STATUS "OpenEye: Using static libraries")
    endif()

    # Extract OEChem version from library name (e.g., liboechem-4.3.0.1.dylib)
    get_filename_component(OECHEM_NAME "${OECHEM_LIBRARY}" NAME)
    string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+" OpenEye_VERSION "${OECHEM_NAME}")
    if(NOT OpenEye_VERSION)
        string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" OpenEye_VERSION "${OECHEM_NAME}")
    endif()
    if(OpenEye_VERSION)
        message(STATUS "OpenEye: OEChem version ${OpenEye_VERSION}")
    else()
        string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" OpenEye_VERSION "${OPENEYE_INCLUDE_DIR}")
        if(OpenEye_VERSION)
            message(STATUS "OpenEye: Toolkit version ${OpenEye_VERSION} (from path)")
        endif()
    endif()
endif()

if(OpenEye_FOUND AND NOT TARGET OpenEye::OEChem)
    # Create imported target for zstd if found
    if(OEZSTD_LIBRARY AND NOT TARGET OpenEye::zstd)
        add_library(OpenEye::zstd UNKNOWN IMPORTED)
        set_target_properties(OpenEye::zstd PROPERTIES
            IMPORTED_LOCATION "${OEZSTD_LIBRARY}"
        )
    endif()

    # Create imported target for OEMath
    add_library(OpenEye::OEMath UNKNOWN IMPORTED)
    set_target_properties(OpenEye::OEMath PROPERTIES
        IMPORTED_LOCATION "${OEMATH_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
    )

    # OEPlatform depends on zlib and zstd
    add_library(OpenEye::OEPlatform UNKNOWN IMPORTED)
    set_target_properties(OpenEye::OEPlatform PROPERTIES
        IMPORTED_LOCATION "${OEPLATFORM_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
    )
    if(OEZSTD_LIBRARY)
        set_property(TARGET OpenEye::OEPlatform APPEND PROPERTY
            INTERFACE_LINK_LIBRARIES "OpenEye::zstd;ZLIB::ZLIB"
        )
    else()
        set_property(TARGET OpenEye::OEPlatform APPEND PROPERTY
            INTERFACE_LINK_LIBRARIES "ZLIB::ZLIB"
        )
    endif()

    add_library(OpenEye::OESystem UNKNOWN IMPORTED)
    set_target_properties(OpenEye::OESystem PROPERTIES
        IMPORTED_LOCATION "${OESYSTEM_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "OpenEye::OEPlatform"
    )

    # OEFizzChem (transitive dependency of OEGrid)
    if(OEFIZZCHEM_LIBRARY)
        add_library(OpenEye::OEFizzChem UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEFizzChem PROPERTIES
            IMPORTED_LOCATION "${OEFIZZCHEM_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem"
        )
    endif()

    # OEGrid depends on OESystem and OEFizzChem
    if(OEGRID_LIBRARY)
        add_library(OpenEye::OEGrid UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEGrid PROPERTIES
            IMPORTED_LOCATION "${OEGRID_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
        )
        if(OEFIZZCHEM_LIBRARY)
            set_property(TARGET OpenEye::OEGrid PROPERTY
                INTERFACE_LINK_LIBRARIES "OpenEye::OESystem;OpenEye::OEFizzChem"
            )
        else()
            set_property(TARGET OpenEye::OEGrid PROPERTY
                INTERFACE_LINK_LIBRARIES "OpenEye::OESystem"
            )
        endif()
    endif()

    # OEChem depends on OESystem and OEMath
    add_library(OpenEye::OEChem UNKNOWN IMPORTED)
    set_target_properties(OpenEye::OEChem PROPERTIES
        IMPORTED_LOCATION "${OECHEM_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "OpenEye::OESystem;OpenEye::OEMath"
    )

    if(OEGRAPHSIM_LIBRARY)
        add_library(OpenEye::OEGraphSim UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEGraphSim PROPERTIES
            IMPORTED_LOCATION "${OEGRAPHSIM_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem"
        )
        set(OpenEye_GraphSim_FOUND TRUE)
    endif()

    # Export the library type for use in other CMake files
    set(OpenEye_LIBRARY_TYPE ${OPENEYE_LIBRARY_TYPE} CACHE STRING "OpenEye library type (SHARED or STATIC)")
endif()

mark_as_advanced(
    OPENEYE_INCLUDE_DIR
    OECHEM_LIBRARY
    OESYSTEM_LIBRARY
    OEPLATFORM_LIBRARY
    OEGRAPHSIM_LIBRARY
    OEMATH_LIBRARY
    OEZSTD_LIBRARY
)
