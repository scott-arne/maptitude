include(OpenEyeSWIG)

openeye_add_swig_module(
    NAME maptitude
    SWIG_FILE maptitude.i
    LINK_LIBS maptitude
    PYTHON_OUTPUT_DIR ${CMAKE_SOURCE_DIR}/python/maptitude
    STABLE_ABI ${MAPTITUDE_USE_STABLE_ABI}
    SWIG_FLAGS "-doxygen"
    EXPECTED_LIB_VARS
        OECHEM_LIBRARY OESYSTEM_LIBRARY OEPLATFORM_LIBRARY
        OEMATH_LIBRARY OEGRID_LIBRARY OEFIZZCHEM_LIBRARY
        OEZSTD_LIBRARY
    INIT_PY ${CMAKE_SOURCE_DIR}/python/maptitude/__init__.py
)
