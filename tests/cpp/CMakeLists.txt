include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.15.2
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(googletest)

add_executable(maptitude_tests
    test_unit_cell.cpp
    test_symop.cpp
    test_residue.cpp
    test_scattering_factors.cpp
    test_density_score_result.cpp
    test_grid_ops.cpp
    test_interpolation.cpp
    test_spatial_index.cpp
    test_scoring_helpers.cpp
    test_fixtures.cpp
    test_metric_analytic.cpp
    test_metric_characterization.cpp
    test_density_calculator_characterization.cpp
    test_grid_ops_characterization.cpp
    test_grid_io.cpp
    test_map_io.cpp
    test_input_validation.cpp
)

target_compile_definitions(maptitude_tests
    PRIVATE
        MAPTITUDE_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../data"
        MAPTITUDE_TEST_ASSET_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../assets/mapq")

target_link_libraries(maptitude_tests
    PRIVATE
        maptitude
        GTest::gtest_main
)

target_include_directories(maptitude_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Standalone: prints pin_values.h. Deliberately not part of maptitude_tests and
# not registered with ctest -- it writes the golden file the tests read.
add_executable(maptitude_pin_generator generate_pins.cpp)
target_link_libraries(maptitude_pin_generator PRIVATE maptitude)
target_include_directories(maptitude_pin_generator PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

include(GoogleTest)
# A per-test timeout belongs here, not only on the CI ctest invocation: a
# deadlocked FFT -- which is what four of these tests did while a second,
# incompatible FFTW was reachable through OpenEye's static liboegrid -- is
# indistinguishable from a slow test, and ctest otherwise waits forever. The
# slowest test in this suite takes 3.3 s in a Debug build.
gtest_discover_tests(maptitude_tests PROPERTIES TIMEOUT 60)
