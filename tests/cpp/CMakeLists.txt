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
    test_spatial_index.cpp
    test_scoring_helpers.cpp
)

target_link_libraries(maptitude_tests
    PRIVATE
        maptitude
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(maptitude_tests)
