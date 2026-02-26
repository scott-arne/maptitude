# tests/CMakeLists.txt
include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.14.0
)
FetchContent_MakeAvailable(googletest)

add_executable(maptitude_tests
    cpp/test_main.cpp
    cpp/test_unit_cell.cpp
    cpp/test_symop.cpp
    cpp/test_residue.cpp
    cpp/test_scattering_factors.cpp
    cpp/test_density_score_result.cpp
    cpp/test_grid_ops.cpp
    cpp/test_spatial_index.cpp
    cpp/test_scoring_helpers.cpp
)

target_link_libraries(maptitude_tests
    PRIVATE
        maptitude
        GTest::gtest
)

include(GoogleTest)
gtest_discover_tests(maptitude_tests)
