/// Tier 2 pins for GridOps functions.
///
/// src/GridOps.cpp had 25% function coverage (three of four functions
/// untested) before this task. These pins make Task 12's error-handling
/// additions verifiable.
#include <gtest/gtest.h>

#include <memory>

#include "maptitude/GridOps.h"

#include "fixtures.h"
#include "grid_summary.h"
#include "pin_values.h"

using namespace Maptitude;
using namespace MaptitudeTest;

TEST(GridOpsCharacterizationTest, ScaleMap) {
    OESystem::OEScalarGrid grid = MakeRampGrid(2.0, 0.5);
    scale_map(grid, 1.5);

    const GridSummary s = Summarize(grid);
    ExpectPinned(s.sum, MaptitudePins::GRIDOPS_SCALE_SUM);
    ExpectPinned(s.sum_sq, MaptitudePins::GRIDOPS_SCALE_SUM_SQ);
    ExpectPinned(s.min, MaptitudePins::GRIDOPS_SCALE_MIN);
    ExpectPinned(s.max, MaptitudePins::GRIDOPS_SCALE_MAX);
    ExpectPinned(s.index_moment, MaptitudePins::GRIDOPS_SCALE_INDEX_MOMENT);
}

TEST(GridOpsCharacterizationTest, CombineAdd) {
    const OESystem::OEScalarGrid lhs = MakeRampGrid(2.0, 0.5);
    const OESystem::OEScalarGrid rhs = MakeUniformGrid(10.0f, 2.0, 0.5);

    std::unique_ptr<OESystem::OEScalarGrid> result(combine_maps(lhs, rhs, MapOp::ADD));
    ASSERT_NE(result, nullptr);

    const GridSummary s = Summarize(*result);
    ExpectPinned(s.sum, MaptitudePins::GRIDOPS_ADD_SUM);
    ExpectPinned(s.sum_sq, MaptitudePins::GRIDOPS_ADD_SUM_SQ);
    ExpectPinned(s.min, MaptitudePins::GRIDOPS_ADD_MIN);
    ExpectPinned(s.max, MaptitudePins::GRIDOPS_ADD_MAX);
    ExpectPinned(s.index_moment, MaptitudePins::GRIDOPS_ADD_INDEX_MOMENT);
}

TEST(GridOpsCharacterizationTest, CombineSubtract) {
    const OESystem::OEScalarGrid lhs = MakeRampGrid(2.0, 0.5);
    const OESystem::OEScalarGrid rhs = MakeUniformGrid(10.0f, 2.0, 0.5);

    std::unique_ptr<OESystem::OEScalarGrid> result(combine_maps(lhs, rhs, MapOp::SUBTRACT));
    ASSERT_NE(result, nullptr);

    const GridSummary s = Summarize(*result);
    ExpectPinned(s.sum, MaptitudePins::GRIDOPS_SUBTRACT_SUM);
    ExpectPinned(s.sum_sq, MaptitudePins::GRIDOPS_SUBTRACT_SUM_SQ);
    ExpectPinned(s.min, MaptitudePins::GRIDOPS_SUBTRACT_MIN);
    ExpectPinned(s.max, MaptitudePins::GRIDOPS_SUBTRACT_MAX);
    ExpectPinned(s.index_moment, MaptitudePins::GRIDOPS_SUBTRACT_INDEX_MOMENT);
}

TEST(GridOpsCharacterizationTest, CombineMin) {
    const OESystem::OEScalarGrid lhs = MakeRampGrid(2.0, 0.5);
    const OESystem::OEScalarGrid rhs = MakeUniformGrid(50.0f, 2.0, 0.5);

    std::unique_ptr<OESystem::OEScalarGrid> result(combine_maps(lhs, rhs, MapOp::MIN));
    ASSERT_NE(result, nullptr);

    const GridSummary s = Summarize(*result);
    ExpectPinned(s.sum, MaptitudePins::GRIDOPS_MIN_SUM);
    ExpectPinned(s.sum_sq, MaptitudePins::GRIDOPS_MIN_SUM_SQ);
    ExpectPinned(s.min, MaptitudePins::GRIDOPS_MIN_MIN);
    ExpectPinned(s.max, MaptitudePins::GRIDOPS_MIN_MAX);
    ExpectPinned(s.index_moment, MaptitudePins::GRIDOPS_MIN_INDEX_MOMENT);
}

TEST(GridOpsCharacterizationTest, CombineMax) {
    const OESystem::OEScalarGrid lhs = MakeRampGrid(2.0, 0.5);
    const OESystem::OEScalarGrid rhs = MakeUniformGrid(-50.0f, 2.0, 0.5);

    std::unique_ptr<OESystem::OEScalarGrid> result(combine_maps(lhs, rhs, MapOp::MAX));
    ASSERT_NE(result, nullptr);

    const GridSummary s = Summarize(*result);
    ExpectPinned(s.sum, MaptitudePins::GRIDOPS_MAX_SUM);
    ExpectPinned(s.sum_sq, MaptitudePins::GRIDOPS_MAX_SUM_SQ);
    ExpectPinned(s.min, MaptitudePins::GRIDOPS_MAX_MIN);
    ExpectPinned(s.max, MaptitudePins::GRIDOPS_MAX_MAX);
    ExpectPinned(s.index_moment, MaptitudePins::GRIDOPS_MAX_INDEX_MOMENT);
}

TEST(GridOpsCharacterizationTest, DiffToCalc) {
    const OESystem::OEScalarGrid obs = MakeRampGrid(2.0, 0.5);
    const OESystem::OEScalarGrid diff = MakeUniformGrid(3.0f, 2.0, 0.5);

    std::unique_ptr<OESystem::OEScalarGrid> result(diff_to_calc(obs, diff));
    ASSERT_NE(result, nullptr);

    const GridSummary s = Summarize(*result);
    ExpectPinned(s.sum, MaptitudePins::GRIDOPS_DIFF_TO_CALC_SUM);
    ExpectPinned(s.sum_sq, MaptitudePins::GRIDOPS_DIFF_TO_CALC_SUM_SQ);
    ExpectPinned(s.min, MaptitudePins::GRIDOPS_DIFF_TO_CALC_MIN);
    ExpectPinned(s.max, MaptitudePins::GRIDOPS_DIFF_TO_CALC_MAX);
    ExpectPinned(s.index_moment, MaptitudePins::GRIDOPS_DIFF_TO_CALC_INDEX_MOMENT);
}
