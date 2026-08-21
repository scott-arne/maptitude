/// Tier 2 pins for the structure-factor pipeline.
///
/// src/DensityCalculator.cpp had no C++ coverage before Phase 1. These pins are
/// what make the FFTW RAII conversion verifiable. Tolerance is relative 1e-6,
/// not 1e-12: FFTW_ESTIMATE may pick different codelets across builds.
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include "maptitude/DensityCalculator.h"
#include "maptitude/SymOp.h"
#include "maptitude/UnitCell.h"

#include "fixtures.h"
#include "pin_values.h"

using namespace Maptitude;
using namespace MaptitudeTest;

namespace {

constexpr double FC_RELATIVE_TOLERANCE = 1e-6;

struct GridSummary {
    double sum = 0.0;
    double sum_sq = 0.0;
    double min = 0.0;
    double max = 0.0;
};

GridSummary Summarize(const OESystem::OEScalarGrid& grid) {
    GridSummary s;
    s.min = grid[0];
    s.max = grid[0];
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        const double v = grid[i];
        s.sum += v;
        s.sum_sq += v * v;
        s.min = std::min(s.min, v);
        s.max = std::max(s.max, v);
    }
    return s;
}

void ExpectPinned(double actual, double pinned) {
    const double scale = std::max(1.0, std::abs(pinned));
    EXPECT_NEAR(actual, pinned, FC_RELATIVE_TOLERANCE * scale);
}

}  // namespace

TEST(DensityCalculatorCharacterizationTest, OrthorhombicP1) {
    UnitCell cell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0);
    std::vector<SymOp> symops = SymOp::ParseAll("x,y,z");
    OEChem::OEGraphMol mol = MakeAtomMol(6, 5.0, 5.0, 5.0);
    OESystem::OEScalarGrid obs = MakeGaussianGrid(5.0, 5.0, 5.0, 1.0, 6.0, 0.5);

    DensityCalculator calc(cell, symops);
    std::unique_ptr<OESystem::OEScalarGrid> fc(calc.Calculate(mol, obs, 2.0));
    ASSERT_NE(fc, nullptr);

    const GridSummary s = Summarize(*fc);
    ExpectPinned(s.sum, MaptitudePins::FC_ORTHORHOMBIC_SUM);
    ExpectPinned(s.sum_sq, MaptitudePins::FC_ORTHORHOMBIC_SUM_SQ);
    ExpectPinned(s.min, MaptitudePins::FC_ORTHORHOMBIC_MIN);
    ExpectPinned(s.max, MaptitudePins::FC_ORTHORHOMBIC_MAX);
}

// Task 9 replaces this test's body with an EXPECT_THROW(..., CellError). Until
// then it pins today's monoclinic output so the regression is deliberate.
TEST(DensityCalculatorCharacterizationTest, MonoclinicProducesNumbersToday) {
    UnitCell cell(20.0, 25.0, 30.0, 90.0, 105.0, 90.0);
    std::vector<SymOp> symops = SymOp::ParseAll("x,y,z");
    OEChem::OEGraphMol mol = MakeAtomMol(6, 5.0, 5.0, 5.0);
    OESystem::OEScalarGrid obs = MakeGaussianGrid(5.0, 5.0, 5.0, 1.0, 6.0, 0.5);

    DensityCalculator calc(cell, symops);
    std::unique_ptr<OESystem::OEScalarGrid> fc(calc.Calculate(mol, obs, 2.0));
    ASSERT_NE(fc, nullptr);

    const GridSummary s = Summarize(*fc);
    ExpectPinned(s.sum, MaptitudePins::FC_MONOCLINIC_SUM);
    ExpectPinned(s.sum_sq, MaptitudePins::FC_MONOCLINIC_SUM_SQ);
    ExpectPinned(s.min, MaptitudePins::FC_MONOCLINIC_MIN);
    ExpectPinned(s.max, MaptitudePins::FC_MONOCLINIC_MAX);
}
