/// Tier 2 pins for the structure-factor pipeline.
///
/// src/DensityCalculator.cpp had no C++ coverage before Phase 1. These pins are
/// what make the FFTW RAII conversion verifiable. Tolerance is two-regime:
/// relative 1e-6 for |pinned| >= 1, absolute 1e-6 below it. The floor avoids
/// cross-machine flakiness when FFTW_ESTIMATE picks different codelets.
#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "maptitude/DensityCalculator.h"
#include "maptitude/SymOp.h"
#include "maptitude/UnitCell.h"

#include "fixtures.h"
#include "grid_summary.h"
#include "pin_values.h"

using namespace Maptitude;
using namespace MaptitudeTest;

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
    ExpectPinned(s.index_moment, MaptitudePins::FC_ORTHORHOMBIC_INDEX_MOMENT);
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
    ExpectPinned(s.index_moment, MaptitudePins::FC_MONOCLINIC_INDEX_MOMENT);
}

// Orthorhombic with n_scale_shells = 4, exercising the per-shell FFT scaling
// branch that Task 15 rewrites.
TEST(DensityCalculatorCharacterizationTest, OrthorhombicShells4) {
    UnitCell cell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0);
    std::vector<SymOp> symops = SymOp::ParseAll("x,y,z");
    OEChem::OEGraphMol mol = MakeAtomMol(6, 5.0, 5.0, 5.0);
    OESystem::OEScalarGrid obs = MakeGaussianGrid(5.0, 5.0, 5.0, 1.0, 6.0, 0.5);

    DensityCalculator calc(cell, symops);
    std::unique_ptr<OESystem::OEScalarGrid> fc(
        calc.Calculate(mol, obs, 2.0, nullptr, 0.35, 46.0, false, 4));
    ASSERT_NE(fc, nullptr);

    const GridSummary s = Summarize(*fc);
    ExpectPinned(s.sum, MaptitudePins::FC_ORTHORHOMBIC_SHELLS4_SUM);
    ExpectPinned(s.sum_sq, MaptitudePins::FC_ORTHORHOMBIC_SHELLS4_SUM_SQ);
    ExpectPinned(s.min, MaptitudePins::FC_ORTHORHOMBIC_SHELLS4_MIN);
    ExpectPinned(s.max, MaptitudePins::FC_ORTHORHOMBIC_SHELLS4_MAX);
    ExpectPinned(s.index_moment, MaptitudePins::FC_ORTHORHOMBIC_SHELLS4_INDEX_MOMENT);
}

// Asymmetric atom position breaks the permutation symmetry that makes the
// (5,5,5) cases above miss axis-order and stride regressions. The three
// distinct coordinates at (5.0, 2.0, -1.0) move the mean-centred index moment
// by ~10^5 times its tolerance under any axis permutation, catching the
// transpositions that the symmetric geometry cannot. Orthorhombic on purpose
// so Task 9 leaves it in place as the surviving permutation guard.
TEST(DensityCalculatorCharacterizationTest, OrthorhombicAsym) {
    UnitCell cell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0);
    std::vector<SymOp> symops = SymOp::ParseAll("x,y,z");
    OEChem::OEGraphMol mol = MakeAtomMol(6, 5.0, 2.0, -1.0);
    OESystem::OEScalarGrid obs = MakeGaussianGrid(5.0, 2.0, -1.0, 1.0, 6.0, 0.5);

    DensityCalculator calc(cell, symops);
    std::unique_ptr<OESystem::OEScalarGrid> fc(calc.Calculate(mol, obs, 2.0));
    ASSERT_NE(fc, nullptr);

    const GridSummary s = Summarize(*fc);
    ExpectPinned(s.sum, MaptitudePins::FC_ORTHORHOMBIC_ASYM_SUM);
    ExpectPinned(s.sum_sq, MaptitudePins::FC_ORTHORHOMBIC_ASYM_SUM_SQ);
    ExpectPinned(s.min, MaptitudePins::FC_ORTHORHOMBIC_ASYM_MIN);
    ExpectPinned(s.max, MaptitudePins::FC_ORTHORHOMBIC_ASYM_MAX);
    ExpectPinned(s.index_moment, MaptitudePins::FC_ORTHORHOMBIC_ASYM_INDEX_MOMENT);
}
