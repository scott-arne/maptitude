/// Tier 2 pins for the structure-factor pipeline and behavioral guards for the
/// FFTW conversion.
///
/// src/DensityCalculator.cpp had no C++ coverage before Phase 1. These pins are
/// what make the FFTW RAII conversion verifiable. Tolerance is two-regime:
/// relative 1e-6 for |pinned| >= 1, absolute 1e-6 below it. The floor avoids
/// cross-machine flakiness when FFTW_ESTIMATE picks different codelets.
#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <thread>
#include <vector>

#include "maptitude/DensityCalculator.h"
#include "maptitude/Error.h"
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
    OESystem::OESkewGrid obs = MakeGaussianGrid(5.0, 5.0, 5.0, 1.0, 6.0, 0.5);

    DensityCalculator calc(cell, symops);
    std::unique_ptr<OESystem::OESkewGrid> fc(calc.Calculate(mol, obs, 2.0));
    ASSERT_NE(fc, nullptr);

    const GridSummary s = Summarize(*fc);
    ExpectPinned(s.sum, MaptitudePins::FC_ORTHORHOMBIC_SUM);
    ExpectPinned(s.sum_sq, MaptitudePins::FC_ORTHORHOMBIC_SUM_SQ);
    ExpectPinned(s.min, MaptitudePins::FC_ORTHORHOMBIC_MIN);
    ExpectPinned(s.max, MaptitudePins::FC_ORTHORHOMBIC_MAX);
    ExpectPinned(s.index_moment, MaptitudePins::FC_ORTHORHOMBIC_INDEX_MOMENT);
}

// Was a pin. Phase 1 narrows the supported domain: a monoclinic cell now
// raises rather than returning a value computed as though it were
// orthorhombic. Recorded in CHANGELOG.md as a capability regression.
TEST(DensityCalculatorCharacterizationTest, MonoclinicIsRejected) {
    std::vector<SymOp> symops = SymOp::ParseAll("x,y,z");
    EXPECT_THROW(DensityCalculator(UnitCell(20.0, 25.0, 30.0, 90.0, 105.0, 90.0), symops), CellError);
}

// Orthorhombic with n_scale_shells = 4, exercising the per-shell FFT scaling
// branch that Task 15 rewrites.
TEST(DensityCalculatorCharacterizationTest, OrthorhombicShells4) {
    UnitCell cell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0);
    std::vector<SymOp> symops = SymOp::ParseAll("x,y,z");
    OEChem::OEGraphMol mol = MakeAtomMol(6, 5.0, 5.0, 5.0);
    OESystem::OESkewGrid obs = MakeGaussianGrid(5.0, 5.0, 5.0, 1.0, 6.0, 0.5);

    DensityCalculator calc(cell, symops);
    std::unique_ptr<OESystem::OESkewGrid> fc(
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
    OESystem::OESkewGrid obs = MakeGaussianGrid(5.0, 2.0, -1.0, 1.0, 6.0, 0.5);

    DensityCalculator calc(cell, symops);
    std::unique_ptr<OESystem::OESkewGrid> fc(calc.Calculate(mol, obs, 2.0));
    ASSERT_NE(fc, nullptr);

    const GridSummary s = Summarize(*fc);
    ExpectPinned(s.sum, MaptitudePins::FC_ORTHORHOMBIC_ASYM_SUM);
    ExpectPinned(s.sum_sq, MaptitudePins::FC_ORTHORHOMBIC_ASYM_SUM_SQ);
    ExpectPinned(s.min, MaptitudePins::FC_ORTHORHOMBIC_ASYM_MIN);
    ExpectPinned(s.max, MaptitudePins::FC_ORTHORHOMBIC_ASYM_MAX);
    ExpectPinned(s.index_moment, MaptitudePins::FC_ORTHORHOMBIC_ASYM_INDEX_MOMENT);
}

// Calculate performs nine FFTW allocations and five plans, all held by the
// FftwBuffer and FftwPlan wrappers. This test does not exercise them: it throws
// from the argument check, before the first fftw_alloc_complex, so what it pins
// is that the guard fires repeatedly and leaves the pipeline usable. The
// allocation paths rest on the phase's one-time leak measurement instead --
// reverting the wrappers to raw allocation plus manual fftw_free on the success
// path leaves this test, and the whole suite, green. See the corresponding
// Known limitations entry in CHANGELOG.md.
TEST(DensityCalculatorCharacterizationTest, ThrowingPathDoesNotDestabilizeTheProcess) {
    UnitCell cell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0);
    std::vector<SymOp> symops = SymOp::ParseAll("x,y,z");
    OEChem::OEGraphMol mol = MakeAtomMol(6, 5.0, 5.0, 5.0);
    OESystem::OESkewGrid obs = MakeGaussianGrid(5.0, 5.0, 5.0, 1.0, 6.0, 0.5);

    DensityCalculator calc(cell, symops);
    for (int i = 0; i < 50; ++i) {
        EXPECT_THROW(calc.Calculate(mol, obs, -1.0), GridError);
    }

    // The pipeline must still work afterwards.
    std::unique_ptr<OESystem::OESkewGrid> fc(calc.Calculate(mol, obs, 2.0));
    ASSERT_NE(fc, nullptr);
    ExpectPinned(Summarize(*fc).sum, MaptitudePins::FC_ORTHORHOMBIC_SUM);
}

// fftw_plan_dft_3d and fftw_destroy_plan mutate global planner state and are
// not thread-safe. fftw_execute on an already-created plan is thread-safe. A
// C++ caller invoking Calculate from multiple threads reaches this directly;
// Python callers are currently serialized by the GIL since the module is not
// built with SWIG threading. The planner mutex guards only plan creation and
// destruction.
TEST(DensityCalculatorCharacterizationTest, ConcurrentCalculateIsSafe) {
    UnitCell cell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0);
    std::vector<SymOp> symops = SymOp::ParseAll("x,y,z");

    std::vector<std::thread> threads;
    std::vector<double> sums(4, 0.0);
    for (int t = 0; t < 4; ++t) {
        threads.emplace_back([&, t]() {
            OEChem::OEGraphMol mol = MakeAtomMol(6, 5.0, 5.0, 5.0);
            OESystem::OESkewGrid obs = MakeGaussianGrid(5.0, 5.0, 5.0, 1.0, 6.0, 0.5);
            DensityCalculator calc(cell, symops);
            std::unique_ptr<OESystem::OESkewGrid> fc(calc.Calculate(mol, obs, 2.0));
            sums[t] = Summarize(*fc).sum;
        });
    }
    for (auto& thread : threads) {
        thread.join();
    }

    for (double sum : sums) {
        ExpectPinned(sum, MaptitudePins::FC_ORTHORHOMBIC_SUM);
    }
}

namespace {

/// Build a 21x21x21 grid spanning [0, 10] on each axis, with no unit cell of
/// its own. Deliberately not a `fixtures.h` helper: those all call
/// `SetUnitCell`, which is the property these two tests need absent.
OESystem::OESkewGrid MakeCelllessGrid(float value) {
    OESystem::OESkewGrid grid;
    // The setters report failure by return value. EXPECT_* would record a
    // failure and carry on into the loop below, which would then write through
    // whatever GetValues() yields for a grid whose geometry was never
    // established. ASSERT_* cannot be used instead: it expands to a bare return
    // and this function returns a grid. So throw, as MakeEmptyGrid in
    // fixtures.h does for the same reason.
    if (!grid.SetDim(21u, 21u, 21u) || !grid.SetSpacing(0.5f) ||
        !grid.SetMid(5.0f, 5.0f, 5.0f) || grid.GetValues() == nullptr) {
        throw std::invalid_argument("MakeCelllessGrid: the skew carrier rejected the geometry");
    }
    float* values = grid.GetValues();
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        values[i] = value;
    }
    return grid;
}

}  // namespace

TEST(DensityCalculatorSampling, SucceedsOnAGridWithNoUnitCell) {
    // The FFT counts come from the constructor's cell and the map's node
    // intervals, not from the grid's own cell, so a grid carrying no cell is
    // not an error. The 0.5 A interval divides all three edges exactly
    // (20/0.5 = 40, 25/0.5 = 50, 30/0.5 = 60), so the divisibility guard
    // Step 7(b) added passes.
    OESystem::OESkewGrid obs = MakeCelllessGrid(1.0f);
    ASSERT_FALSE(obs.HasUnitCell());

    UnitCell cell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0);
    std::vector<SymOp> symops = SymOp::ParseAll("x,y,z");
    OEChem::OEGraphMol mol = MakeAtomMol(6, 5.0, 5.0, 5.0);

    DensityCalculator calc(cell, symops);
    std::unique_ptr<OESystem::OESkewGrid> fc(calc.Calculate(mol, obs, 2.0));
    ASSERT_NE(fc, nullptr);

    // Non-null is not enough: an empty or zeroed allocation would pass it. The
    // result is written into a grid built from the observed grid's geometry, so
    // it carries the input's dimensions rather than the 40x50x60 FFT counts.
    EXPECT_EQ(fc->GetXDim(), 21u);
    EXPECT_EQ(fc->GetYDim(), 21u);
    EXPECT_EQ(fc->GetZDim(), 21u);
    ASSERT_EQ(fc->GetSize(), obs.GetSize());

    const float* fc_values = fc->GetValues();
    ASSERT_NE(fc_values, nullptr);
    bool any_nonzero = false;
    for (unsigned int i = 0; i < fc->GetSize(); ++i) {
        ASSERT_TRUE(std::isfinite(fc_values[i])) << "non-finite density at element " << i;
        if (fc_values[i] != 0.0f) {
            any_nonzero = true;
        }
    }
    EXPECT_TRUE(any_nonzero) << "the returned grid is entirely zero";
}

TEST(DensityCalculatorSampling, UsesTheConstructorCellNotTheGridCell) {
    // Two grids with identical node geometry, one carrying a unit cell that
    // disagrees with the calculator's (10.5 A cubic against 20x25x30), must
    // produce identical output -- the grid's own cell is never read.
    //
    // Measured: layering SetUnitCell(10.5, 10.5, 10.5, 90, 90, 90, 21, 21, 21)
    // onto this grid leaves every ElementToSpatialCoord result bit-identical,
    // so the only difference between the two inputs is HasUnitCell().
    OESystem::OESkewGrid without_cell = MakeCelllessGrid(1.0f);
    OESystem::OESkewGrid with_wrong_cell = MakeCelllessGrid(1.0f);
    ASSERT_TRUE(with_wrong_cell.SetUnitCell(10.5f, 10.5f, 10.5f,
                                            90.0f, 90.0f, 90.0f,
                                            21u, 21u, 21u));
    ASSERT_FALSE(without_cell.HasUnitCell());
    ASSERT_TRUE(with_wrong_cell.HasUnitCell());

    UnitCell cell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0);
    std::vector<SymOp> symops = SymOp::ParseAll("x,y,z");
    OEChem::OEGraphMol mol = MakeAtomMol(6, 5.0, 5.0, 5.0);

    DensityCalculator calc(cell, symops);
    std::unique_ptr<OESystem::OESkewGrid> a(calc.Calculate(mol, without_cell, 2.0));
    std::unique_ptr<OESystem::OESkewGrid> b(calc.Calculate(mol, with_wrong_cell, 2.0));
    ASSERT_NE(a, nullptr);
    ASSERT_NE(b, nullptr);

    const GridSummary sa = Summarize(*a);
    const GridSummary sb = Summarize(*b);
    EXPECT_DOUBLE_EQ(sa.sum, sb.sum);
    EXPECT_DOUBLE_EQ(sa.sum_sq, sb.sum_sq);
    EXPECT_DOUBLE_EQ(sa.min, sb.min);
    EXPECT_DOUBLE_EQ(sa.max, sb.max);
    EXPECT_DOUBLE_EQ(sa.index_moment, sb.index_moment);
}
