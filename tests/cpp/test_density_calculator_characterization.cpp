/// Tier 2 pins for the structure-factor pipeline and behavioral guards for the
/// FFTW conversion.
///
/// src/DensityCalculator.cpp had no C++ coverage before Phase 1. These pins are
/// what make the FFTW RAII conversion verifiable. Tolerance is two-regime:
/// relative 1e-6 for |pinned| >= 1, absolute 1e-6 below it. The floor avoids
/// cross-machine flakiness when FFTW_ESTIMATE picks different codelets.
#include <gtest/gtest.h>

#include <memory>
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

// Calculate performs nine FFTW allocations and five plans. The RAII conversion
// guards against leaks when std::bad_alloc or allocation failure escapes from
// the pipeline. The null checks added during the conversion create new throw
// sites in the middle of the allocation region, which is what makes the RAII
// conversion load-bearing rather than precautionary.
TEST(DensityCalculatorCharacterizationTest, ThrowingPathDoesNotDestabilizeTheProcess) {
    UnitCell cell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0);
    std::vector<SymOp> symops = SymOp::ParseAll("x,y,z");
    OEChem::OEGraphMol mol = MakeAtomMol(6, 5.0, 5.0, 5.0);
    OESystem::OEScalarGrid obs = MakeGaussianGrid(5.0, 5.0, 5.0, 1.0, 6.0, 0.5);

    DensityCalculator calc(cell, symops);
    for (int i = 0; i < 50; ++i) {
        // A non-positive resolution is rejected early (before any FFTW
        // allocation), so this test proves the guard works but exercises no
        // FFTW cleanup. Rely on Step 6's leak measurement for the allocation
        // paths.
        EXPECT_THROW(calc.Calculate(mol, obs, -1.0), GridError);
    }

    // The pipeline must still work afterwards.
    std::unique_ptr<OESystem::OEScalarGrid> fc(calc.Calculate(mol, obs, 2.0));
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
            OESystem::OEScalarGrid obs = MakeGaussianGrid(5.0, 5.0, 5.0, 1.0, 6.0, 0.5);
            DensityCalculator calc(cell, symops);
            std::unique_ptr<OESystem::OEScalarGrid> fc(calc.Calculate(mol, obs, 2.0));
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
