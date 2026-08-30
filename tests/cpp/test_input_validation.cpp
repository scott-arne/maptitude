/// Rejection tests for every validation throw site added in Phase 1.
///
/// Named test_input_validation to avoid confusion with the unrelated
/// tests/python/test_validation.py, which validates scores against reference
/// data.
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <oegrid.h>

#include "maptitude/CoverageOptions.h"
#include "maptitude/DensityCalculator.h"
#include "maptitude/Error.h"
#include "maptitude/Grid.h"
#include "maptitude/GridOps.h"
#include "maptitude/Metric.h"
#include "maptitude/QScoreOptions.h"
#include "maptitude/SymOp.h"
#include "maptitude/UnitCell.h"

#include "fixtures.h"

using namespace Maptitude;
using MaptitudeTest::MakeAtomMol;
using MaptitudeTest::MakeEmptyGrid;
using MaptitudeTest::MakeGaussianGrid;

// ---- Acceptance tests: valid cells must pass ----

TEST(CellValidationTest, AcceptsWellFormedCells) {
    // Cubic
    EXPECT_NO_THROW(UnitCell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0));
    // Orthorhombic
    EXPECT_NO_THROW(UnitCell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0));
    // Monoclinic
    EXPECT_NO_THROW(UnitCell(20.0, 25.0, 30.0, 90.0, 105.0, 90.0));
    // Hexagonal
    EXPECT_NO_THROW(UnitCell(10.0, 10.0, 15.0, 90.0, 90.0, 120.0));
    // Rhombohedral
    EXPECT_NO_THROW(UnitCell(10.0, 10.0, 10.0, 80.0, 80.0, 80.0));
    // Triclinic
    EXPECT_NO_THROW(UnitCell(10.0, 10.0, 10.0, 60.0, 70.0, 80.0));
}

TEST(CellValidationTest, AcceptsNearDegenerateTriclinic) {
    // The most degenerate realistic triclinic tested (angles 30/30/40 degrees,
    // radicand 0.062) must not be refused by the MIN_RADICAND floor.
    EXPECT_NO_THROW(UnitCell(30.0, 30.0, 40.0, 30.0, 30.0, 40.0));
}

// ---- Rejection tests via mutation: exercise the free function ----
// The parameterized constructor itself validates, so calling
// validate_cell(UnitCell(...)) never enters the free function -- the
// constructor throws while the argument is being evaluated. These tests
// construct a valid cell, mutate it to a bad state, then assert validate_cell
// throws. This proves the free function itself works, which is load-bearing
// for DensityCalculator's second check and for the geometry readers.

TEST(CellValidationTest, RejectsNonPositiveLengthsViaMutation) {
    UnitCell cell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0);
    cell.a = 0.0;
    EXPECT_THROW(validate_cell(cell), CellError);

    cell = UnitCell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0);
    cell.b = -1.0;
    EXPECT_THROW(validate_cell(cell), CellError);

    cell = UnitCell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0);
    cell.c = 0.0;
    EXPECT_THROW(validate_cell(cell), CellError);
}

TEST(CellValidationTest, RejectsNonFiniteLengthsViaMutation) {
    const double nan_value = std::numeric_limits<double>::quiet_NaN();
    const double inf_value = std::numeric_limits<double>::infinity();

    UnitCell cell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0);
    cell.a = nan_value;
    EXPECT_THROW(validate_cell(cell), CellError);

    cell = UnitCell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0);
    cell.b = inf_value;
    EXPECT_THROW(validate_cell(cell), CellError);
}

TEST(CellValidationTest, RejectsAnglesOutsideTheOpenIntervalViaMutation) {
    UnitCell cell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0);
    cell.alpha = 0.0;
    EXPECT_THROW(validate_cell(cell), CellError);

    cell = UnitCell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0);
    cell.alpha = 180.0;
    EXPECT_THROW(validate_cell(cell), CellError);

    cell = UnitCell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0);
    cell.beta = -10.0;
    EXPECT_THROW(validate_cell(cell), CellError);

    // cos(200) == cos(160), so the radicand test alone would let this through.
    cell = UnitCell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0);
    cell.gamma = 200.0;
    EXPECT_THROW(validate_cell(cell), CellError);
}

TEST(CellValidationTest, RejectsAngleTriplesWithNoRealLatticeViaMutation) {
    // Each angle is individually in range, but no cell has this geometry.
    UnitCell cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
    cell.alpha = cell.beta = cell.gamma = 150.0;
    EXPECT_THROW(validate_cell(cell), CellError);

    cell = UnitCell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
    cell.alpha = cell.beta = 20.0;
    cell.gamma = 150.0;
    EXPECT_THROW(validate_cell(cell), CellError);
}

TEST(CellValidationTest, RejectsNearDegenerateImpossibleTriples) {
    // F1: impossible cells whose radicand rounds to a small positive value due
    // to floating-point noise. These were accepted before the MIN_RADICAND floor.
    UnitCell cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
    cell.alpha = cell.beta = 0.5;
    cell.gamma = 1.0000000010000001;
    EXPECT_THROW(validate_cell(cell), CellError);

    cell = UnitCell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
    cell.alpha = cell.beta = 0.001;
    cell.gamma = 0.002000001;
    EXPECT_THROW(validate_cell(cell), CellError);
}

TEST(CellValidationTest, RejectsUnderflowCell) {
    // F2: finite positive lengths that underflow a*b*c to zero, producing NaN
    // from the transforms. The volume check catches this.
    UnitCell cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
    cell.a = cell.b = cell.c = 1e-200;
    EXPECT_THROW(validate_cell(cell), CellError);
}

TEST(CellValidationTest, RejectsOverflowCell) {
    // F2: lengths so large that a*b*c overflows to infinity.
    UnitCell cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
    cell.a = cell.b = cell.c = 1e200;
    EXPECT_THROW(validate_cell(cell), CellError);
}

// ---- Constructor tests ----

TEST(CellValidationTest, ParameterizedConstructorRejectsBadInputImmediately) {
    // The parameterized constructor calls validate_cell, so bad input is
    // rejected at construction time.
    EXPECT_THROW(UnitCell(0.0, 25.0, 30.0, 90.0, 90.0, 90.0), CellError);
    EXPECT_THROW(UnitCell(10.0, 10.0, 10.0, 150.0, 150.0, 150.0), CellError);
    EXPECT_THROW(UnitCell(1e-200, 1e-200, 1e-200, 90.0, 90.0, 90.0), CellError);
}

TEST(CellValidationTest, DefaultConstructionStillWorks) {
    // A default cell is all zeros and would fail validate_cell. That is why
    // validation lives in a free function called at consumption points rather
    // than unconditionally in every constructor.
    EXPECT_NO_THROW(UnitCell{});
}

// ---- Reader validation (F3) ----

TEST(CellValidationTest, GeometryReadersValidateBeforeComputing) {
    // F3: Volume(), OrthogonalizationMatrix(), and DeorthogonalizationMatrix()
    // are consumption points and must validate. A cell can be invalidated after
    // construction (members are public), so readers must check.
    UnitCell cell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0);
    cell.a = 0.0;

    EXPECT_THROW(cell.Volume(), CellError);
    EXPECT_THROW(cell.OrthogonalizationMatrix(), CellError);
    EXPECT_THROW(cell.DeorthogonalizationMatrix(), CellError);

    // CartesianToFractional delegates to DeorthogonalizationMatrix, so it is
    // covered transitively.
    EXPECT_THROW(cell.CartesianToFractional(1.0, 2.0, 3.0), CellError);
}

TEST(CellValidationTest, DefaultConstructedCellThrowsFromReaders) {
    // A default-constructed cell has a=b=c=0, which is invalid. The readers
    // must reject it.
    UnitCell cell{};
    EXPECT_THROW(cell.Volume(), CellError);
    EXPECT_THROW(cell.OrthogonalizationMatrix(), CellError);
    EXPECT_THROW(cell.DeorthogonalizationMatrix(), CellError);
}

// ---- Matrix denominator overflow (F2 round 2) ----

TEST(CellValidationTest, RejectsOrthogonalizationMatrixOverflow) {
    // Fix round 2: a*b underflows to zero but large c rescues the volume,
    // leaving vol/(a*b*sg) = inf in the orthogonalization matrix.
    UnitCell cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
    cell.a = cell.b = 3.16228e-162;
    cell.c = 1e280;
    cell.gamma = 0.002;
    EXPECT_THROW(validate_cell(cell), CellError);
}

TEST(CellValidationTest, RejectsDeorthogonalizationMatrixOverflow) {
    // Fix round 2: subnormal a with large b*c causes b*c/vol to overflow
    // in the deorthogonalization matrix even though vol is finite.
    UnitCell cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
    cell.a = 1e-310;
    cell.b = cell.c = 1e155;
    EXPECT_THROW(validate_cell(cell), CellError);
}

// ---- Condition-number bound for coordinate conversion (F5 round 3/5) ----

TEST(CellValidationTest, RejectsWorstRoundTripCorruption) {
    // Extreme length ratios can produce finite matrix entries yet return
    // catastrophically wrong coordinates. This cell has condition number
    // ~1e308 (one row sum overflows to inf), far exceeding the threshold.
    UnitCell cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
    cell.a = 1e-20;
    cell.b = 1e-300;
    cell.c = 1e20;
    cell.beta = 105.0;
    EXPECT_THROW(validate_cell(cell), CellError);
}

TEST(CellValidationTest, RejectsCaseThatIdentityCheckMisses) {
    // The proposed deortho*ortho identity check cannot distinguish this cell
    // from a valid hexagonal cell (both have identity error 1.11e-16). The
    // condition-number bound catches it (kappa ~1e308, row sum overflows).
    UnitCell cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
    cell.a = 1e-200;
    cell.b = cell.c = 1e120;
    EXPECT_THROW(validate_cell(cell), CellError);
}

TEST(CellValidationTest, RejectsExtremeLengthRatioFlaggedByReview) {
    // Interior point (0.25, 0.5, 0.75) round-trips to (-1.66e+276, 0.5, 0.75)
    // while basis vectors round-trip cleanly, which is why deortho*ortho
    // identity check and earlier measurements missed this. Condition number
    // bound catches all cases regardless of which points are probed.
    EXPECT_THROW(UnitCell(1e-308, 1.0, 1.0, 90.0, 90.0, 90.0), CellError);
}

TEST(CellValidationTest, RejectsCellWithExactProbeButWrongInterior) {
    // Fix round 5: short-mantissa probes are structurally blind to catastrophic
    // cancellation. This cell's (0.25, 0.5, 0.75) round trip is *exact* while
    // (0.172..., 0.985..., 0.604...) comes back as (0.125, ...), error 0.047.
    // The probe was chosen as exact binary fractions for reproducibility, which
    // is exactly why it cannot see low-order bit destruction. Condition number
    // is 5.74e226, far exceeding threshold.
    UnitCell cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
    cell.a = 7.241769560140648e110;
    cell.b = 7.532447449398725e125;
    cell.c = 1.9046765643114865e-101;
    cell.alpha = 63.52960342238751;
    cell.beta = 43.63896561425281;
    cell.gamma = 71.17433010912423;
    EXPECT_THROW(validate_cell(cell), CellError);
}

TEST(CellValidationTest, RejectsWorstProbeEscape) {
    // Fix round 5: worst escape from the round-trip probe in a 120k-cell sweep.
    // Probe accepted it while another interior point came back with absolute
    // error 1.227. Condition number 1.77e201 exceeds threshold.
    UnitCell cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
    cell.a = 4.73235506481289e135;
    cell.b = 5.588102857699266e-66;
    cell.c = 2.913459289783779e-50;
    cell.alpha = 136.82442378774567;
    cell.beta = 98.94351481288861;
    cell.gamma = 90.04905757433274;
    EXPECT_THROW(validate_cell(cell), CellError);
}

TEST(CellValidationTest, AcceptsLegitimateExtremes) {
    // Fix round 5: pin the condition-number threshold from below. A thin plate
    // (kappa 500) and a ribosome cell (kappa 3.79) are legitimate despite being
    // at the extreme end of the plausible range. Neither may be refused.
    EXPECT_NO_THROW(UnitCell(2.0, 2.0, 1000.0, 90.0, 90.0, 90.0));
    EXPECT_NO_THROW(UnitCell(500.0, 500.0, 1200.0, 90.0, 90.0, 120.0));
}

TEST(CellValidationTest, AcceptsPhysicalRangeEnds) {
    // 1 Angstrom and 3000 Angstrom cubes span the physically plausible range
    // (small molecules to large virus assemblies). Neither may be refused.
    EXPECT_NO_THROW(UnitCell(1.0, 1.0, 1.0, 90.0, 90.0, 90.0));
    EXPECT_NO_THROW(UnitCell(3000.0, 3000.0, 3000.0, 90.0, 90.0, 90.0));
}

// ---- Inverse-residual bound for coordinate conversion (F5 round 6) ----

TEST(CellValidationTest, RejectsSubnormalVolumeThatBreaksTheInverse) {
    // The condition number is only meaningful when the deorthogonalization
    // matrix really is the inverse of the orthogonalization matrix. Here the
    // volume is subnormal at 8.1e-322 (about seven bits of precision), and
    // every deorthogonalization entry divides by it, so the "inverse" inherits
    // the volume's ~0.6% relative error. kappa reads a healthy 2.37 because it
    // is computed from a matrix that is not the inverse, while deortho*ortho
    // sits 1.38e-3 off the identity and (0.25, 0.5, 0.75) round-trips with
    // error 1.03e-3. Only the residual check sees this.
    EXPECT_THROW(UnitCell(1e-107, 1e-107, 1e-107, 60.0, 70.0, 80.0), CellError);
}

TEST(CellValidationTest, AcceptsUsableCellWithSubnormalVolume) {
    // Pin the residual threshold from below, and guard against anyone later
    // "simplifying" the fix into a subnormal-volume rejection. This cell's
    // volume is 8.138e-316 — also deeply subnormal — yet it round-trips to
    // 5e-11 and is perfectly usable. The check gates on usability, not on
    // subnormality.
    EXPECT_NO_THROW(UnitCell(1e-105, 1e-105, 1e-105, 60.0, 70.0, 80.0));
}

// ---- Orthorhombic lattice constraint (Task 9) ----

TEST(OrthorhombicGuardTest, AcceptsOrthorhombicCells) {
    std::vector<SymOp> symops = SymOp::ParseAll("x,y,z");
    EXPECT_NO_THROW(DensityCalculator(UnitCell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0), symops));
    EXPECT_NO_THROW(DensityCalculator(UnitCell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0), symops));
}

TEST(OrthorhombicGuardTest, RejectsMonoclinicAndTriclinicCells) {
    std::vector<SymOp> symops = SymOp::ParseAll("x,y,z");
    EXPECT_THROW(DensityCalculator(UnitCell(20.0, 25.0, 30.0, 90.0, 105.0, 90.0), symops), CellError);
    EXPECT_THROW(DensityCalculator(UnitCell(10.0, 10.0, 15.0, 90.0, 90.0, 120.0), symops), CellError);
    EXPECT_THROW(DensityCalculator(UnitCell(10.0, 12.0, 14.0, 88.0, 95.0, 101.0), symops), CellError);
}

TEST(OrthorhombicGuardTest, ToleranceIsOnTheCosineNotTheAngle) {
    std::vector<SymOp> symops = SymOp::ParseAll("x,y,z");
    // cos(90 deg) is ~6.1e-17 in IEEE-754: well inside the tolerance.
    EXPECT_NO_THROW(DensityCalculator(UnitCell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0), symops));
    // cos(90.00000001 deg) is ~1.7e-10: inside the cosine tolerance, would fail
    // a degree tolerance. This case distinguishes the two rules.
    EXPECT_NO_THROW(DensityCalculator(UnitCell(20.0, 25.0, 30.0, 90.00000001, 90.0, 90.0), symops));
    // cos(90.0000001 deg) is ~1.7e-9: outside the cosine tolerance (note: one
    // fewer decimal place in the angle, 10x larger cosine).
    EXPECT_THROW(DensityCalculator(UnitCell(20.0, 25.0, 30.0, 90.0000001, 90.0, 90.0), symops),
                 CellError);
}

TEST(OrthorhombicGuardTest, CellValidityIsCheckedBeforeTheLatticeType) {
    // A cell that is both invalid and non-orthorhombic should report the
    // geometric failure, which is the more specific diagnosis. UnitCell's
    // members are public and mutable, so we construct a valid cell then mutate
    // it to reach DensityCalculator's constructor (passing angles=(150,150,150)
    // directly to UnitCell(...) throws before DensityCalculator is entered).
    std::vector<SymOp> symops = SymOp::ParseAll("x,y,z");
    UnitCell cell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0);
    cell.alpha = 150.0;
    cell.beta = 150.0;
    cell.gamma = 150.0;
    try {
        DensityCalculator(cell, symops);
        FAIL() << "expected CellError";
    } catch (const CellError& e) {
        EXPECT_NE(std::string(e.what()).find("geometrically impossible"), std::string::npos)
            << "got: " << e.what();
    }
}

// ---- QScoreOptions validation ----

TEST(QScoreOptionsValidationTest, AcceptsTheDefaults) {
    QScoreOptions options;
    EXPECT_GT(options.GetSigma(), 0.0);
    EXPECT_GT(options.GetRadialStep(), 0.0);
    EXPECT_GT(options.GetMaxRadius(), 0.0);
    EXPECT_GE(options.GetNumPoints(), 1u);
}

TEST(QScoreOptionsValidationTest, RejectsNonPositiveSigma) {
    QScoreOptions options;
    const double default_value = options.GetSigma();
    EXPECT_THROW(options.SetSigma(0.0), std::invalid_argument);
    EXPECT_DOUBLE_EQ(options.GetSigma(), default_value);  // Field unchanged after rejection
    EXPECT_THROW(options.SetSigma(-0.5), std::invalid_argument);
    EXPECT_DOUBLE_EQ(options.GetSigma(), default_value);  // Still unchanged
    EXPECT_NO_THROW(options.SetSigma(0.8));
    EXPECT_DOUBLE_EQ(options.GetSigma(), 0.8);
}

TEST(QScoreOptionsValidationTest, RejectsNonPositiveRadialStep) {
    // A non-positive step makes the radial sweep never advance: the process
    // hangs rather than returning a wrong number.
    QScoreOptions options;
    const double default_value = options.GetRadialStep();
    EXPECT_THROW(options.SetRadialStep(0.0), std::invalid_argument);
    EXPECT_DOUBLE_EQ(options.GetRadialStep(), default_value);  // Field unchanged after rejection
    EXPECT_THROW(options.SetRadialStep(-0.1), std::invalid_argument);
    EXPECT_DOUBLE_EQ(options.GetRadialStep(), default_value);  // Still unchanged
    EXPECT_NO_THROW(options.SetRadialStep(0.25));
    EXPECT_DOUBLE_EQ(options.GetRadialStep(), 0.25);
}

TEST(QScoreOptionsValidationTest, RejectsNonPositiveMaxRadius) {
    QScoreOptions options;
    const double default_value = options.GetMaxRadius();
    EXPECT_THROW(options.SetMaxRadius(0.0), std::invalid_argument);
    EXPECT_DOUBLE_EQ(options.GetMaxRadius(), default_value);  // Field unchanged after rejection
    EXPECT_THROW(options.SetMaxRadius(-2.0), std::invalid_argument);
    EXPECT_DOUBLE_EQ(options.GetMaxRadius(), default_value);  // Still unchanged
    EXPECT_NO_THROW(options.SetMaxRadius(3.0));
    EXPECT_DOUBLE_EQ(options.GetMaxRadius(), 3.0);
}

TEST(QScoreOptionsValidationTest, RejectsZeroSamplePoints) {
    QScoreOptions options;
    const unsigned int default_value = options.GetNumPoints();
    EXPECT_THROW(options.SetNumPoints(0), std::invalid_argument);
    EXPECT_EQ(options.GetNumPoints(), default_value);  // Field unchanged after rejection
    EXPECT_NO_THROW(options.SetNumPoints(1));
    EXPECT_EQ(options.GetNumPoints(), 1u);
    EXPECT_NO_THROW(options.SetNumPoints(16));
    EXPECT_EQ(options.GetNumPoints(), 16u);
}

TEST(QScoreOptionsValidationTest, RejectsNonFiniteValues) {
    const double nan_value = std::numeric_limits<double>::quiet_NaN();
    const double pos_inf = std::numeric_limits<double>::infinity();
    const double neg_inf = -std::numeric_limits<double>::infinity();
    QScoreOptions options;
    EXPECT_THROW(options.SetSigma(nan_value), std::invalid_argument);
    EXPECT_THROW(options.SetSigma(pos_inf), std::invalid_argument);
    EXPECT_THROW(options.SetSigma(neg_inf), std::invalid_argument);
    EXPECT_THROW(options.SetRadialStep(nan_value), std::invalid_argument);
    EXPECT_THROW(options.SetRadialStep(pos_inf), std::invalid_argument);
    EXPECT_THROW(options.SetRadialStep(neg_inf), std::invalid_argument);
    EXPECT_THROW(options.SetMaxRadius(nan_value), std::invalid_argument);
    EXPECT_THROW(options.SetMaxRadius(pos_inf), std::invalid_argument);
    EXPECT_THROW(options.SetMaxRadius(neg_inf), std::invalid_argument);
}

TEST(QScoreOptionsValidationTest, RejectsARadialStepBelowTheShellKeyResolution) {
    // Shell keys are round(R * 1e6), so steps below 1e-6 alias distinct shells onto
    // one key; subnormal steps cannot advance the accumulator at all.
    QScoreOptions options;
    EXPECT_THROW(options.SetRadialStep(std::numeric_limits<double>::denorm_min()),
                 std::invalid_argument);
    EXPECT_THROW(options.SetRadialStep(1e-9), std::invalid_argument);
    EXPECT_NO_THROW(options.SetRadialStep(MIN_RADIAL_STEP));
    EXPECT_NO_THROW(options.SetRadialStep(0.25));
}

TEST(QScoreOptionsValidationTest, RejectsAMaxRadiusThatOverflowsTheShellKey) {
    // (2147.47 + 0.01) * 1e6 is the largest shell key that fits in int.
    QScoreOptions options;
    EXPECT_THROW(options.SetMaxRadius(2148.0), std::invalid_argument);
    EXPECT_THROW(options.SetMaxRadius(std::numeric_limits<double>::max()),
                 std::invalid_argument);
    EXPECT_NO_THROW(options.SetMaxRadius(MAX_RADIUS_LIMIT));
    EXPECT_NO_THROW(options.SetMaxRadius(3.0));
}

TEST(QScoreOptionsValidationTest, RejectsAnAbsurdSamplePointCount) {
    QScoreOptions options;
    EXPECT_THROW(options.SetNumPoints(std::numeric_limits<unsigned int>::max()),
                 std::invalid_argument);
    EXPECT_THROW(options.SetNumPoints(MAX_NUM_POINTS + 1), std::invalid_argument);
    EXPECT_NO_THROW(options.SetNumPoints(MAX_NUM_POINTS));
    EXPECT_NO_THROW(options.SetNumPoints(8));
}

TEST(QScoreOptionsValidationTest, RejectsNegativeZeroAndNegativeSubnormals) {
    // -0.0 compares <= 0.0, so it is already refused; pin that so a future rewrite
    // of the predicate cannot quietly start accepting it.
    QScoreOptions options;
    EXPECT_THROW(options.SetSigma(-0.0), std::invalid_argument);
    EXPECT_THROW(options.SetSigma(-std::numeric_limits<double>::denorm_min()),
                 std::invalid_argument);
    EXPECT_THROW(options.SetMaxRadius(-0.0), std::invalid_argument);
}

// ---- Enum-valued setters reject values their enum does not declare ----
//
// A scoped enum with underlying type `int` can hold any int. `static_cast<AtomRadius>(42)`
// and `static_cast<RadialSampling>(42)` are valid values of their types, and SWIG hands
// one straight through from `SetAtomRadiusMethod(42)` in Python without a cast. Every
// switch over these enums is written without a `default:` label, which is what makes a
// future enumerator a -Wswitch warning; that property says nothing about values the enum
// does not declare, so the setters are where those are closed.

TEST(EnumSetterValidationTest, RsccRejectsAnUndeclaredAtomRadius) {
    // Before this guard the undeclared value matched no arm of the radius switch in
    // `rscc`, leaving `double radius;` indeterminate where get_atom_grid_points consumes
    // it. Measured on a Gaussian grid against a matching calc grid, `overall` came back
    // as a confident 1.0 -- a perfect model-to-map fit -- beside `by_atom` entries that
    // were all NaN.
    RsccOptions options;
    const AtomRadius default_value = options.GetAtomRadiusMethod();
    for (int undeclared : {42, 999, -1}) {
        EXPECT_THROW(options.SetAtomRadiusMethod(static_cast<AtomRadius>(undeclared)),
                     std::invalid_argument)
            << "accepted " << undeclared;
        EXPECT_EQ(options.GetAtomRadiusMethod(), default_value)
            << "field moved after rejecting " << undeclared;
    }
}

TEST(EnumSetterValidationTest, RsrRejectsAnUndeclaredAtomRadius) {
    // RsrOptions carries its own copy of the setter, so it needs its own case: the two
    // classes share the enum but not the code.
    RsrOptions options;
    const AtomRadius default_value = options.GetAtomRadiusMethod();
    for (int undeclared : {42, 999, -1}) {
        EXPECT_THROW(options.SetAtomRadiusMethod(static_cast<AtomRadius>(undeclared)),
                     std::invalid_argument)
            << "accepted " << undeclared;
        EXPECT_EQ(options.GetAtomRadiusMethod(), default_value)
            << "field moved after rejecting " << undeclared;
    }
}

TEST(EnumSetterValidationTest, StillAcceptsEveryDeclaredAtomRadius) {
    // The accepting half. All four enumerators must round-trip through both setters,
    // including ADAPTIVE on RsccOptions: rscc rejects that model when it scores, not when
    // it is configured, and moving the rejection into the setter would be a second
    // behavior change rather than the hole this closes.
    for (AtomRadius declared : {AtomRadius::FIXED, AtomRadius::SCALED, AtomRadius::BINNED,
                                AtomRadius::ADAPTIVE}) {
        RsccOptions rscc_options;
        EXPECT_NO_THROW(rscc_options.SetAtomRadiusMethod(declared));
        EXPECT_EQ(rscc_options.GetAtomRadiusMethod(), declared);

        RsrOptions rsr_options;
        EXPECT_NO_THROW(rsr_options.SetAtomRadiusMethod(declared));
        EXPECT_EQ(rsr_options.GetAtomRadiusMethod(), declared);
    }
}

TEST(EnumSetterValidationTest, ARejectedAtomRadiusLeavesScoringOnTheDeclaredModel) {
    // The proof that no path can reach an indeterminate radius from the public API: the
    // rejected set leaves the options carrying their previous model, so the scores are
    // bit-identical to a pristine object's and no atom scores NaN. The same sequence
    // returned overall = 1.0 with every by_atom entry NaN before the setter validated.
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    const OESystem::OEScalarGrid obs = MakeGaussianGrid(0.0, 0.0, 0.0, 1.0, 3.0, 0.5);
    const OESystem::OEScalarGrid calc = MakeGaussianGrid(0.0, 0.0, 0.0, 1.0, 3.0, 0.5);

    RsccOptions pristine;
    const DensityScoreResult reference = rscc(mol, obs, 2.0, nullptr, &calc, pristine);

    RsccOptions poked;
    EXPECT_THROW(poked.SetAtomRadiusMethod(static_cast<AtomRadius>(42)), std::invalid_argument);
    const DensityScoreResult after = rscc(mol, obs, 2.0, nullptr, &calc, poked);

    EXPECT_DOUBLE_EQ(after.overall, reference.overall);
    ASSERT_EQ(after.by_atom.size(), reference.by_atom.size());
    for (const auto& [idx, value] : reference.by_atom) {
        ASSERT_EQ(after.by_atom.count(idx), 1u) << "atom " << idx << " went missing";
        EXPECT_FALSE(std::isnan(after.by_atom.at(idx))) << "atom " << idx << " scored NaN";
        EXPECT_DOUBLE_EQ(after.by_atom.at(idx), value);
    }
}

TEST(EnumSetterValidationTest, QScoreRejectsAnUndeclaredRadialSampling) {
    QScoreOptions options;
    const RadialSampling default_value = options.GetRadialSampling();
    for (int undeclared : {42, 999, -1}) {
        EXPECT_THROW(options.SetRadialSampling(static_cast<RadialSampling>(undeclared)),
                     std::invalid_argument)
            << "accepted " << undeclared;
        EXPECT_EQ(options.GetRadialSampling(), default_value)
            << "field moved after rejecting " << undeclared;
    }
}

TEST(EnumSetterValidationTest, StillAcceptsEveryDeclaredRadialSampling) {
    for (RadialSampling declared : {RadialSampling::FIXED, RadialSampling::ADAPTIVE}) {
        QScoreOptions options;
        EXPECT_NO_THROW(options.SetRadialSampling(declared));
        EXPECT_EQ(options.GetRadialSampling(), declared);
    }
}

TEST(EnumSetterValidationTest, ARejectedRadialSamplingLeavesTheFixedSweepGuarded) {
    // An undeclared sampling mode took the `else` arm of both branches in `qscore`: it
    // was validated as adaptive, which consults neither the step nor the maximum radius
    // in the options, and then executed as fixed, which uses both. The fixed sweep's
    // guard was therefore never applied to the parameters the sweep actually ran on.
    // Measured with step 1.0 and max_radius 0.1, that returned overall = nan silently
    // where RadialSampling::FIXED raised GridError.
    QScoreOptions options;
    options.SetRadialStep(1.0);
    options.SetMaxRadius(0.1);
    EXPECT_THROW(options.SetRadialSampling(static_cast<RadialSampling>(42)),
                 std::invalid_argument);

    OESystem::OEScalarGrid grid = MakeGaussianGrid(0.0, 0.0, 0.0, options.GetSigma(), 3.0, 0.5);
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    EXPECT_THROW(qscore(mol, grid, 2.0, nullptr, options), GridError);
}

// ---- CoverageOptions validation ----

TEST(CoverageOptionsValidationTest, RejectsNonFiniteSigma) {
    // A NaN sigma makes `mean + sigma * stddev` NaN, every `rho >= threshold`
    // comparison false, and coverage returns a completely plausible 0.0.
    const double nan_value = std::numeric_limits<double>::quiet_NaN();
    CoverageOptions options;
    const double default_value = options.GetSigma();
    EXPECT_THROW(options.SetSigma(nan_value), std::invalid_argument);
    EXPECT_THROW(options.SetSigma(std::numeric_limits<double>::infinity()),
                 std::invalid_argument);
    EXPECT_THROW(options.SetSigma(-std::numeric_limits<double>::infinity()),
                 std::invalid_argument);
    EXPECT_EQ(options.GetSigma(), default_value);  // Field unchanged after rejection
}

TEST(CoverageOptionsValidationTest, StillAcceptsZeroAndNegativeSigma) {
    // The accepting half of the boundary, and the reason this guard is weaker than
    // QScoreOptions::SetSigma. Coverage's sigma multiplies stddev in
    // `mean + sigma * stddev`: zero thresholds at the mean and a negative value
    // thresholds below it. Both are meaningful requests, so narrowing this guard to
    // match Q-score's would be a regression, not a consistency fix.
    CoverageOptions options;
    EXPECT_NO_THROW(options.SetSigma(0.0));
    EXPECT_EQ(options.GetSigma(), 0.0);
    EXPECT_NO_THROW(options.SetSigma(-1.0));
    EXPECT_EQ(options.GetSigma(), -1.0);
    EXPECT_NO_THROW(options.SetSigma(1.5));
    EXPECT_EQ(options.GetSigma(), 1.5);
}

// ---- Resolution validation at the five public entry points ----

TEST(ResolutionValidationTest, EveryEntryPointRejectsNonFiniteResolution) {
    // NaN <= 0.0 and +inf <= 0.0 are both false, so the old sign-only guard passed
    // them through. fc_density then returned an all-zero grid with no exception.
    const double nan_value = std::numeric_limits<double>::quiet_NaN();
    const double pos_inf = std::numeric_limits<double>::infinity();

    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    const OESystem::OEScalarGrid obs = MakeGaussianGrid(0.0, 0.0, 0.0, 1.0, 3.0, 0.5);
    const OESystem::OEScalarGrid calc = MakeGaussianGrid(0.0, 0.0, 0.0, 1.0, 3.0, 0.5);

    for (double bad : {nan_value, pos_inf, -pos_inf, 0.0, -1.0}) {
        EXPECT_THROW(rscc(mol, obs, bad, nullptr, &calc), GridError) << "resolution " << bad;
        EXPECT_THROW(rsr(mol, obs, bad, nullptr, &calc), GridError) << "resolution " << bad;
        EXPECT_THROW(qscore(mol, obs, bad), GridError) << "resolution " << bad;
        EXPECT_THROW(ediam(mol, obs, bad), GridError) << "resolution " << bad;

        DensityCalculator density(UnitCell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0),
                                  SymOp::ParseAll("x,y,z"));
        EXPECT_THROW(density.Calculate(mol, obs, bad), GridError) << "resolution " << bad;
    }
}

TEST(ResolutionValidationTest, StillAcceptsAFinitePositiveResolution) {
    // The accepting half: one shared guard must not narrow what the five entry
    // points take.
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    const OESystem::OEScalarGrid obs = MakeGaussianGrid(0.0, 0.0, 0.0, 1.0, 3.0, 0.5);
    const OESystem::OEScalarGrid calc = MakeGaussianGrid(0.0, 0.0, 0.0, 1.0, 3.0, 0.5);

    EXPECT_NO_THROW(rscc(mol, obs, 2.0, nullptr, &calc));
    EXPECT_NO_THROW(rsr(mol, obs, 2.0, nullptr, &calc));
    EXPECT_NO_THROW(qscore(mol, obs, 2.0));
    EXPECT_NO_THROW(ediam(mol, obs, 2.0));

    DensityCalculator density(UnitCell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0),
                              SymOp::ParseAll("x,y,z"));
    EXPECT_NO_THROW(delete density.Calculate(mol, obs, 2.0));
}

// ---- DensityCalculator numeric argument bounds ----

TEST(DensityCalculatorValidationTest, RejectsAnUnboundedScaleShellCount) {
    // n_scale_shells is unsigned: at UINT_MAX, `n_scale_shells + 1` wraps to 0, the
    // shell-edge vector is empty, and `i <= UINT_MAX` never ends -- a non-terminating
    // loop writing past the end of an empty vector on every iteration.
    OEChem::OEGraphMol mol = MakeAtomMol(6, 5.0, 5.0, 5.0);
    const OESystem::OEScalarGrid obs = MakeGaussianGrid(5.0, 5.0, 5.0, 1.0, 3.0, 1.0);
    DensityCalculator density(UnitCell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0),
                              SymOp::ParseAll("x,y,z"));

    EXPECT_THROW(density.Calculate(mol, obs, 2.0, nullptr, 0.35, 46.0, false, 0), GridError);
    EXPECT_THROW(density.Calculate(mol, obs, 2.0, nullptr, 0.35, 46.0, false,
                                   MAX_SCALE_SHELLS + 1),
                 GridError);
    EXPECT_THROW(density.Calculate(mol, obs, 2.0, nullptr, 0.35, 46.0, false,
                                   std::numeric_limits<unsigned int>::max()),
                 GridError);
}

TEST(DensityCalculatorValidationTest, StillAcceptsScaleShellCountsUpToTheLimit) {
    // The accepting half of the bound. MAX_SCALE_SHELLS itself must work, or the
    // guard has narrowed the usable range rather than closing the wrap.
    OEChem::OEGraphMol mol = MakeAtomMol(6, 5.0, 5.0, 5.0);
    const OESystem::OEScalarGrid obs = MakeGaussianGrid(5.0, 5.0, 5.0, 1.0, 3.0, 1.0);
    DensityCalculator density(UnitCell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0),
                              SymOp::ParseAll("x,y,z"));

    EXPECT_NO_THROW(delete density.Calculate(mol, obs, 2.0, nullptr, 0.35, 46.0, false, 1));
    EXPECT_NO_THROW(delete density.Calculate(mol, obs, 2.0, nullptr, 0.35, 46.0, false, 4));
    EXPECT_NO_THROW(
        delete density.Calculate(mol, obs, 2.0, nullptr, 0.35, 46.0, false, MAX_SCALE_SHELLS));
}

TEST(DensityCalculatorValidationTest, RejectsASpacingCoarserThanTheCell) {
    // round(2.0 / 5.0) is 0, so the FFT grid would be 0 x 0 x 0: `% nx` is a division
    // by zero (SIGFPE on x86-64) and the allocation is zero-sized. Every other guard
    // on this path passes.
    //
    // Assert on the message, not just the type. Without this guard the zero-sized
    // geometry still reaches FFTW and dies there with "FFTW planning failed", so a
    // bare EXPECT_THROW(GridError) passes on arm64 whether the guard exists or not.
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    const OESystem::OEScalarGrid obs = MakeGaussianGrid(0.0, 0.0, 0.0, 1.0, 5.0, 5.0);
    DensityCalculator density(UnitCell(2.0, 2.0, 2.0, 90.0, 90.0, 90.0),
                              SymOp::ParseAll("x,y,z"));

    try {
        delete density.Calculate(mol, obs, 2.0);
        FAIL() << "expected GridError for a spacing coarser than the cell";
    } catch (const GridError& error) {
        EXPECT_NE(std::string(error.what()).find("too coarse"), std::string::npos)
            << "rejected by the wrong branch: " << error.what();
    }
}

TEST(DensityCalculatorValidationTest, StillAcceptsASpacingThatRoundsToOnePoint) {
    // One point per axis is the smallest usable grid; the guard rejects below it, not
    // at it. cell 2 A at 1.5 A spacing rounds to 1.
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    const OESystem::OEScalarGrid obs = MakeGaussianGrid(0.0, 0.0, 0.0, 1.0, 4.5, 1.5);
    DensityCalculator density(UnitCell(2.0, 2.0, 2.0, 90.0, 90.0, 90.0),
                              SymOp::ParseAll("x,y,z"));

    EXPECT_NO_THROW(delete density.Calculate(mol, obs, 2.0));
}

TEST(DensityCalculatorValidationTest, RejectsAResolutionThatExplodesTheMillerBox) {
    // A tiny but finite resolution passes require_usable_resolution and then sizes a
    // triple loop as (2*ceil(a/resolution)+1)^3. Every value here hung the process
    // before the bound existed, and 1e-9 A is an ordinary double rather than a
    // subnormal -- the loop diverges long before 1/resolution overflows, so this is
    // not a subnormal or an overflow special case.
    //
    // Assert on the message, not just the type: several other GridError sites sit on
    // this path, and a bare EXPECT_THROW could not tell the box bound from any of them.
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    const OESystem::OEScalarGrid obs = MakeGaussianGrid(0.0, 0.0, 0.0, 1.0, 3.0, 0.5);
    DensityCalculator density(UnitCell(2.0, 2.0, 2.0, 90.0, 90.0, 90.0),
                              SymOp::ParseAll("x,y,z"));

    for (double tiny : {1e-9, 1e-300, std::numeric_limits<double>::denorm_min()}) {
        try {
            delete density.Calculate(mol, obs, tiny);
            FAIL() << "expected GridError for resolution " << tiny;
        } catch (const GridError& error) {
            EXPECT_NE(std::string(error.what()).find("Miller-index box"), std::string::npos)
                << "rejected by the wrong branch at resolution " << tiny << ": " << error.what();
        }
    }
}

TEST(DensityCalculatorValidationTest, StillAcceptsRealCrystallographicMillerBoxes) {
    // The accepting half of the bound, in two parts.
    //
    // The two reference geometries are asserted as arithmetic on the constant rather
    // than run end to end: a 200 A cell at 1.0 A is a legal request that allocates
    // gigabytes, so calling Calculate on it would make this a memory test rather than a
    // bound test. Narrowing MAX_MILLER_BOX_POINTS to anything under 6.5e7 -- the kind of
    // "tighten it a bit" edit this guard invites -- fails here. The expressions mirror
    // GenerateMillerIndices' own `a * (1.0 / resolution)` so a rounding difference in
    // ceil cannot make the test disagree with the code it pins.
    const auto box_points = [](double edge, double resolution) {
        const double extent = std::ceil(edge * (1.0 / resolution));
        return (2.0 * extent + 1.0) * (2.0 * extent + 1.0) * (2.0 * extent + 1.0);
    };
    EXPECT_LE(box_points(200.0, 1.0), MAX_MILLER_BOX_POINTS);
    EXPECT_LE(box_points(100.0, 0.8), MAX_MILLER_BOX_POINTS);

    // And one cell small enough to run all the way through, so the guard is shown to
    // pass real work rather than merely to hold a large number.
    OEChem::OEGraphMol mol = MakeAtomMol(6, 5.0, 5.0, 5.0);
    const OESystem::OEScalarGrid obs = MakeGaussianGrid(5.0, 5.0, 5.0, 1.0, 3.0, 1.0);
    DensityCalculator density(UnitCell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0),
                              SymOp::ParseAll("x,y,z"));
    EXPECT_NO_THROW(delete density.Calculate(mol, obs, 0.8));
}

TEST(DensityCalculatorValidationTest, KeepsEveryMillerIndexInsideTheResolutionShell) {
    // MAX_MILLER_BOX_POINTS bounds the product of the three extents, not any one of
    // them, so a cell with one long edge and two short ones drives h far past the point
    // where `h * h` fits in an int while the box stays orders of magnitude under the
    // limit. s2 was accumulated from int products; the overflow wrapped it negative,
    // which passes `s2 <= s_max2` and admitted reflections from outside the shell. The
    // grid that came back looked plausible and was wrong, so the failure this pins is
    // silent -- nothing throws, and only the returned indices show it.
    //
    // a = 50000 A at 1.0 A gives extents 50000, 1, 1 and a box of 100001 * 3 * 3 =
    // 900009 points, about 222x under MAX_MILLER_BOX_POINTS. No tightening of that
    // constant reaches this case, which is why the fix is in the arithmetic instead.
    const double a = 50000.0, b = 1.0, c = 1.0, resolution = 1.0;
    const double s_max2 = (1.0 / resolution) * (1.0 / resolution);

    const auto indices = detail::GenerateMillerIndices(a, b, c, resolution);
    ASSERT_FALSE(indices.empty());

    // Largest |h| whose square still fits in an int, floor(sqrt(INT_MAX)) = 46340.
    const int MAX_NON_OVERFLOWING_INDEX =
        static_cast<int>(std::sqrt(static_cast<double>(std::numeric_limits<int>::max())));

    // Recomputed in double from h, k, l rather than read back from `stol2`, which is
    // derived from the very expression under test and would let a wrapped value certify
    // itself.
    int out_of_shell = 0;
    int max_abs_h = 0;
    detail::MillerIndex worst{0, 0, 0, 0.0};
    double worst_s2 = 0.0;
    for (const auto& index : indices) {
        const int abs_h = index.h < 0 ? -index.h : index.h;
        if (abs_h > max_abs_h) max_abs_h = abs_h;

        const double hd = index.h, kd = index.k, ld = index.l;
        const double s2 = (hd * hd) / (a * a) + (kd * kd) / (b * b) + (ld * ld) / (c * c);
        if (s2 > s_max2) {
            if (out_of_shell == 0) {
                worst = index;
                worst_s2 = s2;
            }
            ++out_of_shell;
        }
    }

    // Without an index past the int-square limit the overflow cannot occur and the
    // assertion below would hold for a reason unrelated to the defect. Asserted rather
    // than expected so a cell edited to a safe size fails here instead of passing
    // vacuously downstream.
    ASSERT_GT(max_abs_h, MAX_NON_OVERFLOWING_INDEX)
        << "cell no longer drives an extent past " << MAX_NON_OVERFLOWING_INDEX
        << ", so this test cannot observe the overflow it exists to pin";

    EXPECT_EQ(out_of_shell, 0)
        << out_of_shell << " of " << indices.size() << " returned indices lie outside the "
        << resolution << " A shell; first is (" << worst.h << ", " << worst.k << ", " << worst.l
        << ") at s2 = " << worst_s2 << ", over the limit of " << s_max2;
}

// ---- wrap_and_pad_grid cell-edge validation ----

TEST(WrapAndPadValidationTest, RejectsZeroOrNonFiniteCellEdges) {
    // fmod(x, 0.0) is NaN, so a zero edge filled the padded grid with NaN voxels and
    // returned it. The old ternary guarded the centroid shift only.
    const double nan_value = std::numeric_limits<double>::quiet_NaN();
    const double pos_inf = std::numeric_limits<double>::infinity();

    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    const OESystem::OESkewGrid grid(MakeEmptyGrid(5.0, 1.0));

    EXPECT_THROW(wrap_and_pad_grid(grid, mol, 0.0, 5.0, 5.0), CellError);
    EXPECT_THROW(wrap_and_pad_grid(grid, mol, 5.0, 0.0, 5.0), CellError);
    EXPECT_THROW(wrap_and_pad_grid(grid, mol, 5.0, 5.0, 0.0), CellError);
    EXPECT_THROW(wrap_and_pad_grid(grid, mol, nan_value, 5.0, 5.0), CellError);
    EXPECT_THROW(wrap_and_pad_grid(grid, mol, pos_inf, 5.0, 5.0), CellError);
    EXPECT_THROW(wrap_and_pad_grid(grid, mol, -5.0, 5.0, 5.0), CellError);
}

TEST(WrapAndPadValidationTest, StillAcceptsPositiveCellEdges) {
    // The atom sits within `padding` of the grid edge, so this actually reaches the
    // periodic sampling rather than returning nullptr for "no padding needed" -- the
    // rejection test above throws before that point and cannot cover it. The grid is
    // 11 nodes at a 1.0 A interval, so 11.0 is the cell edge it samples; the centroid
    // sits 4.0 A from the grid centre, inside half a cell, so no shift runs and the
    // atom stays near the edge where padding is genuinely required.
    OEChem::OEGraphMol mol = MakeAtomMol(6, 4.0, 4.0, 4.0);
    OESystem::OESkewGrid grid(MakeEmptyGrid(5.0, 1.0));
    float* values = grid.GetValues();
    for (unsigned int i = 0; i < grid.GetSize(); ++i) values[i] = 1.0f;

    std::unique_ptr<OESystem::OESkewGrid> padded;
    ASSERT_NO_THROW(padded.reset(wrap_and_pad_grid(grid, mol, 11.0, 11.0, 11.0)));
    ASSERT_NE(padded, nullptr) << "expected the padding path, not the nullptr shortcut";
    const float* padded_values = padded->GetValues();
    for (unsigned int i = 0; i < padded->GetSize(); ++i) {
        ASSERT_FALSE(std::isnan(padded_values[i])) << "NaN voxel at " << i;
    }
}

// ---- wrap_and_pad_grid padding and interval-count validation ----

TEST(WrapAndPadValidationTest, RejectsNonFiniteOrNegativePadding) {
    // padding grows the extent the padded grid is sized from, and that extent divides
    // by the node interval into a float-to-unsigned conversion, which is undefined
    // outside the range the type represents. Of the values below only +inf reached
    // that conversion before this guard existed: for this molecule the other four
    // make every needs_pad comparison false, so the function returned nullptr and
    // reported that the atoms already fit. Turning that silent answer into an error
    // is the behaviour this changes.
    //
    // Assert on the message. +inf does reach the sizing loop and is already rejected
    // there on this arm64 host -- see the interval-count case below -- so the
    // exception type alone cannot say which guard fired.
    const double nan_value = std::numeric_limits<double>::quiet_NaN();
    const double pos_inf = std::numeric_limits<double>::infinity();

    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    const OESystem::OESkewGrid grid(MakeEmptyGrid(5.0, 1.0));  // 11 nodes at 1.0 A, cell 11

    // -1e-9 is here because the contract is non-negative rather than "not badly
    // negative": a negative padding shrinks the box the atoms have to fit inside,
    // whatever its magnitude. Zero stays admissible and is covered by
    // GridOpsTest.WrapAndPadGridRejectsAPaddingTooSmallForTheNodeInterval.
    for (const double bad : {nan_value, pos_inf, -pos_inf, -5.0, -1e-9}) {
        try {
            std::unique_ptr<OESystem::OESkewGrid> padded(
                wrap_and_pad_grid(grid, mol, 11.0, 11.0, 11.0, bad));
            FAIL() << "expected GridError for a padding of " << bad;
        } catch (const GridError& error) {
            EXPECT_NE(std::string(error.what()).find("finite non-negative padding"),
                      std::string::npos)
                << "rejected by the wrong branch at padding " << bad << ": " << error.what();
        }
    }
}

TEST(WrapAndPadValidationTest, RejectsAnIntervalCountNoGridDimensionCanHold) {
    // A valid padding is not sufficient: the atom extent is caller-controlled too, and
    // it reaches the same conversion. The two atoms below are 5e9 A apart, a separation
    // the two targets convert differently. Compiling the conversion for both: on this
    // arm64 host 5e9 saturates to UINT_MAX, whose `+ 1u` wraps to 0 and is caught by
    // the two-node guard; built for x86-64 it wraps modularly to 705032704, and the
    // dimension of 705032705 that follows passes that guard and reaches SetDim.
    //
    // So assert on the message: on this host the type is GridError with the bound
    // removed as well as with it, and only the message says which guard fired.
    //
    // The atoms straddle the grid centre, so no coordinate shift runs and the extent
    // under test is the pair's own separation.
    OEChem::OEGraphMol mol = MakeAtomMol(6, -2.5e9, 0.0, 0.0);
    OEChem::OEAtomBase* far_atom = mol.NewAtom(6);
    const double far_coords[3] = {2.5e9, 0.0, 0.0};
    mol.SetCoords(far_atom, far_coords);

    const OESystem::OESkewGrid grid(MakeEmptyGrid(5.0, 1.0));  // 11 nodes at 1.0 A, cell 11

    try {
        std::unique_ptr<OESystem::OESkewGrid> padded(
            wrap_and_pad_grid(grid, mol, 11.0, 11.0, 11.0, 0.0));
        FAIL() << "expected GridError for a 5e9 A atom extent at a 1 A node interval";
    } catch (const GridError& error) {
        EXPECT_NE(std::string(error.what()).find("no grid dimension can hold"),
                  std::string::npos)
            << "rejected by the wrong branch: " << error.what();
    }
}

TEST(WrapAndPadValidationTest, RejectsAPaddedCellEdgeNoFloatCanHold) {
    // The interval count is not the only quantity the sizing narrows. The padded
    // grid's cell edge is pad_dim times the source node interval, and OESkewGrid
    // takes it as a float. A count well inside MAX_PAD_INTERVALS still drives that
    // conversion out of range when the source grid is sampled coarsely enough: the
    // count divides the extent by the interval and the edge multiplies it back, so
    // the count bound does not constrain the edge.
    //
    // The source grid is two nodes at a 2^126 A interval, which is a cell edge of
    // 2^127 A that a float holds exactly. A padding of 2^127 A gives an extent of
    // 2^128 A, four intervals, five nodes, and a padded edge of 5 * 2^126 =
    // 4.25e38 A, past the 3.40e38 a float holds -- while the count of 4 is nine
    // orders inside the interval bound. The cell edge, the node interval and the
    // padding are all powers of two, so the interval the node walk derives is
    // exact and the count is exactly 4.
    //
    // Assert on the message. The two-node guard and the interval-count guard raise
    // GridError from this same loop, so the type alone cannot say which fired.
    const float cell_edge = 0x1p127f;  // 1.70141e38 A over two nodes
    OESystem::OESkewGrid grid;
    ASSERT_TRUE(grid.SetDim(2u, 2u, 2u));
    ASSERT_TRUE(grid.SetUnitCell(cell_edge, cell_edge, cell_edge,
                                 90.0f, 90.0f, 90.0f, 2u, 2u, 2u));
    ASSERT_TRUE(grid.SetMid(0.0f, 0.0f, 0.0f));

    // The atom sits at the grid centre, so no coordinate shift runs and the extent
    // under test is the padding's alone.
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);

    try {
        std::unique_ptr<OESystem::OESkewGrid> padded(
            wrap_and_pad_grid(grid, mol, cell_edge, cell_edge, cell_edge, 0x1p127));
        FAIL() << "expected GridError for a padded cell edge past the float maximum";
    } catch (const GridError& error) {
        EXPECT_NE(std::string(error.what()).find("which no float can hold"),
                  std::string::npos)
            << "rejected by the wrong branch: " << error.what();
    }
}

// ---- wrap_and_pad_grid centroid-shift validation ----

namespace {

/// Every atom's coordinates in iteration order, for asserting that a rejected
/// call left the caller's molecule alone.
std::vector<float> AtomCoords(OEChem::OEMolBase& mol) {
    std::vector<float> out;
    float coords[3];
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        mol.GetCoords(&(*atom), coords);
        out.push_back(coords[0]);
        out.push_back(coords[1]);
        out.push_back(coords[2]);
    }
    return out;
}

/// A two-node cubic grid with the given cell edge, centred at `xmid` on x and
/// at the origin on the other two axes.
OESystem::OESkewGrid MakeHugeCellGrid(const float cell_edge, const float xmid) {
    OESystem::OESkewGrid grid;
    EXPECT_TRUE(grid.SetDim(2u, 2u, 2u));
    EXPECT_TRUE(grid.SetUnitCell(cell_edge, cell_edge, cell_edge,
                                 90.0f, 90.0f, 90.0f, 2u, 2u, 2u));
    EXPECT_TRUE(grid.SetMid(xmid, 0.0f, 0.0f));
    return grid;
}

}  // namespace

TEST(WrapAndPadValidationTest, RejectsACentroidShiftNoFloatCoordinateCanHold) {
    // The shift is computed in double and stored back as a float, and nothing
    // bounded that narrowing. The atom below sits half a cell short of the grid
    // centre, so std::round takes the shift to one full cell edge and the sum
    // lands past the float maximum. Without this guard the conversion saturated
    // on this host, the infinity was written into the caller's molecule, and what
    // the caller saw was the sizing guard rejecting the infinite extent that
    // produced -- an error raised over coordinates the call had already replaced.
    //
    // The grid centre is placed below the float maximum by half the cell edge's
    // node interval, so the grid's own node coordinates stay finite and the shift
    // is the only thing that overflows.
    const float cell_edge = 0x1p127f;  // 1.70141e38 A over two nodes
    const float grid_mid = std::numeric_limits<float>::max() - 0x1p125f;
    const OESystem::OESkewGrid grid = MakeHugeCellGrid(cell_edge, grid_mid);

    const double atom_x =
        static_cast<double>(grid_mid) - 0.5 * static_cast<double>(cell_edge);
    OEChem::OEGraphMol mol = MakeAtomMol(6, atom_x, 0.0, 0.0);
    const std::vector<float> before = AtomCoords(mol);

    // Assert on the message: the sizing guards below raise GridError from this same
    // function, so the type alone cannot say which one fired.
    try {
        std::unique_ptr<OESystem::OESkewGrid> padded(
            wrap_and_pad_grid(grid, mol, cell_edge, cell_edge, cell_edge, 0.0));
        FAIL() << "expected GridError for a centroid shift past the float maximum";
    } catch (const GridError& error) {
        EXPECT_NE(std::string(error.what()).find("would shift atom"), std::string::npos)
            << "rejected by the wrong branch: " << error.what();
    }

    // The reason the guard exists: a rejected call must leave the caller's
    // coordinates as they were rather than saturated.
    EXPECT_EQ(AtomCoords(mol), before) << "the rejected call moved the molecule";
}

TEST(WrapAndPadValidationTest, RejectsTheShiftBeforeMovingAnEarlierAtomThatFits) {
    // A guard that tested each coordinate as it wrote it would pass the single-atom
    // case above and still leave a multi-atom molecule half-shifted. Atom 0 sits at
    // the Cartesian origin and atom 1 at the grid centre, so their centroid is
    // between half a cell and one and a half cells short of that centre and the
    // shift is again one full cell edge. Atom 0 lands exactly on that edge, a value
    // a float holds; atom 1 lands past the float maximum. Atom 0 is written first,
    // so it is what distinguishes rejecting before the first write from rejecting
    // during the walk.
    const float cell_edge = 0x1p127f;
    const float grid_mid = std::numeric_limits<float>::max() - 0x1p125f;
    const OESystem::OESkewGrid grid = MakeHugeCellGrid(cell_edge, grid_mid);

    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    OEChem::OEAtomBase* far_atom = mol.NewAtom(6);
    const double far_coords[3] = {static_cast<double>(grid_mid), 0.0, 0.0};
    mol.SetCoords(far_atom, far_coords);
    const std::vector<float> before = AtomCoords(mol);

    try {
        std::unique_ptr<OESystem::OESkewGrid> padded(
            wrap_and_pad_grid(grid, mol, cell_edge, cell_edge, cell_edge, 0.0));
        FAIL() << "expected GridError for a centroid shift past the float maximum";
    } catch (const GridError& error) {
        EXPECT_NE(std::string(error.what()).find("would shift atom"), std::string::npos)
            << "rejected by the wrong branch: " << error.what();
    }

    EXPECT_EQ(AtomCoords(mol), before)
        << "the rejected call moved at least one atom; the shift is not all-or-nothing";
}

// ---- Per-axis grid geometry derivation (spec §2.2, §4.4) ----

namespace {

/// A skew grid with the given dimensions and an orthorhombic cell whose edges
/// give a 1.0 A node interval on every axis.
OESystem::OESkewGrid MakeSkewGrid(const unsigned int nx,
                                  const unsigned int ny,
                                  const unsigned int nz) {
    OESystem::OESkewGrid grid;
    EXPECT_TRUE(grid.SetDim(nx, ny, nz));
    EXPECT_TRUE(grid.SetUnitCell(static_cast<float>(nx), static_cast<float>(ny),
                                 static_cast<float>(nz), 90.0f, 90.0f, 90.0f,
                                 nx, ny, nz));
    EXPECT_TRUE(grid.SetMid(0.0f, 0.0f, 0.0f));
    return grid;
}

}  // namespace

TEST(GridParamsDerivation, RejectsAnAxisWithASingleNode) {
    // Check 1: one node on an axis leaves no interval to measure.
    const OESystem::OESkewGrid grid = MakeSkewGrid(4u, 4u, 1u);
    EXPECT_THROW(get_grid_params(grid), GridError);
    try {
        get_grid_params(grid);
        FAIL() << "expected GridError";
    } catch (const GridError& err) {
        const std::string what = err.what();
        EXPECT_NE(what.find("axis z"), std::string::npos) << what;
        EXPECT_NE(what.find("1"), std::string::npos) << what;
    }
}

TEST(GridParamsDerivation, RejectsANonFiniteNodeCoordinate) {
    // Check 3: a zero a-edge makes the cell matrix singular, and every node
    // coordinate comes back NaN. ElementToSpatialCoord still reports success,
    // so the finiteness test is what catches this, not the return value.
    OESystem::OESkewGrid grid;
    ASSERT_TRUE(grid.SetDim(4u, 4u, 4u));
    grid.SetUnitCell(0.0f, 4.0f, 4.0f, 90.0f, 90.0f, 90.0f, 4u, 4u, 4u);
    grid.SetMid(0.0f, 0.0f, 0.0f);

    // Assert the state the check needs, so a future toolkit that rejects the
    // degenerate cell fails here with a readable message rather than below.
    float x = 0.0f, y = 0.0f, z = 0.0f;
    ASSERT_TRUE(grid.ElementToSpatialCoord(0u, x, y, z));
    ASSERT_FALSE(std::isfinite(x));

    EXPECT_THROW(get_grid_params(grid), GridError);
    try {
        get_grid_params(grid);
        FAIL() << "expected GridError";
    } catch (const GridError& err) {
        const std::string what = err.what();
        EXPECT_NE(what.find("non-finite"), std::string::npos) << what;
        EXPECT_NE(what.find("element 0"), std::string::npos) << what;
    }
}

TEST(GridParamsDerivation, RejectsAZeroDerivedNodeInterval) {
    // Check 5: a 1e-30 A a-edge is small enough that all four x nodes round to
    // the same float, so the x span is exactly zero -- but the coordinates stay
    // finite, so check 3 does not fire. The off-axis leak is ~1e-15, four orders
    // under the 1e-4 limit, so check 4 does not fire either.
    OESystem::OESkewGrid grid;
    ASSERT_TRUE(grid.SetDim(4u, 4u, 4u));
    grid.SetUnitCell(1e-30f, 4.0f, 4.0f, 90.0f, 90.0f, 90.0f, 4u, 4u, 4u);
    grid.SetMid(0.0f, 0.0f, 0.0f);

    float x0 = 0.0f, y0 = 0.0f, z0 = 0.0f;
    float xf = 0.0f, yf = 0.0f, zf = 0.0f;
    ASSERT_TRUE(grid.ElementToSpatialCoord(0u, x0, y0, z0));
    ASSERT_TRUE(grid.ElementToSpatialCoord(3u, xf, yf, zf));
    ASSERT_TRUE(std::isfinite(x0));
    ASSERT_EQ(x0, xf) << "the x span must be exactly zero for this to reach check 5";

    EXPECT_THROW(get_grid_params(grid), GridError);
    try {
        get_grid_params(grid);
        FAIL() << "expected GridError";
    } catch (const GridError& err) {
        const std::string what = err.what();
        EXPECT_NE(what.find("axis x"), std::string::npos) << what;
        EXPECT_NE(what.find("finite and positive"), std::string::npos) << what;
    }
}

// Derivation check 2 (ElementToSpatialCoord returning false) has no test.
//   - It is the only one of the five with no reachable public construction.
//   - Probed against OpenEye 2026.1.0: the call returned true on a default
//     1x1x1 grid, on a well-formed 4x4x4, and on both degenerate grids above --
//     including the one whose coordinates are all NaN. (Task 1 report.)
//   - The check stays in get_grid_params because a grid arriving from a reader
//     rather than from the setters, or a future toolkit release, could reach it.

TEST(GridParamsDerivation, RejectsANonAxisAlignedCell) {
    // Check 4: a 60-degree gamma tilts the b axis into x, so a per-axis
    // spacing cannot describe the sampling.
    OESystem::OESkewGrid grid;
    ASSERT_TRUE(grid.SetDim(10u, 10u, 10u));
    ASSERT_TRUE(grid.SetUnitCell(10.0f, 10.0f, 10.0f, 90.0f, 90.0f, 60.0f,
                                 10u, 10u, 10u));
    ASSERT_TRUE(grid.SetMid(0.0f, 0.0f, 0.0f));
    EXPECT_THROW(get_grid_params(grid), CellError);
    try {
        get_grid_params(grid);
        FAIL() << "expected CellError";
    } catch (const CellError& err) {
        const std::string what = err.what();
        EXPECT_NE(what.find("axis y"), std::string::npos) << what;
        EXPECT_NE(what.find("into x"), std::string::npos) << what;
    }
}

TEST(GridParamsDerivation, ReportsTheFirstFailingCheckNotTheLast) {
    // A default-constructed grid is 1x1x1, so it fails check 1 on every axis;
    // it would also fail check 5, because a zero span over zero intervals
    // derives a NaN spacing. The check-1 message is the one that must surface.
    const OESystem::OESkewGrid grid;
    try {
        get_grid_params(grid);
        FAIL() << "expected GridError";
    } catch (const GridError& err) {
        const std::string what = err.what();
        EXPECT_NE(what.find("axis x"), std::string::npos) << what;
        EXPECT_NE(what.find("at least 2 nodes"), std::string::npos) << what;
    }
}

TEST(SameGridGeometry, TrueForTwoGridsWithTheSameSampling) {
    const OESystem::OESkewGrid lhs = MakeSkewGrid(4u, 5u, 6u);
    const OESystem::OESkewGrid rhs = MakeSkewGrid(4u, 5u, 6u);
    EXPECT_TRUE(same_grid_geometry(lhs, rhs));
}

TEST(SameGridGeometry, FalseWhenADimensionDiffers) {
    const OESystem::OESkewGrid lhs = MakeSkewGrid(4u, 5u, 6u);
    const OESystem::OESkewGrid rhs = MakeSkewGrid(4u, 5u, 7u);
    EXPECT_FALSE(same_grid_geometry(lhs, rhs));
}

TEST(SameGridGeometry, ThrowsRatherThanReportingDifferentForAnUnderivableGrid) {
    // Both operands go through get_grid_params, so a grid whose geometry cannot
    // be derived is an error rather than a "these differ" answer. This is the
    // documented departure from OEGridSameGeometry, which returned false.
    const OESystem::OESkewGrid lhs = MakeSkewGrid(4u, 5u, 6u);
    const OESystem::OESkewGrid rhs;  // default-constructed: 1x1x1
    EXPECT_THROW(same_grid_geometry(lhs, rhs), GridError);
}

TEST(SameGridGeometry, TrueForOriginDifferenceBelowRelativeTolerance) {
    // The two mids differ by one representable step near 50: 50.00002f is the
    // float 50.0000190734863, and the derivation carries that rigidly to the
    // node origins at 45.5 and 45.5000190734863. The resulting 1.907e-5 A gap
    // sits strictly between the two tolerance regimes -- 19.1x over an absolute
    // 1e-6, and 2.4x inside the relative 1e-6 * 45.5 = 4.55e-5. Measured, not
    // derived: an absolute comparison reports these two grids as different.
    OESystem::OESkewGrid lhs;
    ASSERT_TRUE(lhs.SetDim(10u, 10u, 10u));
    ASSERT_TRUE(lhs.SetUnitCell(10.0f, 10.0f, 10.0f, 90.0f, 90.0f, 90.0f,
                                10u, 10u, 10u));
    ASSERT_TRUE(lhs.SetMid(50.0f, 0.0f, 0.0f));

    OESystem::OESkewGrid rhs;
    ASSERT_TRUE(rhs.SetDim(10u, 10u, 10u));
    ASSERT_TRUE(rhs.SetUnitCell(10.0f, 10.0f, 10.0f, 90.0f, 90.0f, 90.0f,
                                10u, 10u, 10u));
    ASSERT_TRUE(rhs.SetMid(50.00002f, 0.0f, 0.0f));

    EXPECT_TRUE(same_grid_geometry(lhs, rhs));
}

TEST(SameGridGeometry, FalseWhenOnlyNodeOriginDiffers) {
    // Two grids with identical dims and cells, mids at (0,0,0) and (5,0,0).
    // This proves the origin comparison at src/Grid.cpp:189 is reachable: dims,
    // spacings match, but origins differ by 5 A (far outside any tolerance).
    OESystem::OESkewGrid lhs;
    ASSERT_TRUE(lhs.SetDim(10u, 10u, 10u));
    ASSERT_TRUE(lhs.SetUnitCell(10.0f, 10.0f, 10.0f, 90.0f, 90.0f, 90.0f,
                                10u, 10u, 10u));
    ASSERT_TRUE(lhs.SetMid(0.0f, 0.0f, 0.0f));

    OESystem::OESkewGrid rhs;
    ASSERT_TRUE(rhs.SetDim(10u, 10u, 10u));
    ASSERT_TRUE(rhs.SetUnitCell(10.0f, 10.0f, 10.0f, 90.0f, 90.0f, 90.0f,
                                10u, 10u, 10u));
    ASSERT_TRUE(rhs.SetMid(5.0f, 0.0f, 0.0f));

    const GridParams a = get_grid_params(lhs);
    const GridParams b = get_grid_params(rhs);

    // Assert preconditions: dims and spacings match, origins differ substantially.
    ASSERT_EQ(a.x_dim, b.x_dim);
    ASSERT_EQ(a.y_dim, b.y_dim);
    ASSERT_EQ(a.z_dim, b.z_dim);
    ASSERT_DOUBLE_EQ(a.x_spacing, b.x_spacing);
    ASSERT_DOUBLE_EQ(a.y_spacing, b.y_spacing);
    ASSERT_DOUBLE_EQ(a.z_spacing, b.z_spacing);
    ASSERT_GT(std::abs(a.x_origin - b.x_origin), 1.0);

    EXPECT_FALSE(same_grid_geometry(lhs, rhs));
}

TEST(SameGridGeometry, FalseWhenOnlyCellParametersDiffer) {
    // Two grids with identical dims and mid, cells with a-edge 10 vs 20 and
    // matching divisions (10 vs 20). This keeps the x node interval at 1.0 A
    // for both grids, and with the same mid and dims the node origins coincide
    // too. Dims, spacings, origins all match; only cell a differs (10 vs 20).
    // This proves the cell-parameter comparison at src/Grid.cpp:197 is reachable.
    OESystem::OESkewGrid lhs;
    ASSERT_TRUE(lhs.SetDim(10u, 10u, 10u));
    ASSERT_TRUE(lhs.SetUnitCell(10.0f, 10.0f, 10.0f, 90.0f, 90.0f, 90.0f,
                                10u, 10u, 10u));
    ASSERT_TRUE(lhs.SetMid(0.0f, 0.0f, 0.0f));

    OESystem::OESkewGrid rhs;
    ASSERT_TRUE(rhs.SetDim(10u, 10u, 10u));
    ASSERT_TRUE(rhs.SetUnitCell(20.0f, 10.0f, 10.0f, 90.0f, 90.0f, 90.0f,
                                20u, 10u, 10u));
    ASSERT_TRUE(rhs.SetMid(0.0f, 0.0f, 0.0f));

    const GridParams a = get_grid_params(lhs);
    const GridParams b = get_grid_params(rhs);

    // Assert preconditions: dims match, spacings match (both 1.0 A on all axes),
    // origins match, both have unit cells, and cell a differs by ~10.0.
    ASSERT_EQ(a.x_dim, b.x_dim);
    ASSERT_EQ(a.y_dim, b.y_dim);
    ASSERT_EQ(a.z_dim, b.z_dim);
    ASSERT_DOUBLE_EQ(a.x_spacing, b.x_spacing);
    ASSERT_DOUBLE_EQ(a.y_spacing, b.y_spacing);
    ASSERT_DOUBLE_EQ(a.z_spacing, b.z_spacing);
    ASSERT_DOUBLE_EQ(a.x_origin, b.x_origin);
    ASSERT_DOUBLE_EQ(a.y_origin, b.y_origin);
    ASSERT_DOUBLE_EQ(a.z_origin, b.z_origin);
    ASSERT_TRUE(lhs.HasUnitCell());
    ASSERT_TRUE(rhs.HasUnitCell());
    ASSERT_NEAR(std::abs(get_unit_cell(lhs).a - get_unit_cell(rhs).a), 10.0, 0.1);

    EXPECT_FALSE(same_grid_geometry(lhs, rhs));
}

TEST(SameGridGeometry, FalseWhenOnlyNodeSpacingDiffers) {
    // Doubling only the x divisions halves the x node interval while leaving the
    // six cell parameters byte-identical. The compensating mid shift restores the
    // x origin. This proves the spacing comparison at src/Grid.cpp:186 is reachable.
    OESystem::OESkewGrid lhs;
    ASSERT_TRUE(lhs.SetDim(10u, 10u, 10u));
    ASSERT_TRUE(lhs.SetUnitCell(10.0f, 10.0f, 10.0f, 90.0f, 90.0f, 90.0f,
                                10u, 10u, 10u));
    ASSERT_TRUE(lhs.SetMid(0.0f, 0.0f, 0.0f));

    OESystem::OESkewGrid rhs;
    ASSERT_TRUE(rhs.SetDim(10u, 10u, 10u));
    ASSERT_TRUE(rhs.SetUnitCell(10.0f, 10.0f, 10.0f, 90.0f, 90.0f, 90.0f,
                                20u, 10u, 10u));
    ASSERT_TRUE(rhs.SetMid(-2.25f, 0.0f, 0.0f));

    const GridParams a = get_grid_params(lhs);
    const GridParams b = get_grid_params(rhs);

    // Preconditions: dims match, origins match, both have unit cells with matching
    // parameters, but x spacing differs by ~0.5 A.
    ASSERT_EQ(a.x_dim, b.x_dim);
    ASSERT_EQ(a.y_dim, b.y_dim);
    ASSERT_EQ(a.z_dim, b.z_dim);
    ASSERT_DOUBLE_EQ(a.x_origin, b.x_origin);
    ASSERT_DOUBLE_EQ(a.y_origin, b.y_origin);
    ASSERT_DOUBLE_EQ(a.z_origin, b.z_origin);
    ASSERT_TRUE(lhs.HasUnitCell());
    ASSERT_TRUE(rhs.HasUnitCell());
    const UnitCellParams ca = get_unit_cell(lhs);
    const UnitCellParams cb = get_unit_cell(rhs);
    ASSERT_DOUBLE_EQ(ca.a, cb.a);
    ASSERT_DOUBLE_EQ(ca.b, cb.b);
    ASSERT_DOUBLE_EQ(ca.c, cb.c);
    ASSERT_DOUBLE_EQ(ca.alpha, cb.alpha);
    ASSERT_DOUBLE_EQ(ca.beta, cb.beta);
    ASSERT_DOUBLE_EQ(ca.gamma, cb.gamma);
    ASSERT_NEAR(std::abs(a.x_spacing - b.x_spacing), 0.5, 0.01);

    EXPECT_FALSE(same_grid_geometry(lhs, rhs));
}

TEST(SameGridGeometry, FalseWhenOnlyUnitCellPresenceDiffers) {
    // Two grids with matching derived geometry (dims, spacings, origins all equal),
    // but one has a unit cell and the other does not. This proves the HasUnitCell
    // presence check at src/Grid.cpp:193 is reachable.
    //
    // lhs: dim-only construction, measured to have no unit cell and spacings
    // 0.5/0.5/0.5, x origin -0.75.
    OESystem::OESkewGrid lhs;
    ASSERT_TRUE(lhs.SetDim(4u, 5u, 6u));

    // rhs: matching derived geometry with a unit cell.
    OESystem::OESkewGrid rhs;
    ASSERT_TRUE(rhs.SetDim(4u, 5u, 6u));
    ASSERT_TRUE(rhs.SetUnitCell(2.0f, 2.5f, 3.0f, 90.0f, 90.0f, 90.0f,
                                4u, 5u, 6u));
    ASSERT_TRUE(rhs.SetMid(0.0f, 0.0f, 0.0f));

    const GridParams a = get_grid_params(lhs);
    const GridParams b = get_grid_params(rhs);

    // Preconditions: dims, spacings, origins all match, but HasUnitCell differs.
    ASSERT_EQ(a.x_dim, b.x_dim);
    ASSERT_EQ(a.y_dim, b.y_dim);
    ASSERT_EQ(a.z_dim, b.z_dim);
    ASSERT_DOUBLE_EQ(a.x_spacing, b.x_spacing);
    ASSERT_DOUBLE_EQ(a.y_spacing, b.y_spacing);
    ASSERT_DOUBLE_EQ(a.z_spacing, b.z_spacing);
    ASSERT_DOUBLE_EQ(a.x_origin, b.x_origin);
    ASSERT_DOUBLE_EQ(a.y_origin, b.y_origin);
    ASSERT_DOUBLE_EQ(a.z_origin, b.z_origin);
    ASSERT_FALSE(lhs.HasUnitCell());
    ASSERT_TRUE(rhs.HasUnitCell());

    EXPECT_FALSE(same_grid_geometry(lhs, rhs));
}

TEST(SameGridGeometry, FalseWhenOnlyDimensionDiffers) {
    // Two grids differing only in z_dim (6 vs 7), with a compensating mid shift
    // that keeps z origin equal. This proves the dim check at src/Grid.cpp:185
    // is reachable and isolated. FalseWhenADimensionDiffers (:1018) is not
    // mutation-proof because its grids also differ in z origin.
    OESystem::OESkewGrid lhs;
    ASSERT_TRUE(lhs.SetDim(4u, 5u, 6u));
    ASSERT_TRUE(lhs.SetUnitCell(4.0f, 5.0f, 6.0f, 90.0f, 90.0f, 90.0f,
                                4u, 5u, 6u));
    ASSERT_TRUE(lhs.SetMid(0.0f, 0.0f, 0.0f));

    OESystem::OESkewGrid rhs;
    ASSERT_TRUE(rhs.SetDim(4u, 5u, 7u));
    ASSERT_TRUE(rhs.SetUnitCell(4.0f, 5.0f, 6.0f, 90.0f, 90.0f, 90.0f,
                                4u, 5u, 6u));
    ASSERT_TRUE(rhs.SetMid(0.0f, 0.0f, 0.5f));

    const GridParams a = get_grid_params(lhs);
    const GridParams b = get_grid_params(rhs);

    // Preconditions: z_dim differs (6 vs 7), but x_dim and y_dim match; all three
    // spacings match; all three origins match; both have unit cells with matching
    // parameters.
    ASSERT_NE(a.z_dim, b.z_dim);
    ASSERT_EQ(a.x_dim, b.x_dim);
    ASSERT_EQ(a.y_dim, b.y_dim);
    ASSERT_DOUBLE_EQ(a.x_spacing, b.x_spacing);
    ASSERT_DOUBLE_EQ(a.y_spacing, b.y_spacing);
    ASSERT_DOUBLE_EQ(a.z_spacing, b.z_spacing);
    ASSERT_DOUBLE_EQ(a.x_origin, b.x_origin);
    ASSERT_DOUBLE_EQ(a.y_origin, b.y_origin);
    ASSERT_DOUBLE_EQ(a.z_origin, b.z_origin);
    ASSERT_TRUE(lhs.HasUnitCell());
    ASSERT_TRUE(rhs.HasUnitCell());
    const UnitCellParams ca = get_unit_cell(lhs);
    const UnitCellParams cb = get_unit_cell(rhs);
    ASSERT_DOUBLE_EQ(ca.a, cb.a);
    ASSERT_DOUBLE_EQ(ca.b, cb.b);
    ASSERT_DOUBLE_EQ(ca.c, cb.c);
    ASSERT_DOUBLE_EQ(ca.alpha, cb.alpha);
    ASSERT_DOUBLE_EQ(ca.beta, cb.beta);
    ASSERT_DOUBLE_EQ(ca.gamma, cb.gamma);

    EXPECT_FALSE(same_grid_geometry(lhs, rhs));
}

TEST(GridParamsDerivation, DerivesTheExpectedValuesForAKnownGrid) {
    // MakeSkewGrid(4,5,6) sets cell edges equal to the dims with matching divisions,
    // so each axis's node interval is edge/divisions == 1.0 A by construction.
    // This is the independent assertion proving the derivation is correct.
    const OESystem::OESkewGrid grid = MakeSkewGrid(4u, 5u, 6u);
    const GridParams gp = get_grid_params(grid);

    // Dims match the construction.
    EXPECT_EQ(gp.x_dim, 4u);
    EXPECT_EQ(gp.y_dim, 5u);
    EXPECT_EQ(gp.z_dim, 6u);

    // Spacings are 1.0 A on all axes (independent: derived from the fixture, not
    // from reading the implementation).
    EXPECT_NEAR(gp.x_spacing, 1.0, 1e-6);
    EXPECT_NEAR(gp.y_spacing, 1.0, 1e-6);
    EXPECT_NEAR(gp.z_spacing, 1.0, 1e-6);

    // The derivation's arithmetic is self-consistent: the far coordinate on each
    // axis equals origin + (dim - 1) * spacing.
    const unsigned int n[3] = {gp.x_dim, gp.y_dim, gp.z_dim};
    const unsigned int step[3] = {1u, n[0], n[0] * n[1]};  // x-fastest linearization
    const double origin[3] = {gp.x_origin, gp.y_origin, gp.z_origin};
    const double spacing[3] = {gp.x_spacing, gp.y_spacing, gp.z_spacing};

    for (int axis = 0; axis < 3; ++axis) {
        float coord[3] = {0.0f, 0.0f, 0.0f};
        ASSERT_TRUE(grid.ElementToSpatialCoord((n[axis] - 1u) * step[axis],
                                               coord[0], coord[1], coord[2]));
        const double expected_far = origin[axis] + (n[axis] - 1u) * spacing[axis];
        EXPECT_NEAR(static_cast<double>(coord[axis]), expected_far, 1e-5);
    }

    // Origins equal what ElementToSpatialCoord(0) reports. This shares its oracle
    // with the implementation (both read element 0), so it is a consistency check
    // rather than an independent one — the spacings and span identity above are
    // the load-bearing assertions.
    float x0 = 0.0f, y0 = 0.0f, z0 = 0.0f;
    ASSERT_TRUE(grid.ElementToSpatialCoord(0u, x0, y0, z0));
    EXPECT_DOUBLE_EQ(gp.x_origin, static_cast<double>(x0));
    EXPECT_DOUBLE_EQ(gp.y_origin, static_cast<double>(y0));
    EXPECT_DOUBLE_EQ(gp.z_origin, static_cast<double>(z0));

    // Limitation: this fixture is isotropic 1.0 A on all three axes, so it cannot
    // catch an x/y/z divisor swap. Task 6's test_per_axis_spacing_matches_gemmi
    // owns that case against real anisotropic assets.
}
