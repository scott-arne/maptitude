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
    const OESystem::OEScalarGrid grid = MakeEmptyGrid(5.0, 1.0);

    EXPECT_THROW(wrap_and_pad_grid(grid, mol, 0.0, 5.0, 5.0), CellError);
    EXPECT_THROW(wrap_and_pad_grid(grid, mol, 5.0, 0.0, 5.0), CellError);
    EXPECT_THROW(wrap_and_pad_grid(grid, mol, 5.0, 5.0, 0.0), CellError);
    EXPECT_THROW(wrap_and_pad_grid(grid, mol, nan_value, 5.0, 5.0), CellError);
    EXPECT_THROW(wrap_and_pad_grid(grid, mol, pos_inf, 5.0, 5.0), CellError);
    EXPECT_THROW(wrap_and_pad_grid(grid, mol, -5.0, 5.0, 5.0), CellError);
}

TEST(WrapAndPadValidationTest, StillAcceptsPositiveCellEdges) {
    // The atom sits within `padding` of the grid edge, so this actually reaches the
    // fmod wrap rather than returning nullptr for "no padding needed" -- the
    // rejection test above throws before that point and cannot cover it. Cell edges
    // larger than the grid keep the centroid shift at zero, so the atom stays near
    // the edge and padding is genuinely required.
    OEChem::OEGraphMol mol = MakeAtomMol(6, 4.0, 4.0, 4.0);
    OESystem::OEScalarGrid grid = MakeEmptyGrid(5.0, 1.0);
    for (unsigned int i = 0; i < grid.GetSize(); ++i) grid[i] = 1.0f;

    std::unique_ptr<OESystem::OEScalarGrid> padded;
    ASSERT_NO_THROW(padded.reset(wrap_and_pad_grid(grid, mol, 20.0, 20.0, 20.0)));
    ASSERT_NE(padded, nullptr) << "expected the padding path, not the nullptr shortcut";
    for (unsigned int i = 0; i < padded->GetSize(); ++i) {
        ASSERT_FALSE(std::isnan((*padded)[i])) << "NaN voxel at " << i;
    }
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
