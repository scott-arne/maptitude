/// Rejection tests for every validation throw site added in Phase 1.
///
/// Named test_input_validation to avoid confusion with the unrelated
/// tests/python/test_validation.py, which validates scores against reference
/// data.
#include <gtest/gtest.h>
#include <limits>

#include "maptitude/Error.h"
#include "maptitude/UnitCell.h"

using namespace Maptitude;

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
