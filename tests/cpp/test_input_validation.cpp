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
