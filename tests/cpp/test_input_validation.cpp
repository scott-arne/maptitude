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

TEST(CellValidationTest, AcceptsAWellFormedCell) {
    EXPECT_NO_THROW(validate_cell(UnitCell(20.0, 25.0, 30.0, 90.0, 90.0, 90.0)));
    EXPECT_NO_THROW(validate_cell(UnitCell(20.0, 25.0, 30.0, 90.0, 105.0, 90.0)));
    EXPECT_NO_THROW(validate_cell(UnitCell(10.0, 10.0, 10.0, 60.0, 70.0, 80.0)));
}

TEST(CellValidationTest, RejectsNonPositiveLengths) {
    EXPECT_THROW(validate_cell(UnitCell(0.0, 25.0, 30.0, 90.0, 90.0, 90.0)), CellError);
    EXPECT_THROW(validate_cell(UnitCell(20.0, -1.0, 30.0, 90.0, 90.0, 90.0)), CellError);
    EXPECT_THROW(validate_cell(UnitCell(20.0, 25.0, 0.0, 90.0, 90.0, 90.0)), CellError);
}

TEST(CellValidationTest, RejectsNonFiniteLengths) {
    const double nan_value = std::numeric_limits<double>::quiet_NaN();
    const double inf_value = std::numeric_limits<double>::infinity();
    EXPECT_THROW(validate_cell(UnitCell(nan_value, 25.0, 30.0, 90.0, 90.0, 90.0)), CellError);
    EXPECT_THROW(validate_cell(UnitCell(20.0, inf_value, 30.0, 90.0, 90.0, 90.0)), CellError);
}

TEST(CellValidationTest, RejectsAnglesOutsideTheOpenInterval) {
    EXPECT_THROW(validate_cell(UnitCell(20.0, 25.0, 30.0, 0.0, 90.0, 90.0)), CellError);
    EXPECT_THROW(validate_cell(UnitCell(20.0, 25.0, 30.0, 180.0, 90.0, 90.0)), CellError);
    EXPECT_THROW(validate_cell(UnitCell(20.0, 25.0, 30.0, -10.0, 90.0, 90.0)), CellError);
    // cos(200) == cos(160), so the radicand test alone would let this through.
    EXPECT_THROW(validate_cell(UnitCell(20.0, 25.0, 30.0, 200.0, 90.0, 90.0)), CellError);
}

TEST(CellValidationTest, RejectsAngleTriplesWithNoRealLattice) {
    // Each angle is individually in range, but no cell has this geometry.
    EXPECT_THROW(validate_cell(UnitCell(10.0, 10.0, 10.0, 150.0, 150.0, 150.0)), CellError);
    EXPECT_THROW(validate_cell(UnitCell(10.0, 10.0, 10.0, 20.0, 20.0, 150.0)), CellError);
}

TEST(CellValidationTest, ParameterizedConstructorRejectsBadInputImmediately) {
    EXPECT_THROW(UnitCell(0.0, 25.0, 30.0, 90.0, 90.0, 90.0), CellError);
}

TEST(CellValidationTest, DefaultConstructionStillWorks) {
    // A default cell is all zeros and would fail validate_cell. That is why
    // validation lives in a free function called at consumption points rather
    // than unconditionally in every constructor.
    EXPECT_NO_THROW(UnitCell{});
}
