#include <gtest/gtest.h>
#include "maptitude/UnitCell.h"

#include <cmath>

using namespace Maptitude;

TEST(UnitCellTest, OrthorhombicVolume) {
    UnitCell cell(50.0, 60.0, 70.0, 90.0, 90.0, 90.0);
    EXPECT_NEAR(cell.Volume(), 50.0 * 60.0 * 70.0, 1e-6);
}

TEST(UnitCellTest, CubicVolume) {
    UnitCell cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
    EXPECT_NEAR(cell.Volume(), 1000.0, 1e-6);
}

TEST(UnitCellTest, OrthorhombicFractionalRoundTrip) {
    UnitCell cell(50.0, 60.0, 70.0, 90.0, 90.0, 90.0);

    double x = 12.5, y = 30.0, z = 35.0;
    auto frac = cell.CartesianToFractional(x, y, z);
    auto cart = cell.FractionalToCartesian(frac[0], frac[1], frac[2]);

    EXPECT_NEAR(cart[0], x, 1e-10);
    EXPECT_NEAR(cart[1], y, 1e-10);
    EXPECT_NEAR(cart[2], z, 1e-10);
}

TEST(UnitCellTest, MonoclinicFractionalRoundTrip) {
    UnitCell cell(40.0, 50.0, 60.0, 90.0, 110.0, 90.0);

    double x = 10.0, y = 20.0, z = 30.0;
    auto frac = cell.CartesianToFractional(x, y, z);
    auto cart = cell.FractionalToCartesian(frac[0], frac[1], frac[2]);

    EXPECT_NEAR(cart[0], x, 1e-8);
    EXPECT_NEAR(cart[1], y, 1e-8);
    EXPECT_NEAR(cart[2], z, 1e-8);
}

TEST(UnitCellTest, ToString) {
    UnitCell cell(50.0, 60.0, 70.0, 90.0, 90.0, 90.0);
    std::string s = cell.ToString();
    EXPECT_NE(s.find("50"), std::string::npos);
    EXPECT_NE(s.find("UnitCell"), std::string::npos);
}

TEST(UnitCellTest, Equality) {
    UnitCell a(50.0, 60.0, 70.0, 90.0, 90.0, 90.0);
    UnitCell b(50.0, 60.0, 70.0, 90.0, 90.0, 90.0);
    UnitCell c(50.0, 60.0, 80.0, 90.0, 90.0, 90.0);

    EXPECT_EQ(a, b);
    EXPECT_NE(a, c);
}
