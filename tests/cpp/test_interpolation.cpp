/// Tier 1: the trilinear core's domain and blend (spec §3.2).
///
/// The fixture is a 4x4x4 grid with unit spacing whose node origin is the
/// Cartesian origin, so a coordinate IS its fractional index. Values hold the
/// element number, which makes every expectation checkable by hand.
#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <vector>

#include "maptitude/Grid.h"

using namespace Maptitude;

namespace {

constexpr unsigned int N = 4u;
constexpr double OUTSIDE = -99.0;

GridParams UnitGridParams() {
    return GridParams{0.0, 0.0, 0.0, N, N, N, 1.0, 1.0, 1.0};
}

/// values[el] == el, with el = iz*N*N + iy*N + ix.
std::vector<float> ElementNumberValues() {
    std::vector<float> values(N * N * N);
    for (unsigned int i = 0; i < values.size(); ++i) {
        values[i] = static_cast<float>(i);
    }
    return values;
}

}  // namespace

TEST(InterpolateDensityAt, ReturnsTheNodeValueAtANode) {
    const GridParams gp = UnitGridParams();
    const std::vector<float> v = ElementNumberValues();
    // (ix, iy, iz) = (1, 2, 3) -> el = 3*16 + 2*4 + 1 = 57.
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 1.0, 2.0, 3.0, OUTSIDE), 57.0);
}

TEST(InterpolateDensityAt, BlendsAlongEachAxisWithTheRightStride) {
    const GridParams gp = UnitGridParams();
    const std::vector<float> v = ElementNumberValues();
    // Element number is linear in the indices, so a half step along x adds 1/2,
    // along y adds 4/2, along z adds 16/2. A transposed stride moves all three.
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.5, 0.0, 0.0, OUTSIDE), 0.5);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.0, 0.5, 0.0, OUTSIDE), 2.0);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.0, 0.0, 0.5, OUTSIDE), 8.0);
}

TEST(InterpolateDensityAt, InterpolatesOnTheClosedFarFace) {
    const GridParams gp = UnitGridParams();
    const std::vector<float> v = ElementNumberValues();
    // f == n - 1 exactly. OEFloatGridLinearInterpolate returned the default
    // here; the closed interval is a deliberate departure that moves metrics.
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 3.0, 0.0, 0.0, OUTSIDE), 3.0);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 3.0, 3.0, 3.0, OUTSIDE), 63.0);
    // The face is closed, so a face-interior point bilinearly blends.
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 3.0, 0.5, 0.5, OUTSIDE), 13.0);
}

TEST(InterpolateDensityAt, ReturnsTheDefaultOneUlpOutsideEachFace) {
    const GridParams gp = UnitGridParams();
    const std::vector<float> v = ElementNumberValues();
    const double just_over = std::nextafter(3.0, 4.0);
    const double just_under = std::nextafter(0.0, -1.0);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), just_over, 0.0, 0.0, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.0, just_over, 0.0, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.0, 0.0, just_over, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), just_under, 0.0, 0.0, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.0, just_under, 0.0, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.0, 0.0, just_under, OUTSIDE), OUTSIDE);
}

TEST(InterpolateDensityAt, ReturnsTheDefaultInTheOuterHalfSpacingShell) {
    const GridParams gp = UnitGridParams();
    const std::vector<float> v = ElementNumberValues();
    // The scalar carrier's IsInGrid admitted this shell -- the old box ran to
    // 3.5 on each face. The node span does not, so the domain narrows here.
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 3.25, 1.0, 1.0, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), -0.25, 1.0, 1.0, OUTSIDE), OUTSIDE);
}

TEST(InterpolateDensityAt, ReturnsTheDefaultForANonFiniteCoordinate) {
    const GridParams gp = UnitGridParams();
    const std::vector<float> v = ElementNumberValues();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    // A NaN fractional index fails no range comparison, so ContainsFractionalIndex
    // needs its isfinite conjuncts as well as its range ones. The infinities would
    // be caught by the range test alone; the three NaN cases would not.
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), nan, 1.0, 1.0, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 1.0, nan, 1.0, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 1.0, 1.0, nan, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), inf, 1.0, 1.0, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), -inf, 1.0, 1.0, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 1.0, inf, 1.0, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 1.0, -inf, 1.0, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 1.0, 1.0, inf, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 1.0, 1.0, -inf, OUTSIDE), OUTSIDE);
}

TEST(InterpolateDensityAt, MeasuresFromTheGridOriginNotFromZero) {
    // Every other case here sits on a zero origin, which cannot distinguish
    // (x - x_origin) / spacing from x / spacing. This one can: the whole grid
    // is translated, so world (0,0,0) falls outside it entirely.
    const GridParams gp{2.0, -3.0, 10.0, N, N, N, 1.0, 1.0, 1.0};
    const std::vector<float> v = ElementNumberValues();
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 2.0, -3.0, 10.0, OUTSIDE), 0.0);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 3.0, -3.0, 10.0, OUTSIDE), 1.0);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 2.0, -2.0, 10.0, OUTSIDE), 4.0);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 2.0, -3.0, 11.0, OUTSIDE), 16.0);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 2.5, -3.0, 10.0, OUTSIDE), 0.5);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.0, 0.0, 0.0, OUTSIDE), OUTSIDE);
}

TEST(InterpolateDensityAt, UsesThePerAxisSpacingIndependently) {
    const GridParams gp{0.0, 0.0, 0.0, N, N, N, 0.5, 1.0, 2.0};
    const std::vector<float> v = ElementNumberValues();
    // One full node step on each axis: 0.5 A in x, 1.0 A in y, 2.0 A in z.
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.5, 0.0, 0.0, OUTSIDE), 1.0);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.0, 1.0, 0.0, OUTSIDE), 4.0);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.0, 0.0, 2.0, OUTSIDE), 16.0);
    // A grid that collapsed the three spacings onto one would place these
    // elsewhere; 1.0 A along x is two node steps, not one.
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 1.0, 0.0, 0.0, OUTSIDE), 2.0);
}
