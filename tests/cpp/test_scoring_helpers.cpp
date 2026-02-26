#include <gtest/gtest.h>
#include "maptitude/detail/ScoringHelpers.h"

#include <cmath>
#include <numeric>
#include <set>

using namespace Maptitude::detail;

// =====================================================================
// PearsonCorrelation
// =====================================================================

TEST(PearsonCorrelationTest, PerfectPositiveCorrelation) {
    std::vector<double> x = {1, 2, 3, 4, 5};
    std::vector<double> y = {2, 4, 6, 8, 10};
    EXPECT_NEAR(pearson_correlation(x, y), 1.0, 1e-12);
}

TEST(PearsonCorrelationTest, PerfectNegativeCorrelation) {
    std::vector<double> x = {1, 2, 3, 4, 5};
    std::vector<double> y = {10, 8, 6, 4, 2};
    EXPECT_NEAR(pearson_correlation(x, y), -1.0, 1e-12);
}

TEST(PearsonCorrelationTest, ZeroCorrelation) {
    // Orthogonal signals
    std::vector<double> x = {1, 0, -1, 0};
    std::vector<double> y = {0, 1, 0, -1};
    EXPECT_NEAR(pearson_correlation(x, y), 0.0, 1e-12);
}

TEST(PearsonCorrelationTest, KnownValues) {
    // np.corrcoef([1,2,3,4,5], [2,4,5,4,5])[0,1] = 0.7745966692414834
    std::vector<double> x = {1, 2, 3, 4, 5};
    std::vector<double> y = {2, 4, 5, 4, 5};
    EXPECT_NEAR(pearson_correlation(x, y), 0.7745966692414834, 1e-10);
}

TEST(PearsonCorrelationTest, SelfCorrelation) {
    std::vector<double> x = {3.14, 2.72, 1.41, 1.62, 0.58};
    EXPECT_NEAR(pearson_correlation(x, x), 1.0, 1e-12);
}

TEST(PearsonCorrelationTest, ConstantInputReturnsNaN) {
    std::vector<double> x = {5, 5, 5, 5};
    std::vector<double> y = {1, 2, 3, 4};
    EXPECT_TRUE(std::isnan(pearson_correlation(x, y)));
}

TEST(PearsonCorrelationTest, SingleElementReturnsNaN) {
    std::vector<double> x = {1};
    std::vector<double> y = {2};
    EXPECT_TRUE(std::isnan(pearson_correlation(x, y)));
}

TEST(PearsonCorrelationTest, EmptyReturnsNaN) {
    std::vector<double> x, y;
    EXPECT_TRUE(std::isnan(pearson_correlation(x, y)));
}

TEST(PearsonCorrelationTest, MismatchedSizeReturnsNaN) {
    std::vector<double> x = {1, 2, 3};
    std::vector<double> y = {1, 2};
    EXPECT_TRUE(std::isnan(pearson_correlation(x, y)));
}

TEST(PearsonCorrelationTest, ScaleInvariant) {
    std::vector<double> x = {1, 2, 3, 4, 5};
    std::vector<double> y = {2, 4, 5, 4, 5};
    // Scale x by 100 and shift by 1000 — should not change correlation
    std::vector<double> x_scaled(x.size());
    for (size_t i = 0; i < x.size(); ++i) x_scaled[i] = x[i] * 100 + 1000;
    EXPECT_NEAR(pearson_correlation(x, y),
                pearson_correlation(x_scaled, y), 1e-10);
}

// =====================================================================
// FibonacciSpherePoints
// =====================================================================

TEST(FibonacciSphereTest, CorrectPointCount) {
    for (int n : {1, 2, 4, 8, 16, 32, 100}) {
        auto pts = fibonacci_sphere_points(0, 0, 0, 1.0, n);
        EXPECT_EQ(static_cast<int>(pts.size()), n) << "n=" << n;
    }
}

TEST(FibonacciSphereTest, ZeroPointsReturnsEmpty) {
    auto pts = fibonacci_sphere_points(0, 0, 0, 1.0, 0);
    EXPECT_TRUE(pts.empty());
}

TEST(FibonacciSphereTest, NegativePointsReturnsEmpty) {
    auto pts = fibonacci_sphere_points(0, 0, 0, 1.0, -5);
    EXPECT_TRUE(pts.empty());
}

TEST(FibonacciSphereTest, AllPointsAtCorrectRadius) {
    double r = 2.5;
    auto pts = fibonacci_sphere_points(0, 0, 0, r, 100);
    for (const auto& p : pts) {
        double dist = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
        EXPECT_NEAR(dist, r, 1e-10);
    }
}

TEST(FibonacciSphereTest, CenterOffsetApplied) {
    double cx = 10.0, cy = 20.0, cz = 30.0;
    double r = 1.0;
    auto pts = fibonacci_sphere_points(cx, cy, cz, r, 50);
    for (const auto& p : pts) {
        double dx = p[0] - cx;
        double dy = p[1] - cy;
        double dz = p[2] - cz;
        double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
        EXPECT_NEAR(dist, r, 1e-10);
    }
}

TEST(FibonacciSphereTest, SinglePointIsNorthPole) {
    auto pts = fibonacci_sphere_points(5.0, 10.0, 15.0, 3.0, 1);
    ASSERT_EQ(pts.size(), 1u);
    EXPECT_NEAR(pts[0][0], 5.0, 1e-12);
    EXPECT_NEAR(pts[0][1], 13.0, 1e-12);  // cy + radius
    EXPECT_NEAR(pts[0][2], 15.0, 1e-12);
}

TEST(FibonacciSphereTest, PointsAreDistributed) {
    // With enough points, the centroid should be near the center
    auto pts = fibonacci_sphere_points(0, 0, 0, 1.0, 1000);
    double sx = 0, sy = 0, sz = 0;
    for (const auto& p : pts) {
        sx += p[0]; sy += p[1]; sz += p[2];
    }
    sx /= pts.size(); sy /= pts.size(); sz /= pts.size();
    EXPECT_NEAR(sx, 0.0, 0.05);
    EXPECT_NEAR(sy, 0.0, 0.05);
    EXPECT_NEAR(sz, 0.0, 0.05);
}

TEST(FibonacciSphereTest, FirstAndLastPointsArePoles) {
    auto pts = fibonacci_sphere_points(0, 0, 0, 1.0, 50);
    // First point should be south pole (z ≈ -1)
    EXPECT_NEAR(pts[0][2], -1.0, 1e-10);
    // Last point should be north pole (z ≈ +1)
    EXPECT_NEAR(pts.back()[2], 1.0, 1e-10);
}

TEST(FibonacciSphereTest, ConsistencyWithPython) {
    // Verify first few points match the Python _nb_fibonacci_sphere_points
    // for n=8, r=0.5 centered at origin
    auto pts = fibonacci_sphere_points(0, 0, 0, 0.5, 8);
    ASSERT_EQ(pts.size(), 8u);

    // First point: k=1, h=-1, phi=pi, theta=0 -> (0, 0, -0.5)
    EXPECT_NEAR(pts[0][0], 0.0, 1e-10);
    EXPECT_NEAR(pts[0][1], 0.0, 1e-10);
    EXPECT_NEAR(pts[0][2], -0.5, 1e-10);

    // Last point: k=8, h=1, phi=0, theta=0 -> (0, 0, 0.5)
    EXPECT_NEAR(pts[7][0], 0.0, 1e-10);
    EXPECT_NEAR(pts[7][1], 0.0, 1e-10);
    EXPECT_NEAR(pts[7][2], 0.5, 1e-10);
}

// =====================================================================
// EDIAmSigmoid
// =====================================================================

TEST(EDIAmSigmoidTest, AtOneReturnsHalf) {
    EXPECT_NEAR(ediam_sigmoid(1.0), 0.5, 1e-12);
}

TEST(EDIAmSigmoidTest, MonotonicallyIncreasing) {
    double prev = ediam_sigmoid(-10.0);
    for (double x = -9.0; x <= 10.0; x += 0.5) {
        double val = ediam_sigmoid(x);
        EXPECT_GE(val, prev) << "x=" << x;
        prev = val;
    }
}

TEST(EDIAmSigmoidTest, BoundedZeroToOne) {
    for (double x = -100.0; x <= 100.0; x += 1.0) {
        double val = ediam_sigmoid(x);
        EXPECT_GE(val, 0.0) << "x=" << x;
        EXPECT_LE(val, 1.0) << "x=" << x;
    }
}

TEST(EDIAmSigmoidTest, SaturatesHigh) {
    EXPECT_NEAR(ediam_sigmoid(5.0), 1.0, 1e-6);
}

TEST(EDIAmSigmoidTest, SaturatesLow) {
    EXPECT_NEAR(ediam_sigmoid(-3.0), 0.0, 1e-6);
}

TEST(EDIAmSigmoidTest, KnownValue) {
    // At x=1.5: 1/(1+exp(-4*(1.5-1))) = 1/(1+exp(-2))
    double expected = 1.0 / (1.0 + std::exp(-2.0));
    EXPECT_NEAR(ediam_sigmoid(1.5), expected, 1e-12);
}

// =====================================================================
// ScoringRadius
// =====================================================================

TEST(ScoringRadiusTest, DefaultBFactor) {
    // With default B-factor (20) and resolution 2.0:
    // dmin = 2.0, r = (0.5 + 0.3*2) + (0.015 - 0.002*2)*20 = 1.1 + 0.22 = 1.32
    double r = scoring_radius(20.0, 2.0);
    EXPECT_NEAR(r, 1.32, 1e-10);
}

TEST(ScoringRadiusTest, HighBFactorClamped) {
    // Very high B-factor should be clamped to MAX_RADIUS = 2.7
    double r = scoring_radius(200.0, 2.0);
    EXPECT_NEAR(r, 2.7, 1e-10);
}

TEST(ScoringRadiusTest, ZeroBFactorUsesDefault) {
    // B-factor <= 0 should use DEFAULT_BFACTOR = 20
    double r_zero = scoring_radius(0.0, 2.0);
    double r_default = scoring_radius(20.0, 2.0);
    EXPECT_NEAR(r_zero, r_default, 1e-10);
}

TEST(ScoringRadiusTest, NegativeBFactorUsesDefault) {
    double r_neg = scoring_radius(-5.0, 2.0);
    double r_default = scoring_radius(20.0, 2.0);
    EXPECT_NEAR(r_neg, r_default, 1e-10);
}

TEST(ScoringRadiusTest, LowResolutionClamped) {
    // Resolution < 1.0 should be clamped to 1.0
    double r = scoring_radius(20.0, 0.5);
    double r_clamped = scoring_radius(20.0, 1.0);
    EXPECT_NEAR(r, r_clamped, 1e-10);
}

TEST(ScoringRadiusTest, MinimumFloor) {
    // Result should never be below 0.7
    double r = scoring_radius(1.0, 1.0);
    EXPECT_GE(r, 0.7);
}

TEST(ScoringRadiusTest, IncreasingWithBFactor) {
    double r_low = scoring_radius(10.0, 2.0);
    double r_high = scoring_radius(50.0, 2.0);
    EXPECT_GT(r_high, r_low);
}

// =====================================================================
// ComputeMapStats
// =====================================================================

TEST(MapStatsTest, EmptyReturnsZeros) {
    auto stats = compute_map_stats(static_cast<const double*>(nullptr), 0);
    EXPECT_DOUBLE_EQ(stats.mean, 0.0);
    EXPECT_DOUBLE_EQ(stats.stddev, 0.0);
}

TEST(MapStatsTest, SingleValue) {
    std::vector<double> v = {5.0};
    auto stats = compute_map_stats(v);
    EXPECT_NEAR(stats.mean, 5.0, 1e-12);
    EXPECT_NEAR(stats.stddev, 0.0, 1e-12);
    EXPECT_NEAR(stats.min_val, 5.0, 1e-12);
    EXPECT_NEAR(stats.max_val, 5.0, 1e-12);
}

TEST(MapStatsTest, KnownDistribution) {
    // Values: 1, 2, 3, 4, 5
    // Mean: 3.0, Var: 2.0, Std: sqrt(2) ≈ 1.4142
    std::vector<double> v = {1, 2, 3, 4, 5};
    auto stats = compute_map_stats(v);
    EXPECT_NEAR(stats.mean, 3.0, 1e-12);
    EXPECT_NEAR(stats.stddev, std::sqrt(2.0), 1e-12);
    EXPECT_NEAR(stats.min_val, 1.0, 1e-12);
    EXPECT_NEAR(stats.max_val, 5.0, 1e-12);
}

TEST(MapStatsTest, NegativeValues) {
    std::vector<double> v = {-3, -1, 0, 1, 3};
    auto stats = compute_map_stats(v);
    EXPECT_NEAR(stats.mean, 0.0, 1e-12);
    EXPECT_NEAR(stats.min_val, -3.0, 1e-12);
    EXPECT_NEAR(stats.max_val, 3.0, 1e-12);
}

TEST(MapStatsTest, ConstantArray) {
    std::vector<double> v(100, 7.5);
    auto stats = compute_map_stats(v);
    EXPECT_NEAR(stats.mean, 7.5, 1e-12);
    EXPECT_NEAR(stats.stddev, 0.0, 1e-12);
}

// =====================================================================
// GetMapNormalization
// =====================================================================

TEST(MapNormalizationTest, KnownValues) {
    // Values 0..9, mean=4.5, std≈2.872
    std::vector<double> v = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    double A, B;
    get_map_normalization(v.data(), v.size(), A, B);
    // min_d = max(mean - 1*std, min) = max(4.5 - 2.872, 0) = 1.628
    // max_d = min(mean + 10*std, max) = min(4.5 + 28.72, 9) = 9
    // A = 9 - 1.628 = 7.372, B = 1.628
    EXPECT_GT(A, 0.0);
    EXPECT_GT(B, 0.0);
    EXPECT_NEAR(A + B, 9.0, 0.01);  // max_d clamped to 9
}
