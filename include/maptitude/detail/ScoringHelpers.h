#ifndef MAPTITUDE_DETAIL_SCORINGHELPERS_H
#define MAPTITUDE_DETAIL_SCORINGHELPERS_H

#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace Maptitude {
namespace detail {

static constexpr double TWO_PI = 6.283185307179586;

// ---- Map statistics (OE-free overload for raw arrays) ----

struct MapStats {
    double mean = 0.0;
    double stddev = 0.0;
    double min_val = 0.0;
    double max_val = 0.0;
};

inline MapStats compute_map_stats(const double* values, size_t count) {
    MapStats stats;
    if (count == 0) return stats;

    double sum = 0.0;
    double sum2 = 0.0;
    stats.min_val = std::numeric_limits<double>::max();
    stats.max_val = std::numeric_limits<double>::lowest();

    for (size_t i = 0; i < count; ++i) {
        const double v = values[i];
        sum += v;
        sum2 += v * v;
        if (v < stats.min_val) stats.min_val = v;
        if (v > stats.max_val) stats.max_val = v;
    }

    stats.mean = sum / static_cast<double>(count);
    const double variance = (sum2 / static_cast<double>(count)) -
                      (stats.mean * stats.mean);
    stats.stddev = (variance > 0.0) ? std::sqrt(variance) : 0.0;
    return stats;
}

inline MapStats compute_map_stats(const std::vector<double>& values) {
    return compute_map_stats(values.data(), values.size());
}

// ---- Map normalization for Q-score (MinMaxD algorithm, Pintilie 2020) ----

inline void get_map_normalization(const double* values, size_t count,
                                  double& A, double& B) {
    const MapStats stats = compute_map_stats(values, count);
    const double max_d = std::min(stats.mean + stats.stddev * 10.0, stats.max_val);
    const double min_d = std::max(stats.mean - stats.stddev * 1.0, stats.min_val);
    A = max_d - min_d;
    B = min_d;
}

// ---- Pearson correlation ----

inline double pearson_correlation(const std::vector<double>& x,
                                  const std::vector<double>& y) {
    const size_t n = x.size();
    if (n < 2 || y.size() != n)
        return std::numeric_limits<double>::quiet_NaN();

    double sum_x = 0.0, sum_y = 0.0;
    for (size_t i = 0; i < n; ++i) {
        sum_x += x[i];
        sum_y += y[i];
    }
    const double mean_x = sum_x / static_cast<double>(n);
    const double mean_y = sum_y / static_cast<double>(n);

    double sum_xy = 0.0, sum_xx = 0.0, sum_yy = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const double dx = x[i] - mean_x;
        const double dy = y[i] - mean_y;
        sum_xy += dx * dy;
        sum_xx += dx * dx;
        sum_yy += dy * dy;
    }

    if (sum_xx <= 0.0 || sum_yy <= 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    return sum_xy / std::sqrt(sum_xx * sum_yy);
}

// ---- EDIAm sigmoid ----

inline double ediam_sigmoid(double x) {
    return 1.0 / (1.0 + std::exp(-4.0 * (x - 1.0)));
}

// ---- Adaptive scoring radius (Tickle 2012 linear model) ----

inline double scoring_radius(double bfactor, double resolution) {
    constexpr double DEFAULT_BFACTOR = 20.0;
    constexpr double MAX_RADIUS = 2.7;
    if (bfactor <= 0.0) bfactor = DEFAULT_BFACTOR;
    double dmin = std::max(resolution, 1.0);
    double r = (0.5 + 0.3 * dmin) + (0.015 - 0.002 * dmin) * bfactor;
    return std::max(0.7, std::min(r, MAX_RADIUS));
}

// ---- Binned atom radius (Phenix/CCTBX resolution-dependent bins) ----

inline double binned_atom_radius(double resolution) {
    if (resolution < 1.0)       return 1.0;
    else if (resolution < 2.0)  return 1.5;
    else if (resolution < 4.0)  return 2.0;
    else                        return 2.5;
}

// ---- Fibonacci sphere points (Pintilie 2020 SpherePts algorithm) ----

inline std::vector<std::array<double, 3>> fibonacci_sphere_points(
    double cx, double cy, double cz, double radius, int n) {
    if (n < 1) return {};
    if (n == 1) {
        return {{cx, cy + radius, cz}};
    }

    std::vector<double> thetas(n);
    std::vector<double> phis(n);

    for (int k = 1; k <= n; ++k) {
        double h = -1.0 + 2.0 * (k - 1) / (n - 1);
        if (h < -1.0) h = -1.0;
        if (h > 1.0) h = 1.0;
        phis[k - 1] = std::acos(h);
        if (k == 1 || k == n) {
            thetas[k - 1] = 0.0;
        } else {
            const double inc = 3.6 / std::sqrt(n * (1.0 - h * h));
            thetas[k - 1] = std::fmod(thetas[k - 2] + inc, TWO_PI);
        }
    }

    std::vector<std::array<double, 3>> pts(n);
    for (int i = 0; i < n; ++i) {
        const double sin_phi = std::sin(phis[i]);
        pts[i] = {
            cx + radius * sin_phi * std::cos(thetas[i]),
            cy + radius * sin_phi * std::sin(thetas[i]),
            cz + radius * std::cos(phis[i])
        };
    }
    return pts;
}

}  // namespace detail
}  // namespace Maptitude

#endif  // MAPTITUDE_DETAIL_SCORINGHELPERS_H
