#ifndef MAPTITUDE_TEST_GRID_SUMMARY_H
#define MAPTITUDE_TEST_GRID_SUMMARY_H

#include <algorithm>
#include <cmath>

#include <gtest/gtest.h>
#include <oegrid.h>

/// Shared grid-reduction helpers for characterization tests.
///
/// generate_pins.cpp must NOT include this header. It has its own
/// EmitGridSummary that duplicates the reduction logic by transcription, not by
/// sharing. That duplication is the point: one edit cannot move the pin and the
/// assertion in lockstep if they are written independently. Sharing between test
/// files is safe — neither is the generator.

namespace MaptitudeTest {

constexpr double PIN_RELATIVE_TOLERANCE = 1e-6;

struct GridSummary {
    double sum = 0.0;
    double sum_sq = 0.0;
    double min = 0.0;
    double max = 0.0;
    double index_moment = 0.0;
};

/// Reduce a grid to five scalars: sum, sum-of-squares, min, max, and a
/// mean-centred order-sensitive index moment.
///
/// The index moment is sum((i+1) * (v[i] - mean)). Mean-centring removes the
/// permutation-invariant component: the mean contributes uniformly to every
/// index, so including it only inflates the tolerance without adding signal.
/// It catches deterministic axis-order and stride permutations (e.g., a
/// transposed grid with the same value multiset).
inline GridSummary Summarize(const OESystem::OESkewGrid& grid) {
    const unsigned int size = grid.GetSize();
    const float* values = grid.GetValues();

    GridSummary s;
    s.min = values[0];
    s.max = values[0];

    // First pass: sum, sum_sq, min, max
    for (unsigned int i = 0; i < size; ++i) {
        const double v = values[i];
        s.sum += v;
        s.sum_sq += v * v;
        s.min = std::min(s.min, v);
        s.max = std::max(s.max, v);
    }

    // Second pass: mean-centred index moment
    const double mean = s.sum / static_cast<double>(size);
    for (unsigned int i = 0; i < size; ++i) {
        s.index_moment += static_cast<double>(i + 1) * (values[i] - mean);
    }

    return s;
}

/// ExpectPinned with two-regime tolerance: relative 1e-6 for |pinned| >= 1,
/// absolute 1e-6 below it.
///
/// The floor exists to avoid manufacturing cross-machine flakiness. For FC pins,
/// FFTW_ESTIMATE may pick different codelets across builds. For GridOps pins,
/// the floor ensures that near-zero results do not turn into brittle noise-floor
/// assertions. Tightening toward float epsilon would make pins fragile.
inline void ExpectPinned(double actual, double pinned) {
    const double scale = std::max(1.0, std::abs(pinned));
    EXPECT_NEAR(actual, pinned, PIN_RELATIVE_TOLERANCE * scale);
}

}  // namespace MaptitudeTest

#endif  // MAPTITUDE_TEST_GRID_SUMMARY_H
