#ifndef MAPTITUDE_TEST_GRID_SUMMARY_H
#define MAPTITUDE_TEST_GRID_SUMMARY_H

#include <algorithm>
#include <cmath>

#include <oegrid.h>

/// Shared grid-reduction helpers for characterization tests.
///
/// generate_pins.cpp must NOT include this header. It has its own
/// EmitGridSummary that duplicates the reduction logic by transcription, not by
/// sharing. That duplication is the point: one edit cannot move the pin and the
/// assertion in lockstep if they are written independently. Sharing between test
/// files is safe — neither is the generator.

namespace MaptitudeTest {

constexpr double FC_RELATIVE_TOLERANCE = 1e-6;

struct GridSummary {
    double sum = 0.0;
    double sum_sq = 0.0;
    double min = 0.0;
    double max = 0.0;
    double index_moment = 0.0;
};

/// Reduce a grid to five scalars: sum, sum-of-squares, min, max, and an
/// order-sensitive index moment.
///
/// The index moment is sum((i+1) * v[i]). It catches deterministic axis-order
/// and stride permutations (e.g., a transposed grid with the same value
/// multiset). It does not catch value changes that happen to preserve the
/// weighted sum.
inline GridSummary Summarize(const OESystem::OEScalarGrid& grid) {
    GridSummary s;
    s.min = grid[0];
    s.max = grid[0];
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        const double v = grid[i];
        s.sum += v;
        s.sum_sq += v * v;
        s.min = std::min(s.min, v);
        s.max = std::max(s.max, v);
        s.index_moment += static_cast<double>(i + 1) * v;
    }
    return s;
}

/// ExpectPinned with two-regime tolerance: relative 1e-6 for |pinned| >= 1,
/// absolute 1e-6 below it.
///
/// The floor exists to avoid manufacturing cross-machine flakiness when FFTW's
/// FFTW_ESTIMATE picks different codelets. Tightening toward the float noise
/// floor on sub-1.0 values would make the pins brittle, not more precise.
inline void ExpectPinned(double actual, double pinned) {
    const double scale = std::max(1.0, std::abs(pinned));
    EXPECT_NEAR(actual, pinned, FC_RELATIVE_TOLERANCE * scale);
}

}  // namespace MaptitudeTest

#endif  // MAPTITUDE_TEST_GRID_SUMMARY_H
