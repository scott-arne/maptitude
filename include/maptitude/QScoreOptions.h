/**
 * @file QScoreOptions.h
 * @brief Configuration options for Q-score computation.
 */

#ifndef MAPTITUDE_QSCOREOPTIONS_H
#define MAPTITUDE_QSCOREOPTIONS_H

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

namespace Maptitude {

/// The Q-score radial sweep keys each shell as `static_cast<int>(std::round(R * 1e6))`
/// (`src/Metric.cpp`). Two shells closer together than 1e-6 A collide on that key and the
/// sweep silently merges them; below roughly 1e-13 the accumulator `R += step` stops
/// advancing altogether and the loop never terminates.
constexpr double MIN_RADIAL_STEP = 1e-6;

/// The same shell key overflows `int` once `round(R * 1e6)` passes INT_MAX (2147483647),
/// wrapping to a negative key that aliases another shell. The sweep runs while
/// `R < max_radius + 0.01`, so the largest key comes from `max_radius + 0.01`:
/// (2147.47 + 0.01) * 1e6 = 2147480000, which still fits.
constexpr double MAX_RADIUS_LIMIT = 2147.47;

/// A sweep of more than this many shells is the finite-but-unbounded tail of the same
/// failure the step floor closes: at the floor a full-range sweep would be 2.1e9 shells.
/// The defaults give 4 shells; even 0.001 A sampling over 2 A gives 2000.
constexpr int MAX_SHELLS = 1000000;

/// Each shell allocates `num_points` samples per atom, and the sweep replicates
/// `num_points` centre samples besides. Published Q-score sampling uses 8
/// (Pintilie 2020); this bound is far above any real use and keeps the per-atom sample
/// vectors bounded. Values above INT_MAX additionally wrap the `static_cast<int>` at
/// `src/Metric.cpp:446`.
constexpr unsigned int MAX_NUM_POINTS = 10000;

namespace detail {
/// Reject a non-finite or non-positive option value with a message naming the
/// setter, so the caller can find it without a debugger.
inline void RequirePositiveFinite(const char* what, double value) {
    if (!std::isfinite(value) || value <= 0.0) {
        std::ostringstream message;
        message << "QScoreOptions::" << what << " requires a finite positive value (got " << value
                << ")";
        throw std::invalid_argument(message.str());
    }
}

/// Reject an option value outside its usable range, naming the setter and the bound.
inline void RequireAtMost(const char* what, double value, double limit) {
    if (value > limit) {
        std::ostringstream message;
        message << "QScoreOptions::" << what << " requires a value of at most " << limit
                << " (got " << value << ")";
        throw std::invalid_argument(message.str());
    }
}
}  // namespace detail

/**
 * @brief Radial sampling strategy for Q-score computation.
 */
enum class RadialSampling {
    FIXED,     ///< Fixed radial step and max radius
    ADAPTIVE   ///< Grid-spacing and atom-radius dependent
};

/**
 * @brief Configuration for Q-score density scoring (Pintilie et al., 2020).
 *
 * @code
 * QScoreOptions opts;
 * opts.SetSigma(0.6);
 * opts.SetNumPoints(8);
 * auto result = qscore(mol, grid, resolution, nullptr, opts);
 * @endcode
 */
class QScoreOptions {
public:
    void SetSigma(double sigma) {
        detail::RequirePositiveFinite("SetSigma", sigma);
        sigma_ = sigma;
    }
    double GetSigma() const { return sigma_; }

    void SetRadialStep(double d_rad) {
        // A non-positive step makes the radial sweep non-terminating.
        detail::RequirePositiveFinite("SetRadialStep", d_rad);
        if (d_rad < MIN_RADIAL_STEP) {
            std::ostringstream message;
            message << "QScoreOptions::SetRadialStep requires a value of at least "
                    << MIN_RADIAL_STEP << " (MIN_RADIAL_STEP, got " << d_rad
                    << "); the radial sweep's shell key has 1e-6 A resolution, so smaller "
                       "steps alias distinct shells onto one key";
            throw std::invalid_argument(message.str());
        }
        d_rad_ = d_rad;
    }
    double GetRadialStep() const { return d_rad_; }

    void SetMaxRadius(double to_rad) {
        detail::RequirePositiveFinite("SetMaxRadius", to_rad);
        detail::RequireAtMost("SetMaxRadius", to_rad, MAX_RADIUS_LIMIT);
        to_rad_ = to_rad;
    }
    double GetMaxRadius() const { return to_rad_; }

    void SetNumPoints(unsigned int num_points) {
        if (num_points == 0) {
            // num_points is unsigned, so zero is the only reachable bad value; name
            // it anyway, per the "errors carry their values" rule in spec 5.3.
            throw std::invalid_argument("QScoreOptions::SetNumPoints requires at least one point, got " +
                                        std::to_string(num_points));
        }
        if (num_points > MAX_NUM_POINTS) {
            throw std::invalid_argument("QScoreOptions::SetNumPoints requires at most " +
                                        std::to_string(MAX_NUM_POINTS) +
                                        " points (MAX_NUM_POINTS), got " +
                                        std::to_string(num_points));
        }
        num_points_ = num_points;
    }
    unsigned int GetNumPoints() const { return num_points_; }

    void SetNormalizeMap(bool normalize) { normalize_map_ = normalize; }
    bool GetNormalizeMap() const { return normalize_map_; }

    void SetIsolatePoints(bool isolate) { isolate_points_ = isolate; }
    bool GetIsolatePoints() const { return isolate_points_; }

    void SetRadialSampling(RadialSampling method) { radial_sampling_ = method; }
    RadialSampling GetRadialSampling() const { return radial_sampling_; }

private:
    double sigma_ = 0.6;                              ///< Gaussian width in Angstroms
    double d_rad_ = 0.5;                              ///< Radial step size in Angstroms
    double to_rad_ = 2.0;                             ///< Maximum sampling radius in Angstroms
    unsigned int num_points_ = 8;                     ///< Points per radial shell
    bool normalize_map_ = true;                       ///< Normalize map before scoring
    bool isolate_points_ = true;                      ///< Exclude shell points near neighbor atoms
    RadialSampling radial_sampling_ = RadialSampling::FIXED;  ///< Radial sampling strategy
};

}  // namespace Maptitude

#endif  // MAPTITUDE_QSCOREOPTIONS_H
