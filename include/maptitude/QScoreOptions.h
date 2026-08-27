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
        d_rad_ = d_rad;
    }
    double GetRadialStep() const { return d_rad_; }

    void SetMaxRadius(double to_rad) {
        detail::RequirePositiveFinite("SetMaxRadius", to_rad);
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
