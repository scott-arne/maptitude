/**
 * @file CoverageOptions.h
 * @brief Configuration options for Coverage scoring.
 */

#ifndef MAPTITUDE_COVERAGEOPTIONS_H
#define MAPTITUDE_COVERAGEOPTIONS_H

#include <cmath>
#include <sstream>
#include <stdexcept>

namespace Maptitude {

/**
 * @brief Configuration for Coverage density scoring.
 *
 * @code
 * CoverageOptions opts;
 * opts.SetSigma(1.5);
 * auto result = coverage(mol, grid, nullptr, opts);
 * @endcode
 */
class CoverageOptions {
public:
    /// Reject only a non-finite sigma.
    ///
    /// Deliberately weaker than `QScoreOptions::SetSigma`, which also refuses zero and
    /// negatives. The two sigmas are different quantities: Q-score's is a Gaussian width,
    /// which has to be positive, while coverage's is a multiplier in the threshold
    /// `mean + sigma * stddev`. There zero means "threshold at the mean" and a negative
    /// value means "threshold below the mean" -- both meaningful requests, so the
    /// asymmetry is intended and must not be "fixed". NaN is not: it makes every
    /// `rho >= threshold` comparison false and coverage returns a plausible 0.0.
    void SetSigma(double sigma) {
        if (!std::isfinite(sigma)) {
            std::ostringstream message;
            message << "CoverageOptions::SetSigma requires a finite value (got " << sigma << ")";
            throw std::invalid_argument(message.str());
        }
        sigma_ = sigma;
    }
    double GetSigma() const { return sigma_; }

private:
    double sigma_ = 1.0;
};

}  // namespace Maptitude

#endif  // MAPTITUDE_COVERAGEOPTIONS_H
