/**
 * @file CoverageOptions.h
 * @brief Configuration options for Coverage scoring.
 */

#ifndef MAPTITUDE_COVERAGEOPTIONS_H
#define MAPTITUDE_COVERAGEOPTIONS_H

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
    void SetSigma(double sigma) { sigma_ = sigma; }
    double GetSigma() const { return sigma_; }

private:
    double sigma_ = 1.0;
};

}  // namespace Maptitude

#endif  // MAPTITUDE_COVERAGEOPTIONS_H
