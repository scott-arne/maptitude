/**
 * @file RsrOptions.h
 * @brief Configuration options for RSR atom radius computation.
 */

#ifndef MAPTITUDE_RSROPTIONS_H
#define MAPTITUDE_RSROPTIONS_H

#include "maptitude/RsccOptions.h"  // for AtomRadius enum

namespace Maptitude {

/**
 * @brief Configuration for RSR density scoring.
 *
 * @code
 * RsrOptions opts;
 * opts.SetAtomRadiusMethod(AtomRadius::ADAPTIVE);
 * auto result = rsr(mol, grid, resolution, nullptr, nullptr, opts);
 * @endcode
 */
class RsrOptions {
public:
    void SetAtomRadiusMethod(AtomRadius method) { method_ = method; }
    AtomRadius GetAtomRadiusMethod() const { return method_; }

    void SetFixedAtomRadius(double radius) { fixed_radius_ = radius; }
    double GetFixedAtomRadius() const { return fixed_radius_; }

    void SetAtomRadiusScaling(double scaling) { scaling_ = scaling; }
    double GetAtomRadiusScaling() const { return scaling_; }

private:
    AtomRadius method_ = AtomRadius::ADAPTIVE;
    double fixed_radius_ = 1.5;
    double scaling_ = 1.0;
};

}  // namespace Maptitude

#endif  // MAPTITUDE_RSROPTIONS_H
