/**
 * @file RSCCOptions.h
 * @brief Configuration options for RSCC atom radius computation.
 */

#ifndef MAPTITUDE_RSCCOPTIONS_H
#define MAPTITUDE_RSCCOPTIONS_H

namespace Maptitude {

/**
 * @brief Atom radius strategy for RSCC integration sphere.
 */
enum class AtomRadius {
    Fixed,    ///< Same radius for all atoms
    Scaled,   ///< atom->GetRadius() * scaling factor
    Binned    ///< Resolution-dependent bin (Phenix/CCTBX)
};

/**
 * @brief Configuration for RSCC density scoring.
 *
 * @code
 * RSCCOptions opts;
 * opts.SetAtomRadiusMethod(AtomRadius::Fixed);
 * opts.SetFixedAtomRadius(2.0);
 * auto result = RSCC(mol, grid, resolution, nullptr, nullptr, opts);
 * @endcode
 */
class RSCCOptions {
public:
    void SetAtomRadiusMethod(AtomRadius method) { method_ = method; }
    AtomRadius GetAtomRadiusMethod() const { return method_; }

    void SetFixedAtomRadius(double radius) { fixed_radius_ = radius; }
    double GetFixedAtomRadius() const { return fixed_radius_; }

    void SetAtomRadiusScaling(double scaling) { scaling_ = scaling; }
    double GetAtomRadiusScaling() const { return scaling_; }

private:
    AtomRadius method_ = AtomRadius::Binned;
    double fixed_radius_ = 1.5;
    double scaling_ = 1.0;
};

}  // namespace Maptitude

#endif  // MAPTITUDE_RSCCOPTIONS_H
