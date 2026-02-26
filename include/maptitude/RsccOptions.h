/**
 * @file RsccOptions.h
 * @brief Configuration options for RSCC atom radius computation.
 */

#ifndef MAPTITUDE_RSCCOPTIONS_H
#define MAPTITUDE_RSCCOPTIONS_H

namespace Maptitude {

/**
 * @brief Atom radius computation method for density scoring.
 */
enum class AtomRadius {
    FIXED,     ///< Same radius for all atoms
    SCALED,    ///< atom->GetRadius() * scaling factor
    BINNED,    ///< Resolution-dependent bin (Phenix/CCTBX)
    ADAPTIVE   ///< B-factor and resolution dependent (Tickle 2012)
};

/**
 * @brief Configuration for RSCC density scoring.
 *
 * @code
 * RsccOptions opts;
 * opts.SetAtomRadiusMethod(AtomRadius::BINNED);
 * auto result = rscc(mol, grid, resolution, nullptr, nullptr, opts);
 * @endcode
 */
class RsccOptions {
public:
    void SetAtomRadiusMethod(AtomRadius method) { method_ = method; }
    AtomRadius GetAtomRadiusMethod() const { return method_; }

    void SetFixedAtomRadius(double radius) { fixed_radius_ = radius; }
    double GetFixedAtomRadius() const { return fixed_radius_; }

    void SetAtomRadiusScaling(double scaling) { scaling_ = scaling; }
    double GetAtomRadiusScaling() const { return scaling_; }

private:
    AtomRadius method_ = AtomRadius::BINNED;
    double fixed_radius_ = 1.5;
    double scaling_ = 1.0;
};

}  // namespace Maptitude

#endif  // MAPTITUDE_RSCCOPTIONS_H
