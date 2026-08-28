/**
 * @file RsccOptions.h
 * @brief Configuration options for RSCC atom radius computation.
 */

#ifndef MAPTITUDE_RSCCOPTIONS_H
#define MAPTITUDE_RSCCOPTIONS_H

#include <sstream>
#include <stdexcept>

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

namespace detail {
/// Reject an `AtomRadius` value the enum does not declare, naming the setter.
///
/// `AtomRadius` is a scoped enum with underlying type `int`, so
/// `static_cast<AtomRadius>(42)` is a valid value of the type and SWIG passes one
/// through from `SetAtomRadiusMethod(42)` without a cast. It matches no arm of the
/// radius switches in `rscc` and `rsr`, which would leave their `radius` indeterminate.
/// Validating at the setter fails at the point of the mistake and is what makes those
/// switches exhaustive over the values the type can hold rather than only over the ones
/// it declares.
///
/// The switch below deliberately carries no `default:` label, so a new enumerator is a
/// -Wswitch warning here as well; the throw sits after it, where only an undeclared
/// value can arrive.
inline void RequireDeclaredAtomRadius(const char* what, AtomRadius method) {
    switch (method) {
        case AtomRadius::FIXED:
        case AtomRadius::SCALED:
        case AtomRadius::BINNED:
        case AtomRadius::ADAPTIVE:
            return;
    }
    std::ostringstream message;
    message << what << " requires a declared AtomRadius value (got " << static_cast<int>(method)
            << ")";
    throw std::invalid_argument(message.str());
}
}  // namespace detail

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
    void SetAtomRadiusMethod(AtomRadius method) {
        detail::RequireDeclaredAtomRadius("RsccOptions::SetAtomRadiusMethod", method);
        method_ = method;
    }
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
