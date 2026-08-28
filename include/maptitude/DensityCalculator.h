/**
 * @file DensityCalculator.h
 * @brief Electron density computation via Fourier synthesis.
 *
 * DensityCalculator computes model electron density from atomic coordinates
 * using structure factor calculation, FFT, optional bulk solvent correction,
 * and per-shell amplitude scaling.
 */

#ifndef MAPTITUDE_DENSITYCALCULATOR_H
#define MAPTITUDE_DENSITYCALCULATOR_H

#include <memory>
#include <string>
#include <vector>

#include "maptitude/SymOp.h"
#include "maptitude/UnitCell.h"

namespace OEChem {
class OEMolBase;
class OEAtomBase;
}

namespace OESystem {
class OEScalarGrid;
}

namespace OESystem {
template <class T> class OEUnaryPredicate;
}

namespace Maptitude {

/// Per-shell scaling bins partition the resolution range; even a 0.5 A dataset has
/// far fewer independent resolution shells than this. The bound also keeps
/// `n_scale_shells + 1` from wrapping the unsigned addition that sizes `shell_edges`.
constexpr unsigned int MAX_SCALE_SHELLS = 1000;

/// Miller-index generation sweeps `[-h_max, h_max] x [-k_max, k_max] x [-l_max, l_max]`
/// with `h_max = ceil(a / resolution)`, so the loop volume grows as `(2*a/resolution)^3`.
/// A finite positive resolution is not enough to bound it: 1e-9 A is an ordinary double,
/// not a subnormal, and it diverges long before `1.0 / resolution` overflows. The quantity
/// that has to be bounded is the box volume, which depends on the resolution and the three
/// cell edges jointly, so it is bounded here rather than at the resolution check.
///
/// 2e8 points leaves the crystallography this pipeline is for well inside the bound: a
/// 200 A cubic cell at 1.0 A needs 401^3 = 6.4e7 points, 3.1x under it, and a 100 A cubic
/// cell at 0.8 A needs 251^3 = 1.6e7, 12.6x under. At the bound itself the resolution
/// sphere holds pi/6 of the box, about 1.05e8 reflections, and memory rather than time is
/// what binds: `indices` is 24 bytes per reflection for 2.5 GB, `Fc_real` and `Fc_imag`
/// add 1.7 GB, and `f_s_table` adds 0.84 GB per distinct scattering type, against a few
/// seconds for the loop itself. The bound is deliberately placed where an over-large
/// request still fails as an exception rather than as an out-of-memory kill.
///
/// Declared as a double because the check has to run in arithmetic that cannot itself
/// overflow: the box volume is computed and compared in double, before any narrowing to
/// `int`.
constexpr double MAX_MILLER_BOX_POINTS = 2e8;

/**
 * @brief Computes model electron density using Fourier synthesis.
 *
 * This class encapsulates the full density computation pipeline:
 * 1. Structure preparation (extract atoms, assign radii)
 * 2. Cromer-Mann scattering factor lookup
 * 3. Miller index generation within resolution sphere
 * 4. Structure factor accumulation with symmetry expansion (OpenMP-parallel)
 * 5. Inverse FFT to real space via FFTW3
 * 6. Optional flat bulk solvent correction
 * 7. Optional per-shell amplitude scaling against observed density
 * 8. Trilinear interpolation onto output grid
 *
 * @code
 * UnitCell cell(50.0, 60.0, 70.0, 90.0, 90.0, 90.0);
 * auto symops = SymOp::ParseAll("x,y,z\n-x,y+1/2,-z+1/2");
 *
 * DensityCalculator calc(cell, symops);
 * OESystem::OEScalarGrid* result = calc.Calculate(
 *     mol, obs_grid, 2.0);  // 2.0 A resolution
 * @endcode
 */
class DensityCalculator {
public:
    /**
     * @brief Construct a density calculator for a crystal form.
     *
     * @param cell Unit cell parameters.
     * @param symops Symmetry operators (in fractional coordinates).
     */
    DensityCalculator(const UnitCell& cell, const std::vector<SymOp>& symops);

    /// @brief Destructor.
    ~DensityCalculator();

    // Non-copyable, movable
    DensityCalculator(const DensityCalculator&) = delete;
    DensityCalculator& operator=(const DensityCalculator&) = delete;
    DensityCalculator(DensityCalculator&&) noexcept;
    DensityCalculator& operator=(DensityCalculator&&) noexcept;

    /**
     * @brief Compute model electron density.
     *
     * The returned grid has the same dimensions and spacing as obs_grid.
     * The caller owns the returned pointer.
     *
     * @param mol Input molecule (OEMolBase or OEDesignUnit-derived).
     * @param obs_grid Observed electron density grid (defines output geometry).
     * @param resolution Resolution limit in Angstroms.
     * @param mask Optional atom predicate to restrict which atoms contribute.
     * @param k_sol Bulk solvent scale factor (default: 0.35 e/A^3).
     * @param b_sol Bulk solvent B-factor (default: 46.0 A^2).
     * @param include_h Include hydrogen atoms (default: false).
     * @param n_scale_shells Number of per-shell scaling bins, in
     *        [1, MAX_SCALE_SHELLS] (default: 1).
     * @return New OEScalarGrid with computed density. Caller owns the pointer.
     */
    OESystem::OEScalarGrid* Calculate(
        OEChem::OEMolBase& mol,
        const OESystem::OEScalarGrid& obs_grid,
        double resolution,
        const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask = nullptr,
        double k_sol = 0.35,
        double b_sol = 46.0,
        bool include_h = false,
        unsigned int n_scale_shells = 1) const;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

}  // namespace Maptitude

#endif  // MAPTITUDE_DENSITYCALCULATOR_H
