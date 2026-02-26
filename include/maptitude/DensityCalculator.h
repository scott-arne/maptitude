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
     * @param n_scale_shells Number of per-shell scaling bins (default: 1).
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
