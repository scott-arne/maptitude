/**
 * @file GridOps.h
 * @brief Grid manipulation operations (scale, combine, diff-to-calc).
 */

#ifndef MAPTITUDE_GRIDOPS_H
#define MAPTITUDE_GRIDOPS_H

namespace OESystem {
class OEScalarGrid;
}

namespace OEChem {
class OEMolBase;
}

namespace Maptitude {

/**
 * @brief Supported grid combination operations.
 */
enum class MapOp {
    ADD,       ///< Element-wise addition
    SUBTRACT,  ///< Element-wise subtraction
    MIN,       ///< Element-wise minimum
    MAX        ///< Element-wise maximum
};

/**
 * @brief Scale a grid by multiplying all values by a scalar.
 *
 * Modifies the grid in place.
 *
 * @param grid Grid to scale (modified in place).
 * @param factor Scale factor.
 */
void scale_map(OESystem::OEScalarGrid& grid, double factor);

/**
 * @brief Combine two grids element-wise.
 *
 * Both grids must have the same dimensions and spacing. The returned
 * grid has the same geometry as the input grids.
 *
 * @param lhs Left-hand side grid.
 * @param rhs Right-hand side grid.
 * @param op Combination operation.
 * @return New grid with combined values. Caller owns the pointer.
 * @throws GridError if grids have incompatible dimensions.
 */
OESystem::OEScalarGrid* combine_maps(
    const OESystem::OEScalarGrid& lhs,
    const OESystem::OEScalarGrid& rhs,
    MapOp op);

/**
 * @brief Derive calculated density from observed and difference maps.
 *
 * Computes: rho_calc = rho_obs - 2 * rho_diff
 *
 * Both grids must have the same dimensions and spacing.
 *
 * @param obs_grid Observed density map (2mFo-DFc).
 * @param diff_grid Difference density map (mFo-DFc).
 * @return New grid with calculated density. Caller owns the pointer.
 * @throws GridError if grids have incompatible dimensions.
 */
OESystem::OEScalarGrid* diff_to_calc(
    const OESystem::OEScalarGrid& obs_grid,
    const OESystem::OEScalarGrid& diff_grid);

/**
 * @brief Translate molecule into the unit cell and optionally pad the grid.
 *
 * Crystallographic CCP4 maps cover one unit cell. Deposited coordinates may
 * extend beyond the cell boundary. This function:
 * 1. Shifts all atom coordinates by integer multiples of cell vectors to
 *    bring the heavy-atom centroid near the grid center.
 * 2. If all heavy atoms (plus padding) fit within the grid, returns nullptr
 *    (caller should use the original grid).
 * 3. Otherwise, creates a new grid covering the atom range plus padding,
 *    filled by sampling the original grid with periodic wrapping.
 *
 * The molecule is modified in-place (coordinates shifted).
 *
 * @param grid CCP4 unit-cell grid.
 * @param mol Molecule to wrap (modified in-place).
 * @param cell_a Unit cell dimension a (Angstroms).
 * @param cell_b Unit cell dimension b (Angstroms).
 * @param cell_c Unit cell dimension c (Angstroms).
 * @param padding Extra margin around atoms (Angstroms).
 * @return New padded grid if needed (caller owns), or nullptr if original suffices.
 */
OESystem::OEScalarGrid* wrap_and_pad_grid(
    const OESystem::OEScalarGrid& grid,
    OEChem::OEMolBase& mol,
    double cell_a, double cell_b, double cell_c,
    double padding = 3.0);

}  // namespace Maptitude

#endif  // MAPTITUDE_GRIDOPS_H
