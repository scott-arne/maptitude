/**
 * @file GridOps.h
 * @brief Grid manipulation operations (scale, combine, diff-to-calc).
 */

#ifndef MAPTITUDE_GRIDOPS_H
#define MAPTITUDE_GRIDOPS_H

namespace OESystem {
class OESkewGrid;
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
void scale_map(OESystem::OESkewGrid& grid, double factor);

/**
 * @brief Combine two grids element-wise.
 *
 * Both grids must have identical geometry (dimensions, node
 * origin, and per-axis spacing) as determined by same_grid_geometry. The
 * returned grid has the same geometry as the input grids.
 *
 * @param lhs Left-hand side grid.
 * @param rhs Right-hand side grid.
 * @param op Combination operation.
 * @return New grid with combined values. Caller owns the pointer.
 * @throws GridError if grids have different geometry.
 * @throws GridError, CellError If either grid's geometry cannot be derived.
 *         same_grid_geometry raises here where OEGridSameGeometry returned
 *         false for a grid it could not compare.
 */
OESystem::OESkewGrid* combine_maps(
    const OESystem::OESkewGrid& lhs,
    const OESystem::OESkewGrid& rhs,
    MapOp op);

/**
 * @brief Derive calculated density from observed and difference maps.
 *
 * Computes: rho_calc = rho_obs - 2 * rho_diff
 *
 * Both grids must have identical geometry (dimensions, node
 * origin, and per-axis spacing) as determined by same_grid_geometry. The
 * returned grid has the same geometry as the input grids.
 *
 * @param obs_grid Observed density map (2mFo-DFc).
 * @param diff_grid Difference density map (mFo-DFc).
 * @return New grid with calculated density. Caller owns the pointer.
 * @throws GridError if grids have different geometry.
 * @throws GridError, CellError If either grid's geometry cannot be derived.
 *         same_grid_geometry raises here where OEGridSameGeometry returned
 *         false for a grid it could not compare.
 */
OESystem::OESkewGrid* diff_to_calc(
    const OESystem::OESkewGrid& obs_grid,
    const OESystem::OESkewGrid& diff_grid);

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
 * @param cell_a Unit cell dimension a (Angstroms), finite and positive.
 * @param cell_b Unit cell dimension b (Angstroms), finite and positive.
 * @param cell_c Unit cell dimension c (Angstroms), finite and positive.
 * @param padding Extra margin around atoms (Angstroms).
 * @return A newly allocated padded grid, or nullptr if the molecule already
 *         fits and no padding is needed. The caller owns the returned grid.
 * @throws StructureError If the molecule contains no heavy atoms.
 * @throws CellError If any cell dimension is not a finite positive value, or if
 *         a cell dimension is not the extent the grid samples on that axis --
 *         the padded grid is filled by periodic sampling and inherits
 *         interpolate_density_periodic_at's commensurability requirement.
 * @throws GridError As get_grid_params, if an OESkewGrid setter rejects the
 *         padded geometry, or if the atom extent plus padding is too thin on
 *         some axis to give the padded grid two nodes there.
 */
OESystem::OESkewGrid* wrap_and_pad_grid(
    const OESystem::OESkewGrid& grid,
    OEChem::OEMolBase& mol,
    double cell_a, double cell_b, double cell_c,
    double padding = 3.0);

}  // namespace Maptitude

#endif  // MAPTITUDE_GRIDOPS_H
