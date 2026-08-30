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

/// Relative tolerance wrap_and_pad_grid uses to recognise a whole number of node
/// intervals in the extent it must cover.
///
/// The node interval is measured off float node coordinates, so an extent that is
/// an exact multiple of it divides to 10.000000000000002 and rounding up alone
/// would buy a spurious node. A ratio this close to a whole count is taken as
/// that count, which means the padded node span may fall short of the requested
/// extent by up to this fraction of it. That bound is part of the function's
/// contract, not an accident: it is a hundred-thousandth of an Angstrom on a
/// ten-Angstrom extent, and no caller's padding is specified to that precision.
constexpr double PAD_INTERVAL_COUNT_TOL = 1e-6;

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
 * The molecule is modified in-place (coordinates shifted). The cell edges, the
 * padding, the source grid's geometry, the cell's commensurability and the
 * heavy-atom count are all checked before the shift, so a throw from any of
 * those leaves the molecule where it was. A shift the function applies is
 * checked the same way: every atom's shifted coordinate is computed before any
 * of them is written, so a shift whose result on some coordinate is not a
 * finite value a float can hold raises GridError with the molecule untouched
 * rather than leaving it partly moved. A centroid near enough to the grid
 * centre for the rounding to leave it there is not shifted at all, and a
 * molecule this function does not shift it does not modify. The errors raised
 * while sizing and building the padded grid
 * come after the shift, and a molecule that needed one is left where the shift
 * put it. That place is the caller's own coordinates plus whole multiples of
 * the cell vectors, up to the rounding of storing them back as floats, so the
 * molecule sits at a crystallographically equivalent position rather than a
 * corrupted one; the shift is not rolled back.
 *
 * @param grid CCP4 unit-cell grid.
 * @param mol Molecule to wrap (modified in-place).
 * @param cell_a Unit cell dimension a (Angstroms), finite and positive.
 * @param cell_b Unit cell dimension b (Angstroms), finite and positive.
 * @param cell_c Unit cell dimension c (Angstroms), finite and positive.
 * @param padding Extra margin around atoms (Angstroms), finite and
 *        non-negative. Zero is admissible; a negative value would shrink the
 *        box the atoms have to fit inside rather than widen it.
 * @return A newly allocated padded grid, or nullptr if the molecule already
 *         fits and no padding is needed. The caller owns the returned grid.
 * @throws StructureError If the molecule contains no heavy atoms.
 * @throws CellError As get_grid_params, if any cell dimension is not a finite
 *         positive value, or if a cell dimension is not the extent the grid
 *         samples on that axis -- the padded grid is filled by periodic sampling
 *         and inherits interpolate_density_periodic_at's commensurability
 *         requirement.
 * @throws GridError As get_grid_params, if padding is not a finite non-negative
 *         value, if a centroid shift the function applies would leave some atom
 *         at a coordinate that is not a finite value a float can hold, if an
 *         OESkewGrid setter rejects the padded geometry, if the
 *         atom extent plus padding is too thin on some axis to give the padded
 *         grid two nodes there, if that extent needs more node intervals on
 *         some axis than a grid dimension can hold, or if the padded grid's
 *         cell edge on some axis is larger than a float can hold.
 */
OESystem::OESkewGrid* wrap_and_pad_grid(
    const OESystem::OESkewGrid& grid,
    OEChem::OEMolBase& mol,
    double cell_a, double cell_b, double cell_c,
    double padding = 3.0);

}  // namespace Maptitude

#endif  // MAPTITUDE_GRIDOPS_H
