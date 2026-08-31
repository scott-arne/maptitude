/**
 * @file Grid.h
 * @brief Utilities for working with OESkewGrid objects.
 *
 * Provides helper functions for converting between OpenEye grid
 * representations and raw data arrays used by the computation kernels.
 */

#ifndef MAPTITUDE_GRID_H
#define MAPTITUDE_GRID_H

#include <cstddef>
#include <vector>

namespace OESystem {
class OESkewGrid;
}

namespace Maptitude {

/**
 * @brief Parameters describing a grid's geometry.
 *
 * Spacing is per-axis: crystallographic maps commonly sample the three cell
 * edges at different intervals, and collapsing them to one scalar both
 * resamples the density and mispositions the nodes.
 */
struct GridParams {
    double x_origin;      ///< Cartesian x of element 0 (Angstroms)
    double y_origin;      ///< Cartesian y of element 0 (Angstroms)
    double z_origin;      ///< Cartesian z of element 0 (Angstroms)
    unsigned int x_dim;   ///< Number of grid nodes along x
    unsigned int y_dim;   ///< Number of grid nodes along y
    unsigned int z_dim;   ///< Number of grid nodes along z
    double x_spacing;     ///< Node interval along x (Angstroms)
    double y_spacing;     ///< Node interval along y (Angstroms)
    double z_spacing;     ///< Node interval along z (Angstroms)
};

/**
 * @brief Per-axis wrap period, in node intervals, that a cell selects on a grid.
 *
 * The periodic path wraps in index space, so the period it needs is a count of
 * node intervals rather than a Cartesian length. require_commensurate_cell
 * derives it from the caller's cell edge and returns it here so the wrap does
 * not have to guess it from the node count.
 */
struct CellPeriods {
    unsigned int x_period;  ///< Node intervals in one period along x
    unsigned int y_period;  ///< Node intervals in one period along y
    unsigned int z_period;  ///< Node intervals in one period along z
};

/// Unit cell of a skew grid: a, b, c in Angstroms, alpha, beta, gamma in degrees.
struct UnitCellParams {
    double a;      ///< Cell edge a (Angstroms)
    double b;      ///< Cell edge b (Angstroms)
    double c;      ///< Cell edge c (Angstroms)
    double alpha;  ///< Angle between b and c (degrees)
    double beta;   ///< Angle between a and c (degrees)
    double gamma;  ///< Angle between a and b (degrees)
};

/**
 * @brief Derive per-axis geometry from a skew grid.
 *
 * OESkewGrid exposes no spacing getter that survives anisotropy, so the
 * intervals are measured: walk from element 0 to the far node of each axis
 * with ElementToSpatialCoord and divide by the interval count.
 *
 * @param grid Input grid.
 * @return GridParams describing the grid geometry.
 * @throws GridError If an axis has fewer than two nodes, a node has no spatial
 *         coordinate or a non-finite one, or a derived interval is not finite
 *         and positive.
 * @throws CellError If an axis leaks more than 1e-4 Angstroms into another
 *         axis over its full span, i.e. the sampling is not axis-aligned.
 */
GridParams get_grid_params(const OESystem::OESkewGrid& grid);

/**
 * @brief Read a skew grid's unit cell.
 *
 * @param grid Input grid.
 * @return The six cell parameters.
 * @throws CellError If the grid has no unit cell. Returning zeros would let a
 *         caller divide by an edge that was never set.
 */
UnitCellParams get_unit_cell(const OESystem::OESkewGrid& grid);

/// Cartesian position of the first node (element 0).
void grid_node_origin(const GridParams& gp, double& x, double& y, double& z);

/**
 * @brief Fractional node index along each axis.
 *
 * Values outside [0, n_i - 1] denote a point outside the grid.
 */
void grid_fractional_index(const GridParams& gp,
                           double x, double y, double z,
                           double& fx, double& fy, double& fz);

/**
 * @brief True when (x, y, z) lies within the node span on all three axes.
 *
 * That is the domain over which interpolate_density returns interpolated data.
 * It is narrower than the scalar carrier's IsInGrid, which admitted the
 * half-spacing shell outside the outermost nodes.
 */
bool grid_contains(const GridParams& gp, double x, double y, double z);

/**
 * @brief Trilinear interpolation at a Cartesian point, from derived geometry.
 *
 * The caller derives @p gp once and hoists @p values out of its loop. Both are
 * required: deriving geometry per point costs four ElementToSpatialCoord calls
 * and five validation checks, and calling OESkewGrid::GetValues per point
 * re-enters the OpenEye shared library on every element.
 *
 * The node span is closed on both ends, with a small boundary tolerance. A
 * point exactly on the far face interpolates rather than falling out of the
 * grid, and so does one a rounding error beyond it: the span's endpoints are
 * derived from float node coordinates and are not exact. The tolerance is
 * proportional to the largest magnitude the axis's geometry takes -- the further
 * of its two endpoints, or the sampled extent n_i * spacing_i, whichever is
 * greater -- because the error it absorbs is a float rounding of a quantity at
 * that magnitude. Which of those two dominates decides whether the tolerance
 * tracks the grid's distance from the Cartesian origin. While no endpoint is
 * further out than the extent it does not: a five-node unit-spacing grid
 * starting at 0 and the same grid centred on the origin get the identical
 * tolerance. Once an endpoint is further out, it does: the same grid starting
 * at 3000 A gets a tolerance about 600x wider.
 *
 * @p gp must come from get_grid_params, and that is what makes the base-index
 * clamp safe: check 1 there rejects any axis with fewer than two nodes, so the
 * unsigned `n_i - 2` the clamp computes cannot wrap.
 *
 * @param gp Geometry from get_grid_params.
 * @param values The grid's value array, from OESkewGrid::GetValues().
 * @param x Cartesian x coordinate.
 * @param y Cartesian y coordinate.
 * @param z Cartesian z coordinate.
 * @param default_value Value returned for a point outside the node span, or
 *        for a non-finite coordinate.
 * @return Interpolated density value.
 */
double interpolate_density_at(const GridParams& gp, const float* values,
                              double x, double y, double z,
                              double default_value = 0.0);

/**
 * @brief Throw unless each cell edge is a whole number of node intervals, and
 *        that number is the node count or one less.
 *
 * The periodic path wraps in index space, so what it needs from the cell is an
 * interval count: `p_i = round(cell_i / spacing_i)`. Two counts are accepted,
 * because two grid shapes reach this library.
 *
 * A grid can hold `n_i` distinct nodes per period, and then node `n_i - 1`'s
 * upper neighbour is node 0. The four whole-cell crystallographic assets in
 * `tests/assets/mapq` are not that shape. Each stores exactly its `MX,MY,MZ`
 * sections at `NxSTART = 0`, OpenEye's reader appends a closing plane to each,
 * and all four come back with M + 1 nodes per axis -- so one period spans
 * `n_i - 1` intervals. Measured on all four, that appended plane is a bit-exact
 * copy of plane 0 on every axis: 0 differing values, out of 1225 to 7380 per
 * plane. All four report a cell within 8e-06 A of `(n_i - 1) * spacing_i` on all
 * three axes, which is the interval count `llround` selects on each.
 * The reader does not do this unconditionally, and `tests/data/test_map.ccp4`
 * is the case that shows so: it stores 21 of its cell's 42 samples at
 * `NxSTART = -10` and 21 nodes come back, with no period at all.
 *
 * **The caller's cell selects the period; this function does not divine it.**
 * Both counts pass the arithmetic on any grid with two or more nodes per axis,
 * and which of them describes a particular grid turns on whether its closing
 * plane repeats its first -- a property of the values, which the geometry does
 * not record. Passing an edge of `(n_i - 1) * spacing_i` is the caller asserting
 * that it does, exactly as passing `n_i * spacing_i` asserts that every node is
 * distinct. The guard's job is to reject a cell that matches no node lattice at
 * all, and it does not read the values to decide which of the two the caller
 * meant: comparing the end planes costs O(n^2) per call on the hot path, and the
 * batch entry point would pay it for a property of the grid rather than of the
 * batch.
 *
 * An edge that rounds to neither count describes a different lattice, and
 * wrapping onto it would return densities from the wrong place with nothing to
 * mark them as wrong.
 *
 * The comparison is absolute, against an allowance proportional to the largest
 * magnitude the axis's geometry takes -- the further of its two endpoints, or
 * n_i * spacing_i itself, whichever is greater. Both sides of that max matter:
 * the extent is derived from float node coordinates, whose error grows with
 * their magnitude, and on a grid centred on the Cartesian origin the extent is
 * always the larger of the two -- by 2n_i / (n_i - 1), which is fourfold on a
 * two-node axis. Leaving it out of the scale is enough to make a two-node grid
 * at 5.45 A centred on the origin fail on a cell it tiles exactly. A cell that
 * matches neither count for a real reason misses the nearer one by a fraction
 * of a node interval at least, which is orders above the allowance for any grid
 * a crystallographic map produces.
 *
 * Exposed so a caller that is about to sample the grid periodically can reject
 * a bad cell before doing any other work, rather than after. The interval counts
 * come back so the wrap does not recompute them.
 *
 * @param gp Geometry from get_grid_params.
 * @param cell_a Unit cell dimension along x (Angstroms).
 * @param cell_b Unit cell dimension along y (Angstroms).
 * @param cell_c Unit cell dimension along z (Angstroms).
 * @return The interval count each edge selects: n_i, or n_i - 1 for a grid the
 *         caller is declaring to carry a duplicated closing plane.
 * @throws CellError If an edge does not round to a whole number of that axis's
 *         node intervals, if that number is neither n_i nor n_i - 1, or if the
 *         edge or the axis's spacing is not finite and positive.
 */
CellPeriods require_commensurate_cell(const GridParams& gp,
                                      double cell_a, double cell_b, double cell_c);

/**
 * @brief Periodic trilinear interpolation at a point, from derived geometry.
 *
 * The counterpart to interpolate_density_at for a map that tiles space. A point
 * in an axis's final interval blends across the period boundary instead of
 * falling off the end. Which node it blends into node 0 depends on how many
 * intervals the period holds, and that is `p_i = round(cell_i / spacing_i)`, the
 * count require_commensurate_cell derives from the caller's own cell. Both cases
 * occur here and both are definite. Under `p_i == n_i` the grid holds n_i
 * distinct nodes and node `n_i - 1` pairs with node 0. Under `p_i == n_i - 1` --
 * the caller declaring node `n_i - 1` to be a copy of node 0, which is what the
 * whole-cell assets in `tests/assets/mapq` carry -- node `n_i - 1` is never a
 * lower neighbour and it is node `n_i - 2` that pairs with node 0.
 * There is no outside to fall into, so @p default_value is returned
 * only when the point's fractional index is not finite -- which covers a
 * non-finite coordinate and also a finite one large enough that dividing it by
 * the spacing overflows.
 *
 * The wrap runs on the fractional index modulo that integer interval count
 * rather than on the Cartesian coordinate modulo the cell edge. The reduction is
 * then exact -- fmod is exact and the period is an integer -- so the wrapped
 * index carries only the error already in the fractional index, no matter how
 * many cells out the point started. Reducing the coordinate instead would
 * accumulate the cell edge's own rounding once per cell crossed.
 *
 * @param gp Geometry from get_grid_params.
 * @param values The grid's value array, from OESkewGrid::GetValues().
 * @param x Cartesian x coordinate.
 * @param y Cartesian y coordinate.
 * @param z Cartesian z coordinate.
 * @param cell_a Unit cell dimension along x (Angstroms).
 * @param cell_b Unit cell dimension along y (Angstroms).
 * @param cell_c Unit cell dimension along z (Angstroms).
 * @param default_value Value returned when the point's fractional index is not
 *        finite.
 * @return Interpolated density value.
 * @throws CellError If a cell edge does not round to n_i or n_i - 1 of that
 *         axis's node intervals, to within the allowance
 *         require_commensurate_cell makes for float node coordinates. Wrapping
 *         at a period the grid was never sampled on would resample the map onto
 *         a lattice it never had, so the caller is told rather than handed a
 *         plausible wrong number.
 */
double interpolate_density_periodic_at(const GridParams& gp, const float* values,
                                       double x, double y, double z,
                                       double cell_a, double cell_b, double cell_c,
                                       double default_value = 0.0);

/**
 * @brief True when two grids describe the same sampling of the same region.
 *
 * Equal dims, and per-axis spacing, node origin, and cell parameters within a
 * relative tolerance of @p tol (floored at 1.0 to keep sub-Angstrom quantities
 * absolute). Equal unit-cell presence.
 *
 * Both operands go through get_grid_params, so a grid whose geometry cannot be
 * derived throws here where OEGridSameGeometry returned false.
 *
 * @param tol Relative tolerance with unit floor (default 1e-6).
 * @throws GridError, CellError As get_grid_params, for either operand.
 */
bool same_grid_geometry(const OESystem::OESkewGrid& lhs,
                        const OESystem::OESkewGrid& rhs,
                        double tol = 1e-6);

/**
 * @brief Copy grid values to a flat vector (x-fastest order).
 *
 * Element index is iz * x_dim * y_dim + iy * x_dim + ix.
 *
 * @param grid Input grid.
 * @return Vector of grid values.
 */
std::vector<double> grid_to_vector(const OESystem::OESkewGrid& grid);

/**
 * @brief Copy values from a flat vector back into a grid.
 *
 * @param values Input values; the size must equal the grid's element count.
 * @param grid Output grid (modified in place).
 * @throws GridError If values.size() does not equal grid.GetSize(). The
 *         previous behavior silently copied the shorter of the two and left
 *         the rest of the grid holding whatever it held before.
 */
void vector_to_grid(const std::vector<double>& values, OESystem::OESkewGrid& grid);

/**
 * @brief Trilinear interpolation at a Cartesian point.
 *
 * The domain is the node span, `[origin_i, origin_i + (n_i - 1) * spacing_i]` on
 * each axis, plus a boundary tolerance. That is narrower by half a node interval
 * per face than the bounding box the scalar carrier's IsInGrid admitted.
 *
 * @param grid Input grid.
 * @param x Cartesian x coordinate.
 * @param y Cartesian y coordinate.
 * @param z Cartesian z coordinate.
 * @param default_value Value returned for a point outside the node span, or for
 *        a non-finite coordinate.
 * @return Interpolated density value.
 * @throws GridError, CellError As get_grid_params.
 */
double interpolate_density(const OESystem::OESkewGrid& grid,
                           double x, double y, double z,
                           double default_value = 0.0);

/**
 * @brief Batch trilinear interpolation at multiple points.
 *
 * The domain is the node span, as interpolate_density.
 *
 * @param grid Input grid.
 * @param points Flat array of {x0,y0,z0, x1,y1,z1, ...} coordinates.
 * @param num_points Number of points (points.size() / 3).
 * @param default_value Value returned for a point outside the node span, or for
 *        a non-finite coordinate.
 * @return Vector of interpolated values.
 * @throws GridError, CellError As get_grid_params.
 */
std::vector<double> interpolate_density_batch(
    const OESystem::OESkewGrid& grid,
    const std::vector<double>& points,
    size_t num_points,
    double default_value = 0.0);

/**
 * @brief Periodic-aware trilinear interpolation at a Cartesian point.
 *
 * The grid holds one period of a map that tiles space, so node `n_i - 1`'s upper
 * neighbour is node 0 and no point is outside the map. See
 * interpolate_density_periodic_at for the wrap, the one case @p default_value
 * still covers, and the commensurability requirement.
 *
 * @param grid Input grid.
 * @param x Cartesian x coordinate.
 * @param y Cartesian y coordinate.
 * @param z Cartesian z coordinate.
 * @param cell_a Unit cell dimension along x (Angstroms).
 * @param cell_b Unit cell dimension along y (Angstroms).
 * @param cell_c Unit cell dimension along z (Angstroms).
 * @param default_value Value returned when a point's fractional index is not
 *        finite.
 * @return Interpolated density value.
 * @throws GridError As get_grid_params.
 * @throws CellError As get_grid_params, or if a cell edge is not the extent the
 *         grid samples on that axis.
 */
double interpolate_density_periodic(
    const OESystem::OESkewGrid& grid,
    double x, double y, double z,
    double cell_a, double cell_b, double cell_c,
    double default_value = 0.0);

/**
 * @brief Batch periodic-aware trilinear interpolation at multiple points.
 *
 * As interpolate_density_periodic. The cell is checked once for the whole batch.
 *
 * @param grid Input grid.
 * @param points Flat array of {x0,y0,z0, x1,y1,z1, ...} coordinates.
 * @param num_points Number of points (points.size() / 3).
 * @param cell_a Unit cell dimension along x (Angstroms).
 * @param cell_b Unit cell dimension along y (Angstroms).
 * @param cell_c Unit cell dimension along z (Angstroms).
 * @param default_value Value returned when a point's fractional index is not
 *        finite.
 * @return Vector of interpolated values.
 * @throws GridError As get_grid_params.
 * @throws CellError As get_grid_params, or if a cell edge is not the extent the
 *         grid samples on that axis.
 */
std::vector<double> interpolate_density_periodic_batch(
    const OESystem::OESkewGrid& grid,
    const std::vector<double>& points,
    size_t num_points,
    double cell_a, double cell_b, double cell_c,
    double default_value = 0.0);

/**
 * @brief Collect grid element indices within a sphere.
 *
 * Returns the flat indices of all grid elements whose centers lie
 * within the specified radius of the query point.
 *
 * @param grid Input grid.
 * @param x Center x coordinate.
 * @param y Center y coordinate.
 * @param z Center z coordinate.
 * @param radius Search radius in Angstroms.
 * @return Vector of grid element indices.
 */
std::vector<unsigned int> get_atom_grid_points(
    const OESystem::OESkewGrid& grid,
    double x, double y, double z, double radius);

}  // namespace Maptitude

#endif  // MAPTITUDE_GRID_H
