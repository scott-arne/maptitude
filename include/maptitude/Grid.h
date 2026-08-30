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
 * @brief Bounding-box corners, half a spacing outside the first and last nodes.
 *
 * This reproduces what the scalar carrier's GetXMin/GetXMax family returned. It is
 * for diagnostics and geometry descriptions only: it is deliberately NOT the
 * domain over which interpolate_density returns interpolated data. Use
 * grid_contains for that.
 */
void grid_bounds(const GridParams& gp,
                 double& xmin, double& ymin, double& zmin,
                 double& xmax, double& ymax, double& zmax);

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
 * The node span is closed on both ends. A point exactly on the far face
 * interpolates rather than falling out of the grid.
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
 * @param grid Input grid.
 * @param x Cartesian x coordinate.
 * @param y Cartesian y coordinate.
 * @param z Cartesian z coordinate.
 * @param default_value Value to return if point is outside grid.
 * @return Interpolated density value.
 */
double interpolate_density(const OESystem::OESkewGrid& grid,
                           double x, double y, double z,
                           double default_value = 0.0);

/**
 * @brief Batch trilinear interpolation at multiple points.
 *
 * @param grid Input grid.
 * @param points Flat array of {x0,y0,z0, x1,y1,z1, ...} coordinates.
 * @param num_points Number of points (points.size() / 3).
 * @param default_value Value for out-of-bounds points.
 * @return Vector of interpolated values.
 */
std::vector<double> interpolate_density_batch(
    const OESystem::OESkewGrid& grid,
    const std::vector<double>& points,
    size_t num_points,
    double default_value = 0.0);

/**
 * @brief Periodic-aware trilinear interpolation at a Cartesian point.
 *
 * Wraps coordinates modulo the unit cell dimensions before interpolating,
 * so that points outside the grid are mapped back into the unit cell.
 *
 * @param grid Input grid.
 * @param x Cartesian x coordinate.
 * @param y Cartesian y coordinate.
 * @param z Cartesian z coordinate.
 * @param cell_a Unit cell dimension along x (Angstroms).
 * @param cell_b Unit cell dimension along y (Angstroms).
 * @param cell_c Unit cell dimension along z (Angstroms).
 * @param default_value Value to return if point is still outside grid after wrapping.
 * @return Interpolated density value.
 */
double interpolate_density_periodic(
    const OESystem::OESkewGrid& grid,
    double x, double y, double z,
    double cell_a, double cell_b, double cell_c,
    double default_value = 0.0);

/**
 * @brief Batch periodic-aware trilinear interpolation at multiple points.
 *
 * @param grid Input grid.
 * @param points Flat array of {x0,y0,z0, x1,y1,z1, ...} coordinates.
 * @param num_points Number of points (points.size() / 3).
 * @param cell_a Unit cell dimension along x (Angstroms).
 * @param cell_b Unit cell dimension along y (Angstroms).
 * @param cell_c Unit cell dimension along z (Angstroms).
 * @param default_value Value for out-of-bounds points after wrapping.
 * @return Vector of interpolated values.
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
