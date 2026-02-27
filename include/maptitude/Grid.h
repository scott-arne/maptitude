/**
 * @file Grid.h
 * @brief Utilities for working with OEScalarGrid objects.
 *
 * Provides helper functions for converting between OpenEye grid
 * representations and raw data arrays used by the computation kernels.
 */

#ifndef MAPTITUDE_GRID_H
#define MAPTITUDE_GRID_H

#include <cstddef>
#include <vector>

namespace OESystem {
class OEScalarGrid;
}

namespace Maptitude {

/**
 * @brief Parameters describing a grid's geometry.
 */
struct GridParams {
    double x_origin;     ///< X origin (Angstroms)
    double y_origin;     ///< Y origin (Angstroms)
    double z_origin;     ///< Z origin (Angstroms)
    unsigned int x_dim;  ///< Number of grid points in X
    unsigned int y_dim;  ///< Number of grid points in Y
    unsigned int z_dim;  ///< Number of grid points in Z
    double spacing;      ///< Grid spacing (Angstroms)
};

/**
 * @brief Extract geometry parameters from an OEScalarGrid.
 *
 * @param grid Input grid.
 * @return GridParams describing the grid geometry.
 */
GridParams get_grid_params(const OESystem::OEScalarGrid& grid);

/**
 * @brief Copy grid values to a flat vector (z-fastest order).
 *
 * @param grid Input grid.
 * @return Vector of grid values.
 */
std::vector<double> grid_to_vector(const OESystem::OEScalarGrid& grid);

/**
 * @brief Copy values from a flat vector back into a grid.
 *
 * @param values Input values (must match grid size).
 * @param grid Output grid (modified in place).
 */
void vector_to_grid(const std::vector<double>& values, OESystem::OEScalarGrid& grid);

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
double interpolate_density(const OESystem::OEScalarGrid& grid,
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
    const OESystem::OEScalarGrid& grid,
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
    const OESystem::OEScalarGrid& grid,
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
    const OESystem::OEScalarGrid& grid,
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
    const OESystem::OEScalarGrid& grid,
    double x, double y, double z, double radius);

}  // namespace Maptitude

#endif  // MAPTITUDE_GRID_H
