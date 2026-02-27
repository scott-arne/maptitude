#include "maptitude/Grid.h"

#include <oegrid.h>

#include <cmath>

namespace Maptitude {

GridParams get_grid_params(const OESystem::OEScalarGrid& grid) {
    float fx, fy, fz;
    grid.ElementToSpatialCoord(0, fx, fy, fz);
    return GridParams{
        static_cast<double>(fx),
        static_cast<double>(fy),
        static_cast<double>(fz),
        grid.GetXDim(),
        grid.GetYDim(),
        grid.GetZDim(),
        grid.GetSpacing()
    };
}

std::vector<double> grid_to_vector(const OESystem::OEScalarGrid& grid) {
    const size_t size = grid.GetSize();
    std::vector<double> result(size);
    for (size_t i = 0; i < size; ++i) {
        result[i] = grid[static_cast<unsigned int>(i)];
    }
    return result;
}

void vector_to_grid(const std::vector<double>& values, OESystem::OEScalarGrid& grid) {
    const size_t size = grid.GetSize();
    for (size_t i = 0; i < size && i < values.size(); ++i) {
        grid[static_cast<unsigned int>(i)] = static_cast<float>(values[i]);
    }
}

double interpolate_density(const OESystem::OEScalarGrid& grid,
                          const double x, const double y, const double z,
                          const double default_value) {
    if (!grid.IsInGrid(static_cast<float>(x),
                       static_cast<float>(y),
                       static_cast<float>(z))) {
        return default_value;
    }

    // Use OpenEye's built-in trilinear interpolation
    return OESystem::OEFloatGridLinearInterpolate(
        grid,
        static_cast<float>(x),
        static_cast<float>(y),
        static_cast<float>(z),
        static_cast<float>(default_value));
}

std::vector<double> interpolate_density_batch(
    const OESystem::OEScalarGrid& grid,
    const std::vector<double>& points,
    const size_t num_points,
    const double default_value) {
    std::vector<double> result(num_points);
    for (size_t i = 0; i < num_points; ++i) {
        result[i] = interpolate_density(
            grid, points[i * 3], points[i * 3 + 1], points[i * 3 + 2],
            default_value);
    }
    return result;
}

double interpolate_density_periodic(
    const OESystem::OEScalarGrid& grid,
    const double x, const double y, const double z,
    const double cell_a, const double cell_b, const double cell_c,
    const double default_value) {
    // Wrap coordinates modulo unit cell dimensions relative to grid origin
    const double orig_x = grid.GetXMin();
    const double orig_y = grid.GetYMin();
    const double orig_z = grid.GetZMin();

    double wx = orig_x + std::fmod(x - orig_x, cell_a);
    if (wx < orig_x) wx += cell_a;

    double wy = orig_y + std::fmod(y - orig_y, cell_b);
    if (wy < orig_y) wy += cell_b;

    double wz = orig_z + std::fmod(z - orig_z, cell_c);
    if (wz < orig_z) wz += cell_c;

    return interpolate_density(grid, wx, wy, wz, default_value);
}

std::vector<double> interpolate_density_periodic_batch(
    const OESystem::OEScalarGrid& grid,
    const std::vector<double>& points,
    const size_t num_points,
    const double cell_a, const double cell_b, const double cell_c,
    const double default_value) {
    std::vector<double> result(num_points);
    for (size_t i = 0; i < num_points; ++i) {
        result[i] = interpolate_density_periodic(
            grid, points[i * 3], points[i * 3 + 1], points[i * 3 + 2],
            cell_a, cell_b, cell_c, default_value);
    }
    return result;
}

std::vector<unsigned int> get_atom_grid_points(
    const OESystem::OEScalarGrid& grid,
    const double x, const double y, const double z, const double radius) {
    std::vector<unsigned int> result;

    const double spacing = grid.GetSpacing();
    const double r2 = radius * radius;

    // Grid node origin (ElementToSpatialCoord gives the actual first node
    // position, whereas GetXMin returns the bounding-box edge which is
    // offset by half a spacing)
    float fx0, fy0, fz0;
    grid.ElementToSpatialCoord(0, fx0, fy0, fz0);
    const double x0 = fx0;
    const double y0 = fy0;
    const double z0 = fz0;
    const unsigned int xdim = grid.GetXDim();
    const unsigned int ydim = grid.GetYDim();
    const unsigned int zdim = grid.GetZDim();

    // Compute index range to search
    const int ix_min = std::max(0, static_cast<int>(std::floor((x - radius - x0) / spacing)));
    const int ix_max = std::min(static_cast<int>(xdim) - 1,
                          static_cast<int>(std::ceil((x + radius - x0) / spacing)));
    const int iy_min = std::max(0, static_cast<int>(std::floor((y - radius - y0) / spacing)));
    const int iy_max = std::min(static_cast<int>(ydim) - 1,
                          static_cast<int>(std::ceil((y + radius - y0) / spacing)));
    const int iz_min = std::max(0, static_cast<int>(std::floor((z - radius - z0) / spacing)));
    const int iz_max = std::min(static_cast<int>(zdim) - 1,
                          static_cast<int>(std::ceil((z + radius - z0) / spacing)));

    for (int ix = ix_min; ix <= ix_max; ++ix) {
        const double gx = x0 + ix * spacing;
        const double dx = gx - x;
        for (int iy = iy_min; iy <= iy_max; ++iy) {
            const double gy = y0 + iy * spacing;
            const double dy = gy - y;
            for (int iz = iz_min; iz <= iz_max; ++iz) {
                const double gz = z0 + iz * spacing;
                const double dz = gz - z;
                if (dx * dx + dy * dy + dz * dz <= r2) {
                    const unsigned int idx = static_cast<unsigned int>(
                        iz * xdim * ydim + iy * xdim + ix);
                    result.push_back(idx);
                }
            }
        }
    }

    return result;
}

}  // namespace Maptitude
