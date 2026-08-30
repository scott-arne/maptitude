#include "maptitude/Grid.h"
#include "maptitude/Error.h"

#include <oegrid.h>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>

namespace Maptitude {

namespace {

/// Read one node's Cartesian coordinate, enforcing derivation checks 2 and 3.
/// @p axis names the axis whose walk needs this node, for the error text.
void ReadNodeCoord(const OESystem::OESkewGrid& grid, const unsigned int element,
                   const char* axis, double out[3]) {
    static const char* const COMPONENT[3] = {"x", "y", "z"};
    float f[3] = {0.0f, 0.0f, 0.0f};

    if (!grid.ElementToSpatialCoord(element, f[0], f[1], f[2])) {
        std::ostringstream message;
        message << "Grid element " << element << " (walking axis " << axis
                << ") has no spatial coordinate; the geometry cannot be derived";
        throw GridError(message.str());
    }

    for (int j = 0; j < 3; ++j) {
        if (!std::isfinite(f[j])) {
            std::ostringstream message;
            message << "Grid element " << element << " (walking axis " << axis
                    << ") has a non-finite " << COMPONENT[j] << " coordinate ("
                    << f[j] << "); the geometry cannot be derived";
            throw GridError(message.str());
        }
        out[j] = f[j];
    }
}

/// The containment predicate, taking an already-computed fractional index.
///
/// The spec requires `interpolate_density` to share `grid_contains`'s predicate
/// rather than restate it, so the prechecks and the interpolator cannot drift
/// apart. A direct call would not serve: `interpolate_density_at` (Task 2, same
/// translation unit) needs the fractional index itself for the blend, and
/// calling `grid_contains(gp, x, y, z)` would recompute `grid_fractional_index`
/// a second time in the hottest loop in the library. Taking the index instead
/// gives both callers one predicate and one index computation each.
///
/// Written positively: a NaN fractional index must report "outside", and a
/// negated range comparison would pass it.
bool ContainsFractionalIndex(const GridParams& gp,
                             const double fx, const double fy, const double fz) {
    return std::isfinite(fx) && std::isfinite(fy) && std::isfinite(fz) &&
           fx >= 0.0 && fx <= gp.x_dim - 1.0 &&
           fy >= 0.0 && fy <= gp.y_dim - 1.0 &&
           fz >= 0.0 && fz <= gp.z_dim - 1.0;
}

/// True when two quantities agree to within @p tol scaled by their magnitude.
///
/// The operands descend from OpenEye's float grid coordinates, where one ulp
/// exceeds an absolute 1e-6 above roughly 8.4 Angstroms — so an absolute
/// comparison is exact float equality on every real crystallographic cell.
/// The unit floor keeps the tolerance absolute for sub-Angstrom quantities
/// such as node spacings, and gives the degree-valued cell angles a scale
/// their own magnitude supplies.
bool NearlyEqual(const double p, const double q, const double tol) {
    return std::abs(p - q) <= tol * std::max({1.0, std::abs(p), std::abs(q)});
}

}  // namespace

GridParams get_grid_params(const OESystem::OESkewGrid& grid) {
    static const char* const AXIS[3] = {"x", "y", "z"};
    // Off-axis drift above this over an axis's full span means the sampling is
    // not axis-aligned, so no per-axis spacing describes it.
    constexpr double MAX_OFF_AXIS_LEAK = 1e-4;

    const unsigned int n[3] = {grid.GetXDim(), grid.GetYDim(), grid.GetZDim()};

    // Check 1 runs over all three axes before any walk, so a grid that fails
    // both check 1 and check 5 reports check 1.
    for (int i = 0; i < 3; ++i) {
        if (n[i] < 2) {
            std::ostringstream message;
            message << "Grid axis " << AXIS[i] << " has dimension " << n[i]
                    << "; deriving a node interval needs at least 2 nodes on every axis";
            throw GridError(message.str());
        }
    }

    // Elements linearize x-fastest: el = iz*nx*ny + iy*nx + ix.
    const unsigned int step[3] = {1u, n[0], n[0] * n[1]};

    double origin[3];
    ReadNodeCoord(grid, 0u, AXIS[0], origin);

    double spacing[3];
    for (int i = 0; i < 3; ++i) {
        double node[3];
        ReadNodeCoord(grid, (n[i] - 1u) * step[i], AXIS[i], node);

        for (int j = 0; j < 3; ++j) {
            if (j == i) continue;
            const double leak = std::abs(node[j] - origin[j]);
            if (leak > MAX_OFF_AXIS_LEAK) {
                std::ostringstream message;
                message << "Grid axis " << AXIS[i] << " leaks " << leak
                        << " A into " << AXIS[j] << " over its full span (limit "
                        << MAX_OFF_AXIS_LEAK << " A); maptitude requires axis-aligned sampling";
                throw CellError(message.str());
            }
        }

        // Averaging over every interval on the axis, rather than measuring one
        // adjacent step, divides ElementToSpatialCoord's float noise by n_i - 1.
        // On 1d26 that is the difference between 0.902082443 (single step) and
        // 0.902083317 (full walk) against a true 0.902083333.
        spacing[i] = (node[i] - origin[i]) / (n[i] - 1u);

        if (!std::isfinite(spacing[i]) || spacing[i] <= 0.0) {
            std::ostringstream message;
            message << "Grid axis " << AXIS[i] << " has a derived node interval of "
                    << spacing[i] << " A; the interval must be finite and positive";
            throw GridError(message.str());
        }
    }

    return GridParams{origin[0], origin[1], origin[2],
                      n[0],      n[1],      n[2],
                      spacing[0], spacing[1], spacing[2]};
}

UnitCellParams get_unit_cell(const OESystem::OESkewGrid& grid) {
    float a = 0.0f, b = 0.0f, c = 0.0f, alpha = 0.0f, beta = 0.0f, gamma = 0.0f;
    if (!grid.HasUnitCell() || !grid.GetUnitCell(a, b, c, alpha, beta, gamma)) {
        throw CellError("Grid has no unit cell; its cell parameters cannot be read");
    }
    return UnitCellParams{a, b, c, alpha, beta, gamma};
}

void grid_node_origin(const GridParams& gp, double& x, double& y, double& z) {
    x = gp.x_origin;
    y = gp.y_origin;
    z = gp.z_origin;
}

void grid_bounds(const GridParams& gp,
                 double& xmin, double& ymin, double& zmin,
                 double& xmax, double& ymax, double& zmax) {
    // The scalar carrier's box ran half a spacing outside the first and last
    // nodes on each face: GetXMin() == GetXMid() - n_x*s_x/2 while element 0
    // sits at GetXMid() - (n_x-1)*s_x/2.
    xmin = gp.x_origin - gp.x_spacing / 2.0;
    ymin = gp.y_origin - gp.y_spacing / 2.0;
    zmin = gp.z_origin - gp.z_spacing / 2.0;
    xmax = gp.x_origin + (gp.x_dim - 0.5) * gp.x_spacing;
    ymax = gp.y_origin + (gp.y_dim - 0.5) * gp.y_spacing;
    zmax = gp.z_origin + (gp.z_dim - 0.5) * gp.z_spacing;
}

void grid_fractional_index(const GridParams& gp,
                           const double x, const double y, const double z,
                           double& fx, double& fy, double& fz) {
    fx = (x - gp.x_origin) / gp.x_spacing;
    fy = (y - gp.y_origin) / gp.y_spacing;
    fz = (z - gp.z_origin) / gp.z_spacing;
}

bool grid_contains(const GridParams& gp,
                   const double x, const double y, const double z) {
    double fx = 0.0, fy = 0.0, fz = 0.0;
    grid_fractional_index(gp, x, y, z, fx, fy, fz);
    return ContainsFractionalIndex(gp, fx, fy, fz);
}

double interpolate_density_at(const GridParams& gp, const float* values,
                              const double x, const double y, const double z,
                              const double default_value) {
    double f[3] = {0.0, 0.0, 0.0};
    grid_fractional_index(gp, x, y, z, f[0], f[1], f[2]);

    // Steps 2 and 3 of the algorithm, through the one predicate grid_contains
    // also calls, so the containment prechecks and the interpolator cannot drift
    // apart. Passing the index rather than the coordinate keeps this to one
    // grid_fractional_index per point in the batch loops.
    if (!ContainsFractionalIndex(gp, f[0], f[1], f[2])) {
        return default_value;
    }

    // Clamping the base index keeps a point exactly on the far face inside the
    // last cell with t == 1.0, rather than indexing one node past the end. The
    // lower clamp is belt-and-braces: the containment check already rejects f < 0.
    const unsigned int n[3] = {gp.x_dim, gp.y_dim, gp.z_dim};
    unsigned int i0[3];
    double t[3];
    for (int i = 0; i < 3; ++i) {
        const double base = std::floor(f[i]);
        const double clamped = std::min(std::max(base, 0.0),
                                        static_cast<double>(n[i] - 2u));
        i0[i] = static_cast<unsigned int>(clamped);
        t[i] = f[i] - clamped;
    }

    const unsigned int stride_y = gp.x_dim;
    const unsigned int stride_z = gp.x_dim * gp.y_dim;
    const unsigned int base = i0[2] * stride_z + i0[1] * stride_y + i0[0];

    const double c000 = values[base];
    const double c100 = values[base + 1u];
    const double c010 = values[base + stride_y];
    const double c110 = values[base + stride_y + 1u];
    const double c001 = values[base + stride_z];
    const double c101 = values[base + stride_z + 1u];
    const double c011 = values[base + stride_z + stride_y];
    const double c111 = values[base + stride_z + stride_y + 1u];

    const double c00 = c000 * (1.0 - t[0]) + c100 * t[0];
    const double c10 = c010 * (1.0 - t[0]) + c110 * t[0];
    const double c01 = c001 * (1.0 - t[0]) + c101 * t[0];
    const double c11 = c011 * (1.0 - t[0]) + c111 * t[0];

    const double c0 = c00 * (1.0 - t[1]) + c10 * t[1];
    const double c1 = c01 * (1.0 - t[1]) + c11 * t[1];

    return c0 * (1.0 - t[2]) + c1 * t[2];
}

bool same_grid_geometry(const OESystem::OESkewGrid& lhs,
                        const OESystem::OESkewGrid& rhs,
                        const double tol) {
    const GridParams a = get_grid_params(lhs);
    const GridParams b = get_grid_params(rhs);

    if (a.x_dim != b.x_dim || a.y_dim != b.y_dim || a.z_dim != b.z_dim) return false;
    if (!NearlyEqual(a.x_spacing, b.x_spacing, tol)) return false;
    if (!NearlyEqual(a.y_spacing, b.y_spacing, tol)) return false;
    if (!NearlyEqual(a.z_spacing, b.z_spacing, tol)) return false;
    if (!NearlyEqual(a.x_origin, b.x_origin, tol)) return false;
    if (!NearlyEqual(a.y_origin, b.y_origin, tol)) return false;
    if (!NearlyEqual(a.z_origin, b.z_origin, tol)) return false;

    if (lhs.HasUnitCell() != rhs.HasUnitCell()) return false;
    if (lhs.HasUnitCell()) {
        const UnitCellParams ca = get_unit_cell(lhs);
        const UnitCellParams cb = get_unit_cell(rhs);
        if (!NearlyEqual(ca.a, cb.a, tol)) return false;
        if (!NearlyEqual(ca.b, cb.b, tol)) return false;
        if (!NearlyEqual(ca.c, cb.c, tol)) return false;
        if (!NearlyEqual(ca.alpha, cb.alpha, tol)) return false;
        if (!NearlyEqual(ca.beta, cb.beta, tol)) return false;
        if (!NearlyEqual(ca.gamma, cb.gamma, tol)) return false;
    }
    return true;
}

std::vector<double> grid_to_vector(const OESystem::OESkewGrid& grid) {
    const unsigned int size = grid.GetSize();
    const float* values = grid.GetValues();
    std::vector<double> result(size);
    for (unsigned int i = 0; i < size; ++i) {
        result[i] = values[i];
    }
    return result;
}

void vector_to_grid(const std::vector<double>& values, OESystem::OESkewGrid& grid) {
    const unsigned int size = grid.GetSize();
    if (values.size() != size) {
        std::ostringstream message;
        message << "vector_to_grid needs exactly " << size << " values for this grid, got "
                << values.size() << "; a short vector previously left the tail of the grid "
                   "holding stale density";
        throw GridError(message.str());
    }
    float* out = grid.GetValues();
    for (unsigned int i = 0; i < size; ++i) {
        out[i] = static_cast<float>(values[i]);
    }
}

double interpolate_density(const OESystem::OESkewGrid& grid,
                           const double x, const double y, const double z,
                           const double default_value) {
    const GridParams gp = get_grid_params(grid);
    return interpolate_density_at(gp, grid.GetValues(), x, y, z, default_value);
}

std::vector<double> interpolate_density_batch(
    const OESystem::OESkewGrid& grid,
    const std::vector<double>& points,
    const size_t num_points,
    const double default_value) {
    const GridParams gp = get_grid_params(grid);
    const float* values = grid.GetValues();
    std::vector<double> result(num_points);
    for (size_t i = 0; i < num_points; ++i) {
        result[i] = interpolate_density_at(
            gp, values, points[i * 3], points[i * 3 + 1], points[i * 3 + 2],
            default_value);
    }
    return result;
}

namespace {

/// Wrap into the cell relative to the node origin, then interpolate.
/// The wrap origin was GetXMin(), half a spacing below the first node; it is
/// the node origin now, so wrapped points land on the sampled lattice.
double InterpolateDensityPeriodicAt(
    const GridParams& gp, const float* values,
    const double x, const double y, const double z,
    const double cell_a, const double cell_b, const double cell_c,
    const double default_value) {
    double wx = gp.x_origin + std::fmod(x - gp.x_origin, cell_a);
    if (wx < gp.x_origin) wx += cell_a;
    double wy = gp.y_origin + std::fmod(y - gp.y_origin, cell_b);
    if (wy < gp.y_origin) wy += cell_b;
    double wz = gp.z_origin + std::fmod(z - gp.z_origin, cell_c);
    if (wz < gp.z_origin) wz += cell_c;
    return interpolate_density_at(gp, values, wx, wy, wz, default_value);
}

}  // namespace

double interpolate_density_periodic(
    const OESystem::OESkewGrid& grid,
    const double x, const double y, const double z,
    const double cell_a, const double cell_b, const double cell_c,
    const double default_value) {
    const GridParams gp = get_grid_params(grid);
    return InterpolateDensityPeriodicAt(gp, grid.GetValues(), x, y, z,
                                        cell_a, cell_b, cell_c, default_value);
}

std::vector<double> interpolate_density_periodic_batch(
    const OESystem::OESkewGrid& grid,
    const std::vector<double>& points,
    const size_t num_points,
    const double cell_a, const double cell_b, const double cell_c,
    const double default_value) {
    const GridParams gp = get_grid_params(grid);
    const float* values = grid.GetValues();
    std::vector<double> result(num_points);
    for (size_t i = 0; i < num_points; ++i) {
        result[i] = InterpolateDensityPeriodicAt(
            gp, values, points[i * 3], points[i * 3 + 1], points[i * 3 + 2],
            cell_a, cell_b, cell_c, default_value);
    }
    return result;
}

std::vector<unsigned int> get_atom_grid_points(
    const OESystem::OESkewGrid& grid,
    const double x, const double y, const double z, const double radius) {
    std::vector<unsigned int> result;

    const GridParams gp = get_grid_params(grid);
    const double r2 = radius * radius;
    const double origin[3] = {gp.x_origin, gp.y_origin, gp.z_origin};
    const double spacing[3] = {gp.x_spacing, gp.y_spacing, gp.z_spacing};
    const unsigned int dim[3] = {gp.x_dim, gp.y_dim, gp.z_dim};
    const double centre[3] = {x, y, z};

    int lo[3], hi[3];
    for (int i = 0; i < 3; ++i) {
        lo[i] = std::max(0, static_cast<int>(
            std::floor((centre[i] - radius - origin[i]) / spacing[i])));
        hi[i] = std::min(static_cast<int>(dim[i]) - 1, static_cast<int>(
            std::ceil((centre[i] + radius - origin[i]) / spacing[i])));
    }

    for (int ix = lo[0]; ix <= hi[0]; ++ix) {
        const double dx = origin[0] + ix * spacing[0] - x;
        for (int iy = lo[1]; iy <= hi[1]; ++iy) {
            const double dy = origin[1] + iy * spacing[1] - y;
            for (int iz = lo[2]; iz <= hi[2]; ++iz) {
                const double dz = origin[2] + iz * spacing[2] - z;
                if (dx * dx + dy * dy + dz * dz <= r2) {
                    result.push_back(static_cast<unsigned int>(
                        iz * static_cast<int>(dim[0]) * static_cast<int>(dim[1]) +
                        iy * static_cast<int>(dim[0]) + ix));
                }
            }
        }
    }

    return result;
}

}  // namespace Maptitude
