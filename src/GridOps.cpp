#include "maptitude/GridOps.h"
#include "maptitude/Error.h"
#include "maptitude/Grid.h"

#include <oechem.h>
#include <oegrid.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

namespace Maptitude {

namespace {

/// Throw naming the setter and its arguments when an OESkewGrid setter fails,
/// deleting the partially built grid first. A false return leaves the grid
/// holding its previous geometry, which would otherwise be filled with density
/// sampled for a different box.
template <typename... Args>
void RequireSetter(OESystem::OESkewGrid* grid, const bool ok,
                   const char* setter, Args... args) {
    if (ok) return;
    std::ostringstream message;
    message << "OESkewGrid::" << setter << " rejected (";
    const double values[] = {static_cast<double>(args)...};
    for (size_t i = 0; i < sizeof...(args); ++i) {
        if (i) message << ", ";
        message << values[i];
    }
    message << ") while building the padded grid";
    delete grid;
    throw GridError(message.str());
}

}  // namespace

/// Render a grid's geometry for an error message: dimensions, centre, spacing.
static std::string DescribeGeometry(const OESystem::OESkewGrid& grid) {
    const GridParams gp = get_grid_params(grid);
    std::ostringstream out;
    out << std::setprecision(std::numeric_limits<float>::max_digits10);
    out << gp.x_dim << "x" << gp.y_dim << "x" << gp.z_dim
        << " centred at (" << grid.GetXMid() << ", " << grid.GetYMid() << ", "
        << grid.GetZMid() << ") spacing (" << gp.x_spacing << ", "
        << gp.y_spacing << ", " << gp.z_spacing << ")";
    return out.str();
}

void scale_map(OESystem::OESkewGrid& grid, const double factor) {
    const unsigned int size = grid.GetSize();
    float* values = grid.GetValues();
    for (unsigned int i = 0; i < size; ++i) {
        values[i] = static_cast<float>(values[i] * factor);
    }
}

OESystem::OESkewGrid* combine_maps(
    const OESystem::OESkewGrid& lhs,
    const OESystem::OESkewGrid& rhs,
    const MapOp op) {
    // same_grid_geometry compares dimensions, node origin, per-axis spacing and
    // the unit cell. The previous hand-rolled check ignored the origin, so grids
    // of the same shape at different positions were combined element-wise --
    // mixing densities from different points in space.
    if (!same_grid_geometry(lhs, rhs)) {
        throw GridError("Grids must have identical geometry for combination: left is " +
                        DescribeGeometry(lhs) + ", right is " + DescribeGeometry(rhs));
    }

    auto* result = new OESystem::OESkewGrid(lhs);
    const unsigned int size = result->GetSize();
    const float* lv = lhs.GetValues();
    const float* rv = rhs.GetValues();
    float* out = result->GetValues();

    for (unsigned int i = 0; i < size; ++i) {
        const float lval = lv[i];
        const float rval = rv[i];
        float combined = 0.0f;

        switch (op) {
            case MapOp::ADD:      combined = lval + rval;          break;
            case MapOp::SUBTRACT: combined = lval - rval;          break;
            case MapOp::MIN:      combined = std::min(lval, rval); break;
            case MapOp::MAX:      combined = std::max(lval, rval); break;
        }

        out[i] = combined;
    }

    return result;
}

OESystem::OESkewGrid* diff_to_calc(
    const OESystem::OESkewGrid& obs_grid,
    const OESystem::OESkewGrid& diff_grid) {
    if (!same_grid_geometry(obs_grid, diff_grid)) {
        throw GridError("Observed and difference grids must have identical geometry: observed is " +
                        DescribeGeometry(obs_grid) + ", difference is " + DescribeGeometry(diff_grid));
    }

    auto* result = new OESystem::OESkewGrid(obs_grid);
    const unsigned int size = result->GetSize();
    const float* obs = obs_grid.GetValues();
    const float* diff = diff_grid.GetValues();
    float* out = result->GetValues();

    for (unsigned int i = 0; i < size; ++i) {
        // rho_calc = rho_obs - 2 * rho_diff
        out[i] = obs[i] - 2.0f * diff[i];
    }

    return result;
}

OESystem::OESkewGrid* wrap_and_pad_grid(
    const OESystem::OESkewGrid& grid,
    OEChem::OEMolBase& mol,
    const double cell_a, const double cell_b, const double cell_c,
    const double padding) {
    // The periodic wrap below is `std::fmod(offset, cell_edge)`, which is NaN for a zero
    // divisor and meaningless for a non-finite one, so every voxel of the padded grid
    // comes back NaN with no error. Validate all three edges up front rather than
    // guarding only the centroid shift.
    const double edges[3] = {cell_a, cell_b, cell_c};
    const char* edge_names[3] = {"cell_a", "cell_b", "cell_c"};
    for (int i = 0; i < 3; ++i) {
        if (!std::isfinite(edges[i]) || edges[i] <= 0.0) {
            std::ostringstream message;
            message << "wrap_and_pad_grid requires a finite positive " << edge_names[i]
                    << " (got " << edges[i] << ")";
            throw CellError(message.str());
        }
    }

    const GridParams gp = get_grid_params(grid);

    // Compute heavy-atom centroid
    double cx = 0.0, cy = 0.0, cz = 0.0;
    int n = 0;
    float coords[3];
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(OEChem::OEIsHeavy()); atom; ++atom) {
        mol.GetCoords(&(*atom), coords);
        cx += coords[0];
        cy += coords[1];
        cz += coords[2];
        ++n;
    }
    if (n == 0) {
        // nullptr from this function means "no padding was needed". An empty heavy
        // atom set is a different condition entirely and must not share that
        // signal -- the Python wrapper maps nullptr to "return the grid unchanged".
        throw StructureError("wrap_and_pad_grid requires at least one heavy atom; the molecule has " +
                             std::to_string(mol.NumAtoms()) + " atom(s), none of them heavy");
    }
    cx /= n;
    cy /= n;
    cz /= n;

    // Shift centroid to grid centre using integer unit-cell vectors
    const double grid_xmid = grid.GetXMid();
    const double grid_ymid = grid.GetYMid();
    const double grid_zmid = grid.GetZMid();

    const double shift_x = std::round((grid_xmid - cx) / cell_a) * cell_a;
    const double shift_y = std::round((grid_ymid - cy) / cell_b) * cell_b;
    const double shift_z = std::round((grid_zmid - cz) / cell_c) * cell_c;

    if (std::abs(shift_x) > 0.01 || std::abs(shift_y) > 0.01 || std::abs(shift_z) > 0.01) {
        for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
            mol.GetCoords(&(*atom), coords);
            float shifted[3] = {
                static_cast<float>(coords[0] + shift_x),
                static_cast<float>(coords[1] + shift_y),
                static_cast<float>(coords[2] + shift_z)
            };
            mol.SetCoords(&(*atom), shifted);
        }
    }

    // Check whether all heavy atoms fall inside the grid (with padding)
    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    double min_z = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double max_y = std::numeric_limits<double>::lowest();
    double max_z = std::numeric_limits<double>::lowest();

    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(OEChem::OEIsHeavy()); atom; ++atom) {
        mol.GetCoords(&(*atom), coords);
        min_x = std::min(min_x, static_cast<double>(coords[0]));
        min_y = std::min(min_y, static_cast<double>(coords[1]));
        min_z = std::min(min_z, static_cast<double>(coords[2]));
        max_x = std::max(max_x, static_cast<double>(coords[0]));
        max_y = std::max(max_y, static_cast<double>(coords[1]));
        max_z = std::max(max_z, static_cast<double>(coords[2]));
    }

    // The interpolatable domain is the node span, not the bounding box, so the
    // padding test asks whether the atoms fit inside the nodes. Relative to the
    // old GetXMin/GetXMax box each edge moves inward by half a spacing, so the
    // tested region shrinks by one spacing per axis and needs_pad can only flip
    // false to true, never the reverse.
    const double grid_xmin = gp.x_origin;
    const double grid_ymin = gp.y_origin;
    const double grid_zmin = gp.z_origin;
    const double grid_xmax = gp.x_origin + (gp.x_dim - 1) * gp.x_spacing;
    const double grid_ymax = gp.y_origin + (gp.y_dim - 1) * gp.y_spacing;
    const double grid_zmax = gp.z_origin + (gp.z_dim - 1) * gp.z_spacing;

    const bool needs_pad =
        (min_x - padding < grid_xmin) ||
        (max_x + padding > grid_xmax) ||
        (min_y - padding < grid_ymin) ||
        (max_y + padding > grid_ymax) ||
        (min_z - padding < grid_zmin) ||
        (max_z + padding > grid_zmax);

    if (!needs_pad) return nullptr;

    // The scalar carrier's extents-box constructor built the padded grid;
    // OESkewGrid has no equivalent, so reproduce it explicitly. Probed against
    // 2026.1.0: for minmax {0,0,0, 9.5,9.5,9.5} at spacing 1.0 this gives dims
    // 10^3, mid 4.75 and node 0 at 0.25 -- the node origin is NOT minmax[0].
    const double minmax[6] = {
        min_x - padding, min_y - padding, min_z - padding,
        max_x + padding, max_y + padding, max_z + padding
    };
    const double src_spacing[3] = {gp.x_spacing, gp.y_spacing, gp.z_spacing};
    unsigned int pad_dim[3];
    double pad_mid[3];
    for (int i = 0; i < 3; ++i) {
        const double extent = minmax[i + 3] - minmax[i];
        pad_dim[i] = static_cast<unsigned int>(extent / src_spacing[i]) + 1u;
        pad_mid[i] = (minmax[i] + minmax[i + 3]) / 2.0;
    }

    auto* padded = new OESystem::OESkewGrid();
    RequireSetter(padded, padded->SetDim(pad_dim[0], pad_dim[1], pad_dim[2]),
                  "SetDim", pad_dim[0], pad_dim[1], pad_dim[2]);
    RequireSetter(padded, padded->SetUnitCell(
                      static_cast<float>(pad_dim[0] * src_spacing[0]),
                      static_cast<float>(pad_dim[1] * src_spacing[1]),
                      static_cast<float>(pad_dim[2] * src_spacing[2]),
                      90.0f, 90.0f, 90.0f,
                      pad_dim[0], pad_dim[1], pad_dim[2]),
                  "SetUnitCell", pad_dim[0] * src_spacing[0],
                  pad_dim[1] * src_spacing[1], pad_dim[2] * src_spacing[2]);
    RequireSetter(padded, padded->SetMid(static_cast<float>(pad_mid[0]),
                                         static_cast<float>(pad_mid[1]),
                                         static_cast<float>(pad_mid[2])),
                  "SetMid", pad_mid[0], pad_mid[1], pad_mid[2]);

    const GridParams pad_gp = get_grid_params(*padded);
    const float* src_values = grid.GetValues();
    float* out = padded->GetValues();
    const unsigned int size = padded->GetSize();

    for (unsigned int i = 0; i < size; ++i) {
        const unsigned int ix = i % pad_gp.x_dim;
        const unsigned int iy = (i / pad_gp.x_dim) % pad_gp.y_dim;
        const unsigned int iz = i / (pad_gp.x_dim * pad_gp.y_dim);
        const double sx = pad_gp.x_origin + ix * pad_gp.x_spacing;
        const double sy = pad_gp.y_origin + iy * pad_gp.y_spacing;
        const double sz = pad_gp.z_origin + iz * pad_gp.z_spacing;

        double wx = gp.x_origin + std::fmod(sx - gp.x_origin, cell_a);
        if (wx < gp.x_origin) wx += cell_a;
        double wy = gp.y_origin + std::fmod(sy - gp.y_origin, cell_b);
        if (wy < gp.y_origin) wy += cell_b;
        double wz = gp.z_origin + std::fmod(sz - gp.z_origin, cell_c);
        if (wz < gp.z_origin) wz += cell_c;

        out[i] = static_cast<float>(
            interpolate_density_at(gp, src_values, wx, wy, wz, 0.0));
    }

    return padded;
}

}  // namespace Maptitude
