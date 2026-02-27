#include "maptitude/GridOps.h"
#include "maptitude/Error.h"

#include <oechem.h>
#include <oegrid.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace Maptitude {

void scale_map(OESystem::OEScalarGrid& grid, const double factor) {
    const unsigned int size = grid.GetSize();
    for (unsigned int i = 0; i < size; ++i) {
        grid[i] = static_cast<float>(grid[i] * factor);
    }
}

OESystem::OEScalarGrid* combine_maps(
    const OESystem::OEScalarGrid& lhs,
    const OESystem::OEScalarGrid& rhs,
    const MapOp op) {
    if (lhs.GetXDim() != rhs.GetXDim() ||
        lhs.GetYDim() != rhs.GetYDim() ||
        lhs.GetZDim() != rhs.GetZDim()) {
        throw GridError("Grid dimensions must match for combination");
    }

    if (std::abs(lhs.GetSpacing() - rhs.GetSpacing()) > 1e-6) {
        throw GridError("Grid spacings must match for combination");
    }

    auto* result = new OESystem::OEScalarGrid(lhs);
    const unsigned int size = result->GetSize();

    for (unsigned int i = 0; i < size; ++i) {
        const float lval = lhs[i];
        const float rval = rhs[i];
        float combined = 0.0f;

        switch (op) {
            case MapOp::ADD:
                combined = lval + rval;
                break;
            case MapOp::SUBTRACT:
                combined = lval - rval;
                break;
            case MapOp::MIN:
                combined = std::min(lval, rval);
                break;
            case MapOp::MAX:
                combined = std::max(lval, rval);
                break;
        }

        (*result)[i] = combined;
    }

    return result;
}

OESystem::OEScalarGrid* diff_to_calc(
    const OESystem::OEScalarGrid& obs_grid,
    const OESystem::OEScalarGrid& diff_grid) {
    if (obs_grid.GetXDim() != diff_grid.GetXDim() ||
        obs_grid.GetYDim() != diff_grid.GetYDim() ||
        obs_grid.GetZDim() != diff_grid.GetZDim()) {
        throw GridError("Grid dimensions must match for diff_to_calc");
    }

    auto* result = new OESystem::OEScalarGrid(obs_grid);
    const unsigned int size = result->GetSize();

    for (unsigned int i = 0; i < size; ++i) {
        // rho_calc = rho_obs - 2 * rho_diff
        const float calc = obs_grid[i] - 2.0f * diff_grid[i];
        (*result)[i] = calc;
    }

    return result;
}

OESystem::OEScalarGrid* wrap_and_pad_grid(
    const OESystem::OEScalarGrid& grid,
    OEChem::OEMolBase& mol,
    const double cell_a, const double cell_b, const double cell_c,
    const double padding) {
    const double sp = grid.GetSpacing();

    // Compute heavy-atom centroid
    double cx = 0.0, cy = 0.0, cz = 0.0;
    int n = 0;
    float coords[3];
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() == 1) continue;
        mol.GetCoords(&(*atom), coords);
        cx += coords[0];
        cy += coords[1];
        cz += coords[2];
        ++n;
    }
    if (n == 0) return nullptr;
    cx /= n;
    cy /= n;
    cz /= n;

    // Shift centroid to grid centre using integer unit-cell vectors
    const double grid_xmid = grid.GetXMin() + (grid.GetXDim() - 1) * sp / 2.0;
    const double grid_ymid = grid.GetYMin() + (grid.GetYDim() - 1) * sp / 2.0;
    const double grid_zmid = grid.GetZMin() + (grid.GetZDim() - 1) * sp / 2.0;

    const double shift_x = (cell_a > 0) ? std::round((grid_xmid - cx) / cell_a) * cell_a : 0.0;
    const double shift_y = (cell_b > 0) ? std::round((grid_ymid - cy) / cell_b) * cell_b : 0.0;
    const double shift_z = (cell_c > 0) ? std::round((grid_zmid - cz) / cell_c) * cell_c : 0.0;

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

    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() == 1) continue;
        mol.GetCoords(&(*atom), coords);
        min_x = std::min(min_x, static_cast<double>(coords[0]));
        min_y = std::min(min_y, static_cast<double>(coords[1]));
        min_z = std::min(min_z, static_cast<double>(coords[2]));
        max_x = std::max(max_x, static_cast<double>(coords[0]));
        max_y = std::max(max_y, static_cast<double>(coords[1]));
        max_z = std::max(max_z, static_cast<double>(coords[2]));
    }

    const double grid_xmin = grid.GetXMin();
    const double grid_ymin = grid.GetYMin();
    const double grid_zmin = grid.GetZMin();
    const double grid_xmax = grid_xmin + (grid.GetXDim() - 1) * sp;
    const double grid_ymax = grid_ymin + (grid.GetYDim() - 1) * sp;
    const double grid_zmax = grid_zmin + (grid.GetZDim() - 1) * sp;

    const bool needs_pad =
        (min_x - padding < grid_xmin) ||
        (max_x + padding > grid_xmax) ||
        (min_y - padding < grid_ymin) ||
        (max_y + padding > grid_ymax) ||
        (min_z - padding < grid_zmin) ||
        (max_z + padding > grid_zmax);

    if (!needs_pad) return nullptr;

    // Build a padded grid using periodic wrapping of the unit-cell density
    double minmax[6] = {
        min_x - padding, min_y - padding, min_z - padding,
        max_x + padding, max_y + padding, max_z + padding
    };
    auto* padded = new OESystem::OEScalarGrid(minmax, sp);

    const double orig_xmin = grid.GetXMin();
    const double orig_ymin = grid.GetYMin();
    const double orig_zmin = grid.GetZMin();

    const unsigned int size = padded->GetSize();
    for (unsigned int i = 0; i < size; ++i) {
        float sx, sy, sz;
        padded->ElementToSpatialCoord(i, sx, sy, sz);

        double wx = orig_xmin + std::fmod(static_cast<double>(sx) - orig_xmin, cell_a);
        if (wx < orig_xmin) wx += cell_a;

        double wy = orig_ymin + std::fmod(static_cast<double>(sy) - orig_ymin, cell_b);
        if (wy < orig_ymin) wy += cell_b;

        double wz = orig_zmin + std::fmod(static_cast<double>(sz) - orig_zmin, cell_c);
        if (wz < orig_zmin) wz += cell_c;

        (*padded)[i] = OESystem::OEFloatGridLinearInterpolate(
            grid,
            static_cast<float>(wx),
            static_cast<float>(wy),
            static_cast<float>(wz),
            0.0f);
    }

    return padded;
}

}  // namespace Maptitude
