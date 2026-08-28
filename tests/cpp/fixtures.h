#ifndef MAPTITUDE_TEST_FIXTURES_H
#define MAPTITUDE_TEST_FIXTURES_H

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include <oechem.h>
#include <oegrid.h>
#include <oesystem.h>

/// Deterministic builders for metric and grid tests.
///
/// Every builder returns a fresh object by value. Metric entry points take a
/// non-const molecule reference and PrepareStructure assigns radii in place, so
/// a shared fixture would make results depend on test execution order.
namespace MaptitudeTest {

/// Build a single-atom molecule with one residue, positioned at (x, y, z).
static inline OEChem::OEGraphMol MakeAtomMol(unsigned int atomic_num, double x, double y, double z) {
    OEChem::OEGraphMol mol;
    OEChem::OEAtomBase* atom = mol.NewAtom(atomic_num);
    const double coords[3] = {x, y, z};
    mol.SetCoords(atom, coords);

    OEChem::OEResidue residue;
    residue.SetName("LIG");
    residue.SetResidueNumber(1);
    residue.SetChainID('A');
    residue.SetBFactor(0.0);
    OEChem::OEAtomSetResidue(atom, residue);

    return mol;
}

/// Build a cubic grid spanning [-half_width, +half_width] on each axis.
static inline OESystem::OEScalarGrid MakeEmptyGrid(double half_width, double spacing) {
    double minmax[6] = {-half_width, -half_width, -half_width, half_width, half_width, half_width};
    return OESystem::OEScalarGrid(minmax, spacing);
}

/// Build a grid holding an isotropic Gaussian centred at (cx, cy, cz), peak 1.0.
static inline OESystem::OEScalarGrid MakeGaussianGrid(double cx, double cy, double cz, double sigma,
                                                      double half_width, double spacing) {
    OESystem::OEScalarGrid grid = MakeEmptyGrid(half_width, spacing);
    const double two_sigma_sq = 2.0 * sigma * sigma;
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        float gx, gy, gz;
        grid.ElementToSpatialCoord(i, gx, gy, gz);
        const double dx = gx - cx, dy = gy - cy, dz = gz - cz;
        grid[i] = static_cast<float>(std::exp(-(dx * dx + dy * dy + dz * dz) / two_sigma_sq));
    }
    return grid;
}

/// Build a grid where every element holds the same value.
static inline OESystem::OEScalarGrid MakeUniformGrid(float value, double half_width, double spacing) {
    OESystem::OEScalarGrid grid = MakeEmptyGrid(half_width, spacing);
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        grid[i] = value;
    }
    return grid;
}

/// Build a grid whose value increases linearly along x, y, and z with distinct
/// weights, so no symmetry can mask an index bug.
///
/// Values are guaranteed unique. The decade weights separate ten grid points
/// per axis; beyond that they wrap — ten steps in y offset exactly one step in
/// z — so geometries with more than ten points on any axis produce colliding
/// values and are rejected with `std::invalid_argument`. Example: (4.0, 0.5)
/// produces 17 points per axis and would yield 1777 distinct values across 4913
/// elements, so it is rejected rather than silently returned.
static inline OESystem::OEScalarGrid MakeRampGrid(double half_width, double spacing) {
    OESystem::OEScalarGrid grid = MakeEmptyGrid(half_width, spacing);
    const unsigned int max_dim = std::max({grid.GetXDim(), grid.GetYDim(), grid.GetZDim()});
    if (max_dim > 10) {
        throw std::invalid_argument(
            "MakeRampGrid: the decade weights only separate ten points per axis; "
            "this geometry produces more, so values would collide");
    }
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        float gx, gy, gz;
        grid.ElementToSpatialCoord(i, gx, gy, gz);
        grid[i] = gx + 10.0f * gy + 100.0f * gz;
    }
    return grid;
}

/// Return (grid, -grid). Used to assert that RSCC of exact anticorrelation is -1.
static inline std::pair<OESystem::OEScalarGrid, OESystem::OEScalarGrid>
MakeNegatedPair(const OESystem::OEScalarGrid& grid) {
    OESystem::OEScalarGrid negated(grid);
    for (unsigned int i = 0; i < negated.GetSize(); ++i) {
        negated[i] = -negated[i];
    }
    return {OESystem::OEScalarGrid(grid), negated};
}

}  // namespace MaptitudeTest

#endif  // MAPTITUDE_TEST_FIXTURES_H
