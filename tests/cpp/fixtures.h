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

/// Build a cubic grid of `lround(2 * half_width / spacing) + 1` nodes per axis,
/// the first at `-half_width` and the last at `-half_width + (n - 1) * spacing`.
///
/// That last node is `+half_width` only when the spacing divides `2 *
/// half_width`. Every half_width/spacing pair reaching this builder under
/// tests/cpp satisfies that except one, passed to
/// `MakeGaussianGrid(..., HALF_WIDTH, c.spacing)` at
/// test_metric_analytic.cpp:204 with `HALF_WIDTH = 4.0` and `c.spacing = 3.5`:
/// it gives 3 nodes running -4.0 to 3.0, so the far face is never reached. That
/// geometry is deliberate -- the test asserts that an adaptive Q-score sweep at
/// a step this coarse raises -- so do not "fix" it here.
///
/// The skew carrier has no extents-box constructor, so the geometry is set
/// explicitly, as `MakeCubicGrid` in test_grid_ops.cpp does. The node count is
/// the interval count plus one, and SetMid takes the grid centre rather than a
/// corner; for a span symmetric about the origin that centre is the origin.
/// SetDim must precede any GetValues() call: the value array does not exist
/// until the dimensions are known.
static inline OESystem::OESkewGrid MakeEmptyGrid(double half_width, double spacing) {
    const unsigned int n =
        static_cast<unsigned int>(std::lround(2.0 * half_width / spacing)) + 1u;
    OESystem::OESkewGrid grid;
    const float edge = static_cast<float>(n * spacing);
    const float mid = static_cast<float>(-half_width + (n - 1) * spacing / 2.0);
    // A silently rejected geometry call would leave a grid that still reads
    // plausibly but sits somewhere else, which would move every pin at once with
    // nothing pointing at the cause. Fail loudly instead.
    if (!grid.SetDim(n, n, n) || !grid.SetUnitCell(edge, edge, edge, 90.0f, 90.0f, 90.0f, n, n, n) ||
        !grid.SetMid(mid, mid, mid)) {
        throw std::invalid_argument("MakeEmptyGrid: the skew carrier rejected the geometry");
    }
    // The extents-box constructor this replaced returned a zeroed grid, and
    // MakeEmptyGrid's callers rely on that; do not depend on SetDim's allocation
    // happening to zero.
    std::fill_n(grid.GetValues(), grid.GetSize(), 0.0f);
    return grid;
}

/// Build a grid holding an isotropic Gaussian centred at (cx, cy, cz), peak 1.0.
static inline OESystem::OESkewGrid MakeGaussianGrid(double cx, double cy, double cz, double sigma,
                                                    double half_width, double spacing) {
    OESystem::OESkewGrid grid = MakeEmptyGrid(half_width, spacing);
    const double two_sigma_sq = 2.0 * sigma * sigma;
    float* values = grid.GetValues();
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        float gx, gy, gz;
        grid.ElementToSpatialCoord(i, gx, gy, gz);
        const double dx = gx - cx, dy = gy - cy, dz = gz - cz;
        values[i] = static_cast<float>(std::exp(-(dx * dx + dy * dy + dz * dz) / two_sigma_sq));
    }
    return grid;
}

/// Build a grid where every element holds the same value.
static inline OESystem::OESkewGrid MakeUniformGrid(float value, double half_width, double spacing) {
    OESystem::OESkewGrid grid = MakeEmptyGrid(half_width, spacing);
    float* values = grid.GetValues();
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        values[i] = value;
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
static inline OESystem::OESkewGrid MakeRampGrid(double half_width, double spacing) {
    OESystem::OESkewGrid grid = MakeEmptyGrid(half_width, spacing);
    const unsigned int max_dim = std::max({grid.GetXDim(), grid.GetYDim(), grid.GetZDim()});
    if (max_dim > 10) {
        throw std::invalid_argument(
            "MakeRampGrid: the decade weights only separate ten points per axis; "
            "this geometry produces more, so values would collide");
    }
    float* values = grid.GetValues();
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        float gx, gy, gz;
        grid.ElementToSpatialCoord(i, gx, gy, gz);
        values[i] = gx + 10.0f * gy + 100.0f * gz;
    }
    return grid;
}

/// Build a grid holding an isotropic Gaussian centred at (cx, cy, cz), peak 1.0, on a
/// lattice whose node interval differs from axis to axis.
///
/// The builders above are all cubic, where one scalar spacing and each axis's own node
/// interval are the same number, so nothing they produce can distinguish per-axis sampling
/// from scalar sampling. This one can.
///
/// `SetUnitCell` takes one cubic edge and the three sample counts. For the parameter sets
/// used under tests/cpp that makes the node interval on axis i `edge / n_i`. Rather than
/// rest on that arithmetic, each caller checks the intervals it depends on: two read all
/// three back through `get_grid_params` and assert them, and the third reads the interval
/// it names out of the error message it expects.
///
/// The geometry is set explicitly and the setters are checked for the same reasons as in
/// `MakeEmptyGrid` above, including throwing rather than using `ASSERT_*`.
static inline OESystem::OESkewGrid MakeAnisotropicGaussianGrid(
    double cx, double cy, double cz, double sigma, double edge,
    unsigned int nx, unsigned int ny, unsigned int nz) {
    OESystem::OESkewGrid grid;
    const float cell_edge = static_cast<float>(edge);
    if (!grid.SetDim(nx, ny, nz) ||
        !grid.SetUnitCell(cell_edge, cell_edge, cell_edge, 90.0f, 90.0f, 90.0f, nx, ny, nz) ||
        !grid.SetMid(0.0f, 0.0f, 0.0f)) {
        throw std::invalid_argument(
            "MakeAnisotropicGaussianGrid: the skew carrier rejected the geometry");
    }
    const double two_sigma_sq = 2.0 * sigma * sigma;
    float* values = grid.GetValues();
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        float gx, gy, gz;
        grid.ElementToSpatialCoord(i, gx, gy, gz);
        const double dx = gx - cx, dy = gy - cy, dz = gz - cz;
        values[i] = static_cast<float>(std::exp(-(dx * dx + dy * dy + dz * dz) / two_sigma_sq));
    }
    return grid;
}

/// Return (grid, -grid). Used to assert that RSCC of exact anticorrelation is -1.
static inline std::pair<OESystem::OESkewGrid, OESystem::OESkewGrid>
MakeNegatedPair(const OESystem::OESkewGrid& grid) {
    OESystem::OESkewGrid negated(grid);
    float* values = negated.GetValues();
    for (unsigned int i = 0; i < negated.GetSize(); ++i) {
        values[i] = -values[i];
    }
    return {OESystem::OESkewGrid(grid), negated};
}

}  // namespace MaptitudeTest

#endif  // MAPTITUDE_TEST_FIXTURES_H
