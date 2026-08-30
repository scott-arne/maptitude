#include <gtest/gtest.h>
#include "maptitude/Grid.h"
#include "maptitude/GridOps.h"
#include "maptitude/Error.h"

#include <oechem.h>
#include <oegrid.h>

#include <cmath>
#include <memory>
#include <vector>

using namespace Maptitude;

// Placeholder: verify MapOp enum values compile
TEST(GridOpsTest, MapOpEnumValues) {
    EXPECT_NE(static_cast<int>(MapOp::ADD), static_cast<int>(MapOp::SUBTRACT));
    EXPECT_NE(static_cast<int>(MapOp::MIN), static_cast<int>(MapOp::MAX));
}

// --- Helper: create a small grid with a known pattern ---

/// Build a cubic skew grid: `n` nodes per axis at `spacing`, first node at
/// `node0` on all three axes.
///
/// The skew carrier has no extents-box constructor, so the geometry is set
/// explicitly. SetMid takes the grid centre, which sits (n - 1) * spacing / 2
/// above the first node. SetDim must precede any GetValues() call: the value
/// array does not exist until the dimensions are known.
static OESystem::OESkewGrid MakeCubicGrid(const unsigned int n, const double spacing,
                                          const double node0) {
    OESystem::OESkewGrid grid;
    EXPECT_TRUE(grid.SetDim(n, n, n));
    const float edge = static_cast<float>(n * spacing);
    EXPECT_TRUE(grid.SetUnitCell(edge, edge, edge, 90.0f, 90.0f, 90.0f, n, n, n));
    const float mid = static_cast<float>(node0 + (n - 1) * spacing / 2.0);
    EXPECT_TRUE(grid.SetMid(mid, mid, mid));
    return grid;
}

static OESystem::OESkewGrid MakeTestGrid() {
    // 10x10x10 grid, spacing 1.0, first node at (0,0,0)
    OESystem::OESkewGrid grid = MakeCubicGrid(10u, 1.0, 0.0);

    // Fill with a pattern: value = x + 10*y + 100*z at grid points
    float* values = grid.GetValues();
    const unsigned int size = grid.GetSize();
    for (unsigned int i = 0; i < size; ++i) {
        float x, y, z;
        grid.ElementToSpatialCoord(i, x, y, z);
        values[i] = x + 10.0f * y + 100.0f * z;
    }
    return grid;
}

// --- InterpolateDensityPeriodic tests ---

TEST(GridOpsTest, InterpolateDensityPeriodicWraps) {
    auto grid = MakeTestGrid();
    double cell_a = 10.0, cell_b = 10.0, cell_c = 10.0;

    // Query at (2, 3, 4) directly
    double val_direct = interpolate_density_periodic(
        grid, 2.0, 3.0, 4.0, cell_a, cell_b, cell_c);

    // Query at (2+10, 3+10, 4+10) should wrap to the same point
    double val_wrapped = interpolate_density_periodic(
        grid, 12.0, 13.0, 14.0, cell_a, cell_b, cell_c);

    EXPECT_NEAR(val_direct, val_wrapped, 1e-4);

    // Query at (2+20, 3+30, 4+40) — multiple periods
    double val_multi = interpolate_density_periodic(
        grid, 22.0, 33.0, 44.0, cell_a, cell_b, cell_c);

    EXPECT_NEAR(val_direct, val_multi, 1e-4);
}

TEST(GridOpsTest, InterpolateDensityPeriodicNegativeWrap) {
    auto grid = MakeTestGrid();
    double cell_a = 10.0, cell_b = 10.0, cell_c = 10.0;

    // Query at (5, 5, 5) directly
    double val_direct = interpolate_density_periodic(
        grid, 5.0, 5.0, 5.0, cell_a, cell_b, cell_c);

    // Query at (5-10, 5-10, 5-10) = (-5, -5, -5) should wrap back
    double val_neg = interpolate_density_periodic(
        grid, -5.0, -5.0, -5.0, cell_a, cell_b, cell_c);

    EXPECT_NEAR(val_direct, val_neg, 1e-4);

    // Query at (5-20, 5-30, 5-40) — multiple negative periods
    double val_multi_neg = interpolate_density_periodic(
        grid, -15.0, -25.0, -35.0, cell_a, cell_b, cell_c);

    EXPECT_NEAR(val_direct, val_multi_neg, 1e-4);
}

TEST(GridOpsTest, InterpolateDensityPeriodicBatchConsistency) {
    auto grid = MakeTestGrid();
    double cell_a = 10.0, cell_b = 10.0, cell_c = 10.0;

    std::vector<double> points = {
        2.0, 3.0, 4.0,     // in-bounds
        12.0, 13.0, 14.0,  // wrapped +1 period
        -5.0, -5.0, -5.0   // negative wrap
    };

    auto results = interpolate_density_periodic_batch(
        grid, points, 3, cell_a, cell_b, cell_c);

    ASSERT_EQ(results.size(), 3u);

    // First two should be the same (one period apart)
    EXPECT_NEAR(results[0], results[1], 1e-4);

    // Each should match single-point call
    for (size_t i = 0; i < 3; ++i) {
        double single = interpolate_density_periodic(
            grid, points[i*3], points[i*3+1], points[i*3+2],
            cell_a, cell_b, cell_c);
        EXPECT_NEAR(results[i], single, 1e-10);
    }
}

// --- WrapAndPadGrid tests ---

static OEChem::OEGraphMol MakeTestMol(double cx, double cy, double cz) {
    // Create a small molecule with a single heavy atom at (cx, cy, cz)
    OEChem::OEGraphMol mol;
    OEChem::OEAtomBase* atom = mol.NewAtom(6);  // Carbon
    float coords[3] = {
        static_cast<float>(cx),
        static_cast<float>(cy),
        static_cast<float>(cz)
    };
    mol.SetCoords(atom, coords);
    return mol;
}

TEST(GridOpsTest, WrapAndPadGridNoShiftNeeded) {
    auto grid = MakeTestGrid();
    // Molecule centroid at (5, 5, 5) — right at grid center, within padding
    auto mol = MakeTestMol(5.0, 5.0, 5.0);

    OESystem::OESkewGrid* result = wrap_and_pad_grid(
        grid, mol, 10.0, 10.0, 10.0, 3.0);

    // Atom is well within grid, no padding needed → nullptr
    EXPECT_EQ(result, nullptr);

    // Verify atom was not shifted (centroid already near grid center)
    float coords[3];
    OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms();
    mol.GetCoords(&(*atom), coords);
    EXPECT_NEAR(coords[0], 5.0, 0.5);
    EXPECT_NEAR(coords[1], 5.0, 0.5);
    EXPECT_NEAR(coords[2], 5.0, 0.5);
}

TEST(GridOpsTest, WrapAndPadGridShiftsCoordinates) {
    auto grid = MakeTestGrid();
    // Molecule at (25, 35, 45) — far from grid center, needs shifting
    auto mol = MakeTestMol(25.0, 35.0, 45.0);

    OESystem::OESkewGrid* result = wrap_and_pad_grid(
        grid, mol, 10.0, 10.0, 10.0, 3.0);

    // After shifting, atom should be near grid center
    float coords[3];
    OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms();
    mol.GetCoords(&(*atom), coords);

    // Should have been shifted by integer multiples of cell dimensions
    // Grid center is approximately (4.5, 4.5, 4.5)
    // Shift = round((4.5 - 25)/10) * 10 = round(-2.05) * 10 = -20
    // New x = 25 - 20 = 5.0
    EXPECT_NEAR(coords[0], 5.0, 0.5);
    EXPECT_NEAR(coords[1], 5.0, 0.5);
    EXPECT_NEAR(coords[2], 5.0, 0.5);

    // Since atom is within grid after shifting, result should be nullptr
    // (no padding needed)
    if (result) {
        delete result;
    }
}

TEST(GridOpsTest, WrapAndPadGridCreatesPaddedGrid) {
    // Small 5x5x5 grid, spacing 1.0, first node at (0,0,0)
    OESystem::OESkewGrid grid = MakeCubicGrid(5u, 1.0, 0.0);

    // Fill with constant value 42.0
    float* values = grid.GetValues();
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        values[i] = 42.0f;
    }

    // Molecule at (2, 2, 2) with padding 5.0 will exceed the 5x5x5 grid
    auto mol = MakeTestMol(2.0, 2.0, 2.0);

    OESystem::OESkewGrid* result = wrap_and_pad_grid(
        grid, mol, 5.0, 5.0, 5.0, 5.0);  // large padding forces pad

    // Padded grid should have been created
    ASSERT_NE(result, nullptr);

    // The padded grid should be larger than the original 5x5x5 grid
    EXPECT_GT(result->GetXDim(), grid.GetXDim());
    EXPECT_GT(result->GetYDim(), grid.GetYDim());
    EXPECT_GT(result->GetZDim(), grid.GetZDim());

    // Values in the padded grid should be ~42.0 (filled from periodic sampling)
    float sx, sy, sz;
    result->ElementToSpatialCoord(0, sx, sy, sz);
    double val = interpolate_density_periodic(
        grid, sx, sy, sz, 5.0, 5.0, 5.0);
    EXPECT_NEAR(result->GetValues()[0], static_cast<float>(val), 0.1);

    delete result;
}

static OESystem::OESkewGrid MakeShiftedGrid(double shift) {
    OESystem::OESkewGrid grid = MakeCubicGrid(10u, 1.0, shift);
    float* values = grid.GetValues();
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        values[i] = 1.0f;
    }
    return grid;
}

TEST(GridOpsTest, CombineRejectsGridsWithDifferentOrigins) {
    // Same dims, same spacing, different origin. Element-wise combination would
    // mix densities from different points in space.
    OESystem::OESkewGrid a = MakeShiftedGrid(0.0);
    OESystem::OESkewGrid b = MakeShiftedGrid(5.0);
    EXPECT_THROW(combine_maps(a, b, MapOp::ADD), GridError);
}

TEST(GridOpsTest, CombineRejectsGridsWhoseSpacingDiffersBelowTheOldTolerance) {
    // The other half of the exact-geometry change, and the reachable half: the old
    // hand-rolled comparison tested dimensions exactly but spacing only to within 1e-6.
    // 0.5 + 1e-7 lands on the next float32 value, 0.5000001192092896, a delta of
    // 1.19e-7 -- inside that tolerance and outside exact equality. Without this case,
    // restoring the tolerance reverts a documented behavior change with the suite still
    // green, because the origin tests above pass either way.
    //
    // The dimensions are given explicitly rather than derived from a bounding box. From
    // a box the two spacings yield 19 and 18 points per axis, so the old dimension test
    // would reject them and the spacing comparison would never be reached.
    //
    // The perturbation is 4e-7, not the 1e-7 this case used against OEGridSameGeometry.
    // same_grid_geometry compares per-axis spacing and node origin at a relative
    // tolerance floored at 1.0, so it accepts a spacing delta of 1.06e-7 outright and
    // the case stopped rejecting anything. The rejection available to it is the
    // origin's: both grids are centred on the same midpoint, so a spacing delta moves
    // the first node by nine half-steps and clears 1e-6 well before the spacing itself
    // does. The guarded band is spacing deltas in roughly (1.1e-7, 1e-6) -- inside the
    // old tolerance, outside the new one by way of the origin they move. 4e-7 sits with
    // margin on both sides: 4.2e-7 of spacing against a 1e-6 ceiling, 3.8e-6 of origin
    // against a 1e-6 floor.
    //
    // Both grids are built around midpoint 4.5, the way the old extents-box constructor
    // worked: it stored a midpoint and derived the first node from it.
    const double mid = 4.5;
    const double sp_a = 0.5;
    const double sp_b = static_cast<float>(0.5 + 4e-7);
    OESystem::OESkewGrid a = MakeCubicGrid(19u, sp_a, mid - 9 * sp_a);
    OESystem::OESkewGrid b = MakeCubicGrid(19u, sp_b, mid - 9 * sp_b);

    const GridParams ga = get_grid_params(a);
    const GridParams gb = get_grid_params(b);
    ASSERT_EQ(ga.x_dim, gb.x_dim) << "the dimensions must match or this pins nothing";
    ASSERT_EQ(a.GetXMid(), b.GetXMid()) << "the midpoints must match or this pins nothing";
    ASSERT_NE(ga.x_spacing, gb.x_spacing) << "the two spacings collapsed to one float";
    ASSERT_LT(std::fabs(ga.x_spacing - gb.x_spacing), 1e-6)
        << "the spacing delta must sit inside the old tolerance or this pins nothing";
    ASSERT_GT(std::fabs(ga.x_origin - gb.x_origin), 1e-6)
        << "the origin delta must sit outside the new tolerance or nothing can reject";
    EXPECT_THROW(combine_maps(a, b, MapOp::ADD), GridError);
}

TEST(GridOpsTest, CombineStillAcceptsIdenticalGeometry) {
    OESystem::OESkewGrid a = MakeTestGrid();
    OESystem::OESkewGrid b = MakeTestGrid();
    EXPECT_NO_THROW({
        std::unique_ptr<OESystem::OESkewGrid> result(combine_maps(a, b, MapOp::ADD));
        ASSERT_NE(result, nullptr);
    });
}

TEST(GridOpsTest, DiffToCalcRejectsMismatchedGeometry) {
    OESystem::OESkewGrid obs = MakeShiftedGrid(0.0);
    OESystem::OESkewGrid diff = MakeShiftedGrid(5.0);
    EXPECT_THROW(diff_to_calc(obs, diff), GridError);
}

TEST(GridOpsTest, WrapAndPadThrowsWhenTheMoleculeHasNoHeavyAtoms) {
    OEChem::OEGraphMol mol;  // empty
    OESystem::OESkewGrid grid = MakeTestGrid();
    EXPECT_THROW(wrap_and_pad_grid(grid, mol, 20.0, 25.0, 30.0), StructureError);
}

TEST(GridOpsTest, WrapAndPadThrowsForAMoleculeOfOnlyDummyAtoms) {
    OEChem::OEGraphMol mol;
    OEChem::OEAtomBase* atom = mol.NewAtom(0);  // Dummy atom (Z=0)
    const float coords[3] = {4.5f, 4.5f, 4.5f};
    mol.SetCoords(atom, coords);

    OESystem::OESkewGrid grid = MakeTestGrid();
    EXPECT_THROW(wrap_and_pad_grid(grid, mol, 20.0, 25.0, 30.0), StructureError);
}

TEST(GridOpsTest, WrapAndPadThrowsForAnAllHydrogenMolecule) {
    OEChem::OEGraphMol mol;
    OEChem::OEAtomBase* atom = mol.NewAtom(1);  // Hydrogen (Z=1)
    const float coords[3] = {4.5f, 4.5f, 4.5f};
    mol.SetCoords(atom, coords);

    OESystem::OESkewGrid grid = MakeTestGrid();
    EXPECT_THROW(wrap_and_pad_grid(grid, mol, 20.0, 25.0, 30.0), StructureError);
}

TEST(GridOpsTest, WrapAndPadReturnsNullptrOnlyWhenNoPaddingIsNeeded) {
    // A molecule already well inside the grid needs no padding: nullptr means
    // "unchanged", and nothing else.
    OEChem::OEGraphMol mol;
    OEChem::OEAtomBase* atom = mol.NewAtom(6);
    const double coords[3] = {4.5, 4.5, 4.5};
    mol.SetCoords(atom, coords);

    OESystem::OESkewGrid grid = MakeTestGrid();
    std::unique_ptr<OESystem::OESkewGrid> result(
        wrap_and_pad_grid(grid, mol, 20.0, 25.0, 30.0));
    EXPECT_EQ(result, nullptr);
}

// The extents-box constructor is the one piece of scalar-carrier geometry the
// skew carrier cannot express, so this test pins the reproduction directly
// rather than trusting it.
TEST(WrapAndPadGrid, ReproducesTheExtentsBoxConstructor) {
    struct Case {
        double minmax[6];
        double spacing;
        unsigned int dim[3];
        double mid[3];
        double node0[3];
    };
    const Case CASES[] = {
        {{0.0, 0.0, 0.0, 9.0, 9.0, 9.0}, 1.0,
         {10u, 10u, 10u}, {4.5, 4.5, 4.5}, {0.0, 0.0, 0.0}},
        {{0.0, 0.0, 0.0, 9.5, 9.5, 9.5}, 1.0,
         {10u, 10u, 10u}, {4.75, 4.75, 4.75}, {0.25, 0.25, 0.25}},
        {{-1.3, 2.7, 0.4, 8.2, 11.1, 5.9}, 0.7,
         {14u, 12u, 8u}, {3.45, 6.90, 3.15}, {-1.10, 3.05, 0.70}},
    };

    for (const Case& c : CASES) {
        // The reference this reproduction is measured against is the scalar
        // carrier's extents-box constructor itself.
        double minmax[6];
        for (int i = 0; i < 6; ++i) minmax[i] = c.minmax[i];
        const OESystem::OEScalarGrid reference(minmax, c.spacing);  // OE-SCALARGRID-OK: the reproduction's reference

        EXPECT_EQ(reference.GetXDim(), c.dim[0]);
        EXPECT_EQ(reference.GetYDim(), c.dim[1]);
        EXPECT_EQ(reference.GetZDim(), c.dim[2]);
        EXPECT_NEAR(reference.GetXMid(), c.mid[0], 1e-5);
        EXPECT_NEAR(reference.GetYMid(), c.mid[1], 1e-5);
        EXPECT_NEAR(reference.GetZMid(), c.mid[2], 1e-5);

        const OESystem::OESkewGrid converted(reference);
        const GridParams gp = get_grid_params(converted);
        EXPECT_NEAR(gp.x_origin, c.node0[0], 1e-5);
        EXPECT_NEAR(gp.y_origin, c.node0[1], 1e-5);
        EXPECT_NEAR(gp.z_origin, c.node0[2], 1e-5);
    }
}
