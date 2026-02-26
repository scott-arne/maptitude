#include <gtest/gtest.h>
#include "maptitude/Grid.h"
#include "maptitude/GridOps.h"
#include "maptitude/Error.h"

#include <oechem.h>
#include <oegrid.h>

#include <cmath>

using namespace Maptitude;

// Placeholder: verify MapOp enum values compile
TEST(GridOpsTest, MapOpEnumValues) {
    EXPECT_NE(static_cast<int>(MapOp::Add), static_cast<int>(MapOp::Subtract));
    EXPECT_NE(static_cast<int>(MapOp::Min), static_cast<int>(MapOp::Max));
}

// --- Helper: create a small grid with a known pattern ---

static OESystem::OEScalarGrid MakeTestGrid() {
    // 10x10x10 grid, spacing 1.0, origin at (0,0,0)
    // minmax: [0, 0, 0, 9, 9, 9]
    double minmax[6] = {0.0, 0.0, 0.0, 9.0, 9.0, 9.0};
    OESystem::OEScalarGrid grid(minmax, 1.0);

    // Fill with a pattern: value = x + 10*y + 100*z at grid points
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        float x, y, z;
        grid.ElementToSpatialCoord(i, x, y, z);
        grid[i] = x + 10.0f * y + 100.0f * z;
    }
    return grid;
}

// --- InterpolateDensityPeriodic tests ---

TEST(GridOpsTest, InterpolateDensityPeriodicWraps) {
    auto grid = MakeTestGrid();
    double cell_a = 10.0, cell_b = 10.0, cell_c = 10.0;

    // Query at (2, 3, 4) directly
    double val_direct = InterpolateDensityPeriodic(
        grid, 2.0, 3.0, 4.0, cell_a, cell_b, cell_c);

    // Query at (2+10, 3+10, 4+10) should wrap to the same point
    double val_wrapped = InterpolateDensityPeriodic(
        grid, 12.0, 13.0, 14.0, cell_a, cell_b, cell_c);

    EXPECT_NEAR(val_direct, val_wrapped, 1e-4);

    // Query at (2+20, 3+30, 4+40) — multiple periods
    double val_multi = InterpolateDensityPeriodic(
        grid, 22.0, 33.0, 44.0, cell_a, cell_b, cell_c);

    EXPECT_NEAR(val_direct, val_multi, 1e-4);
}

TEST(GridOpsTest, InterpolateDensityPeriodicNegativeWrap) {
    auto grid = MakeTestGrid();
    double cell_a = 10.0, cell_b = 10.0, cell_c = 10.0;

    // Query at (5, 5, 5) directly
    double val_direct = InterpolateDensityPeriodic(
        grid, 5.0, 5.0, 5.0, cell_a, cell_b, cell_c);

    // Query at (5-10, 5-10, 5-10) = (-5, -5, -5) should wrap back
    double val_neg = InterpolateDensityPeriodic(
        grid, -5.0, -5.0, -5.0, cell_a, cell_b, cell_c);

    EXPECT_NEAR(val_direct, val_neg, 1e-4);

    // Query at (5-20, 5-30, 5-40) — multiple negative periods
    double val_multi_neg = InterpolateDensityPeriodic(
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

    auto results = InterpolateDensityPeriodicBatch(
        grid, points, 3, cell_a, cell_b, cell_c);

    ASSERT_EQ(results.size(), 3u);

    // First two should be the same (one period apart)
    EXPECT_NEAR(results[0], results[1], 1e-4);

    // Each should match single-point call
    for (size_t i = 0; i < 3; ++i) {
        double single = InterpolateDensityPeriodic(
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

    OESystem::OEScalarGrid* result = WrapAndPadGrid(
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

    OESystem::OEScalarGrid* result = WrapAndPadGrid(
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
    // Small 5x5x5 grid, spacing 1.0, origin at (0,0,0)
    double minmax[6] = {0.0, 0.0, 0.0, 4.0, 4.0, 4.0};
    OESystem::OEScalarGrid grid(minmax, 1.0);

    // Fill with constant value 42.0
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        grid[i] = 42.0f;
    }

    // Molecule at (2, 2, 2) with padding 5.0 will exceed the 5x5x5 grid
    auto mol = MakeTestMol(2.0, 2.0, 2.0);

    OESystem::OEScalarGrid* result = WrapAndPadGrid(
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
    double val = InterpolateDensityPeriodic(
        grid, sx, sy, sz, 5.0, 5.0, 5.0);
    EXPECT_NEAR((*result)[0], static_cast<float>(val), 0.1);

    delete result;
}
