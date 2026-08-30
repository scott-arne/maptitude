#include <gtest/gtest.h>
#include "maptitude/Grid.h"
#include "maptitude/GridOps.h"
#include "maptitude/Error.h"

#include <oechem.h>
#include <oegrid.h>

#include <algorithm>
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

TEST(GridOpsTest, InterpolateDensityPeriodicReturnsAnAbsoluteValueAtAWrappedNode) {
    // The three tests above each compare one wrapped query against another, so any
    // wrap anchor satisfies them: shifting the anchor shifts both sides equally.
    // MakeTestGrid stores x + 10y + 100z at every node, which gives a wrapped query
    // one correct answer and pins where the anchor actually is.
    auto grid = MakeTestGrid();
    EXPECT_NEAR(interpolate_density_periodic(grid, 12.0, 3.0, 4.0, 10.0, 10.0, 10.0),
                432.0, 1e-4);
    EXPECT_NEAR(interpolate_density_periodic(grid, 2.0, 3.0, 4.0, 10.0, 10.0, 10.0),
                432.0, 1e-4);
}

TEST(GridOpsTest, InterpolateDensityPeriodicBlendsAcrossTheCellBoundary) {
    // The node span is [origin, origin + (n-1)s] but the cell is n*s wide, so the
    // final spacing of every cell has no upper node of its own. It is ordinary
    // interior space: its upper neighbour is node 0 of the next image, and a
    // periodic query there must blend, not return the caller's default.
    auto grid = MakeTestGrid();
    constexpr double DEFAULT = -99.0;
    // At y = 3, z = 4 node 9 holds 439 and node 0 holds 430.
    EXPECT_NEAR(interpolate_density_periodic(grid, 9.5, 3.0, 4.0, 10.0, 10.0, 10.0, DEFAULT),
                434.5, 1e-4);
    EXPECT_NEAR(interpolate_density_periodic(grid, 9.2, 3.0, 4.0, 10.0, 10.0, 10.0, DEFAULT),
                437.2, 1e-4);
    // The origin's own periodic image, and a point half a spacing below it.
    EXPECT_NEAR(interpolate_density_periodic(grid, 10.0, 3.0, 4.0, 10.0, 10.0, 10.0, DEFAULT),
                430.0, 1e-4);
    EXPECT_NEAR(interpolate_density_periodic(grid, -0.5, 3.0, 4.0, 10.0, 10.0, 10.0, DEFAULT),
                434.5, 1e-4);
}

TEST(GridOpsTest, InterpolateDensityPeriodicRejectsAnIncommensurateCell) {
    // Wrapping modulo the cell only lands on the sampled lattice when the cell is
    // the sampled extent. A cell that is not n*s makes node n-1's periodic
    // neighbour something other than node 0, and no wrap can recover the density
    // that was never sampled.
    auto grid = MakeTestGrid();  // 10 nodes at spacing 1.0, so the cell must be 10.
    EXPECT_THROW(interpolate_density_periodic(grid, 2.0, 3.0, 4.0, 12.0, 10.0, 10.0), CellError);
    EXPECT_THROW(interpolate_density_periodic(grid, 2.0, 3.0, 4.0, 10.0, 12.0, 10.0), CellError);
    EXPECT_THROW(interpolate_density_periodic(grid, 2.0, 3.0, 4.0, 10.0, 10.0, 12.0), CellError);

    const std::vector<double> points = {2.0, 3.0, 4.0};
    EXPECT_THROW(interpolate_density_periodic_batch(grid, points, 1, 12.0, 10.0, 10.0), CellError);
    EXPECT_NO_THROW(interpolate_density_periodic_batch(grid, points, 1, 10.0, 10.0, 10.0));
}

TEST(GridOpsTest, InterpolateDensityAcceptsTheNominalNodeOrigin) {
    // The node origin is derived from ElementToSpatialCoord, and on a 90-degree
    // cell the fractional-to-Cartesian matrix carries a cos(90 deg) ~ 6e-17 term,
    // so a grid whose first node is nominally at 0.0 reports it a femtometre above.
    // Without a boundary tolerance the caller's own construction coordinate falls
    // outside the node span while a coordinate round-tripped through the carrier
    // does not.
    OESystem::OESkewGrid grid = MakeCubicGrid(5u, 1.0, 0.0);
    float* values = grid.GetValues();
    for (unsigned int i = 0; i < grid.GetSize(); ++i) values[i] = 42.0f;

    const GridParams gp = get_grid_params(grid);
    ASSERT_GT(gp.x_origin, 0.0)
        << "the derived origin no longer sits above its nominal 0.0, so this pins nothing";

    EXPECT_DOUBLE_EQ(interpolate_density(grid, 0.0, 0.0, 0.0, -99.0), 42.0);
    EXPECT_DOUBLE_EQ(interpolate_density(grid, 4.0, 4.0, 4.0, -99.0), 42.0);
    EXPECT_DOUBLE_EQ(interpolate_density_periodic(grid, 0.0, 0.0, 0.0, 5.0, 5.0, 5.0, -99.0),
                     42.0);
}

// --- grid_to_vector / vector_to_grid ---

TEST(GridOpsTest, GridToVectorAndBackRoundTripsEveryElement) {
    OESystem::OESkewGrid grid = MakeTestGrid();
    const std::vector<double> values = grid_to_vector(grid);
    ASSERT_EQ(values.size(), grid.GetSize());

    std::vector<double> doubled(values.size());
    for (size_t i = 0; i < values.size(); ++i) doubled[i] = values[i] * 2.0;

    vector_to_grid(doubled, grid);
    const float* out = grid.GetValues();
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        EXPECT_FLOAT_EQ(out[i], static_cast<float>(values[i] * 2.0)) << "element " << i;
    }
}

TEST(GridOpsTest, VectorToGridRejectsAVectorOfTheWrongLength) {
    // Copying the shorter of the two left the tail of the grid holding whatever
    // density it held before, which reads as a successful write.
    OESystem::OESkewGrid grid = MakeTestGrid();
    const unsigned int size = grid.GetSize();

    EXPECT_THROW(vector_to_grid(std::vector<double>(size - 1u, 1.0), grid), GridError);
    EXPECT_THROW(vector_to_grid(std::vector<double>(size + 1u, 1.0), grid), GridError);
    EXPECT_THROW(vector_to_grid(std::vector<double>(), grid), GridError);
    EXPECT_NO_THROW(vector_to_grid(std::vector<double>(size, 1.0), grid));
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

    // The shifted atom is inside the grid with room for the padding, so no padded
    // grid is built and there is nothing to own.
    EXPECT_EQ(result, nullptr);
    delete result;
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

    // The source is uniform, so every padded node -- interior, wrapped, or on the
    // cell face -- must carry the source value. Recomputing the expectation with
    // interpolate_density_periodic would put the function under test on both sides
    // of the comparison, and checking element 0 alone left a third of the fill
    // unexamined.
    const float* padded_values = result->GetValues();
    unsigned int wrong = 0u;
    float first_wrong = 0.0f;
    unsigned int first_wrong_index = 0u;
    for (unsigned int i = 0; i < result->GetSize(); ++i) {
        if (std::fabs(padded_values[i] - 42.0f) > 1e-3f) {
            if (wrong == 0u) {
                first_wrong = padded_values[i];
                first_wrong_index = i;
            }
            ++wrong;
        }
    }
    EXPECT_EQ(wrong, 0u) << wrong << " of " << result->GetSize()
                         << " padded nodes disagree with the uniform source; element "
                         << first_wrong_index << " holds " << first_wrong;

    delete result;
}

/// Build a skew grid with a distinct node interval on each axis, first node at the
/// Cartesian origin, filled with @p fill.
static OESystem::OESkewGrid MakeAnisotropicGrid(const unsigned int nx, const unsigned int ny,
                                                const unsigned int nz, const double sx,
                                                const double sy, const double sz,
                                                const float fill) {
    OESystem::OESkewGrid grid;
    EXPECT_TRUE(grid.SetDim(nx, ny, nz));
    EXPECT_TRUE(grid.SetUnitCell(static_cast<float>(nx * sx), static_cast<float>(ny * sy),
                                 static_cast<float>(nz * sz), 90.0f, 90.0f, 90.0f, nx, ny, nz));
    EXPECT_TRUE(grid.SetMid(static_cast<float>((nx - 1) * sx / 2.0),
                            static_cast<float>((ny - 1) * sy / 2.0),
                            static_cast<float>((nz - 1) * sz / 2.0)));
    float* values = grid.GetValues();
    for (unsigned int i = 0; i < grid.GetSize(); ++i) values[i] = fill;
    return grid;
}

TEST(GridOpsTest, WrapAndPadGridSizesEachAxisFromItsOwnNodeInterval) {
    // Per-axis sampling is the reason this carrier exists, and every other padding
    // test uses an isotropic grid, where collapsing all three intervals onto x's is
    // invisible. Distinct intervals give the three axes three different node counts
    // from one common extent.
    OESystem::OESkewGrid grid = MakeAnisotropicGrid(8u, 6u, 5u, 0.5, 1.0, 2.0, 7.0f);
    const GridParams gp = get_grid_params(grid);
    ASSERT_NEAR(gp.x_spacing, 0.5, 1e-5);
    ASSERT_NEAR(gp.y_spacing, 1.0, 1e-5);
    ASSERT_NEAR(gp.z_spacing, 2.0, 1e-5);

    // The atom sits at the grid centre, so no centroid shift runs; padding 2.0 puts
    // it outside the node span on x and forces the padding path.
    auto mol = MakeTestMol(1.75, 2.5, 4.0);
    std::unique_ptr<OESystem::OESkewGrid> result(
        wrap_and_pad_grid(grid, mol, 4.0, 6.0, 10.0, 2.0));
    ASSERT_NE(result, nullptr) << "expected the padding path, not the nullptr shortcut";

    // A 4.0 A extent on every axis, divided by that axis's own node interval:
    // ceil(4/0.5) + 1 = 9, ceil(4/1) + 1 = 5, ceil(4/2) + 1 = 3.
    EXPECT_EQ(result->GetXDim(), 9u);
    EXPECT_EQ(result->GetYDim(), 5u);
    EXPECT_EQ(result->GetZDim(), 3u);
    EXPECT_NEAR(result->GetXMid(), 1.75, 1e-5);
    EXPECT_NEAR(result->GetYMid(), 2.5, 1e-5);
    EXPECT_NEAR(result->GetZMid(), 4.0, 1e-5);

    const GridParams pad_gp = get_grid_params(*result);
    EXPECT_NEAR(pad_gp.x_spacing, 0.5, 1e-5);
    EXPECT_NEAR(pad_gp.y_spacing, 1.0, 1e-5);
    EXPECT_NEAR(pad_gp.z_spacing, 2.0, 1e-5);
    EXPECT_NEAR(pad_gp.x_origin, -0.25, 1e-5);
    EXPECT_NEAR(pad_gp.y_origin, 0.5, 1e-5);
    EXPECT_NEAR(pad_gp.z_origin, 2.0, 1e-5);

    const float* padded_values = result->GetValues();
    for (unsigned int i = 0; i < result->GetSize(); ++i) {
        EXPECT_NEAR(padded_values[i], 7.0f, 1e-3f) << "element " << i;
    }
}

TEST(GridOpsTest, WrapAndPadGridRejectsAPaddingTooSmallForTheNodeInterval) {
    // A padded axis needs two nodes before any interval can be derived from it. The
    // shortfall is the caller's padding against this grid's node interval, and the
    // message has to say so: deriving geometry from the one-node result blames the
    // source grid instead.
    OESystem::OESkewGrid grid = MakeTestGrid();  // spacing 1.0, node span [0, 9]
    auto mol = MakeTestMol(9.4, 4.5, 4.5);       // outside the span, so padding runs

    try {
        std::unique_ptr<OESystem::OESkewGrid> result(
            wrap_and_pad_grid(grid, mol, 10.0, 10.0, 10.0, 0.0));
        FAIL() << "expected a zero padding around a single atom to be rejected";
    } catch (const GridError& e) {
        const std::string message = e.what();
        EXPECT_NE(message.find("padding"), std::string::npos) << message;
        EXPECT_NE(message.find("node interval"), std::string::npos) << message;
    }
}

TEST(GridOpsTest, WrapAndPadGridRejectsACellThatIsNotTheSampledExtent) {
    // The padded grid is filled by periodic sampling, so it inherits the periodic
    // path's precondition: the cell has to be the extent the grid samples.
    OESystem::OESkewGrid grid = MakeTestGrid();  // 10 nodes at spacing 1.0
    auto mol = MakeTestMol(4.5, 4.5, 4.5);
    EXPECT_THROW(wrap_and_pad_grid(grid, mol, 12.0, 10.0, 10.0, 4.75), CellError);
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

/// Add a heavy atom at (x, y, z) to @p mol.
static void AddCarbon(OEChem::OEGraphMol& mol, const double x, const double y, const double z) {
    OEChem::OEAtomBase* atom = mol.NewAtom(6);
    const float coords[3] = {
        static_cast<float>(x), static_cast<float>(y), static_cast<float>(z)
    };
    mol.SetCoords(atom, coords);
}

// The padded grid's geometry is built by hand -- the skew carrier has no
// extents-box constructor -- so it is asserted through wrap_and_pad_grid rather
// than by re-deriving it. Asserting against a separately constructed reference
// grid would pin the reference's behaviour and leave this arithmetic free.
TEST(WrapAndPadGrid, SizesThePaddedGridFromTheAtomExtent) {
    struct Case {
        const char* label;
        double atom_lo[3];   ///< heavy atom at the low corner of the atom extent
        double atom_hi[3];   ///< heavy atom at the high corner
        double padding;
        unsigned int dim[3];
        double mid[3];
        double node0[3];
    };
    // Every case runs against MakeTestGrid: 10 nodes at spacing 1.0, node span
    // [0, 9], centre 4.5, cell 10. Each atom pair keeps its centroid inside half a
    // cell of that centre, so no coordinate shift runs and the extent is the pair's
    // bounding box grown by the padding.
    const Case CASES[] = {
        // A whole number of node intervals: 10.0 A -> 11 nodes spanning 10.0 A.
        {"integral extent",
         {4.5, 4.5, 4.5}, {4.5, 4.5, 4.5}, 5.0,
         {11u, 11u, 11u}, {4.5, 4.5, 4.5}, {-0.5, -0.5, -0.5}},
        // Nine and a half node intervals. A truncating count gives 10 nodes
        // spanning 9.0 A, a quarter of an Angstrom short at each face.
        {"fractional extent",
         {4.5, 4.5, 4.5}, {4.5, 4.5, 4.5}, 4.75,
         {11u, 11u, 11u}, {4.5, 4.5, 4.5}, {-0.5, -0.5, -0.5}},
        // A different extent on each axis, so no axis can borrow another's count.
        {"per-axis extents",
         {1.0, 2.0, 3.0}, {8.0, 6.0, 4.0}, 3.0,
         {14u, 11u, 8u}, {4.5, 4.0, 3.5}, {-2.0, -1.0, 0.0}},
    };

    for (const Case& c : CASES) {
        SCOPED_TRACE(c.label);
        OESystem::OESkewGrid grid = MakeTestGrid();
        OEChem::OEGraphMol mol;
        AddCarbon(mol, c.atom_lo[0], c.atom_lo[1], c.atom_lo[2]);
        AddCarbon(mol, c.atom_hi[0], c.atom_hi[1], c.atom_hi[2]);

        std::unique_ptr<OESystem::OESkewGrid> padded(
            wrap_and_pad_grid(grid, mol, 10.0, 10.0, 10.0, c.padding));
        ASSERT_NE(padded, nullptr) << "expected the padding path, not the nullptr shortcut";

        EXPECT_EQ(padded->GetXDim(), c.dim[0]);
        EXPECT_EQ(padded->GetYDim(), c.dim[1]);
        EXPECT_EQ(padded->GetZDim(), c.dim[2]);
        EXPECT_NEAR(padded->GetXMid(), c.mid[0], 1e-5);
        EXPECT_NEAR(padded->GetYMid(), c.mid[1], 1e-5);
        EXPECT_NEAR(padded->GetZMid(), c.mid[2], 1e-5);

        const GridParams gp = get_grid_params(*padded);
        EXPECT_NEAR(gp.x_origin, c.node0[0], 1e-5);
        EXPECT_NEAR(gp.y_origin, c.node0[1], 1e-5);
        EXPECT_NEAR(gp.z_origin, c.node0[2], 1e-5);

        // The property the node count exists to satisfy: the padded node span
        // contains every atom with its padding. A truncated count breaks this on
        // the fractional case while still producing a plausible grid.
        const double origin[3] = {gp.x_origin, gp.y_origin, gp.z_origin};
        const double spacing[3] = {gp.x_spacing, gp.y_spacing, gp.z_spacing};
        const unsigned int dim[3] = {gp.x_dim, gp.y_dim, gp.z_dim};
        for (int i = 0; i < 3; ++i) {
            EXPECT_LE(origin[i], std::min(c.atom_lo[i], c.atom_hi[i]) - c.padding + 1e-9)
                << "axis " << i << " node span starts inside the required extent";
            EXPECT_GE(origin[i] + (dim[i] - 1) * spacing[i],
                      std::max(c.atom_lo[i], c.atom_hi[i]) + c.padding - 1e-9)
                << "axis " << i << " node span ends inside the required extent";
        }
    }
}
