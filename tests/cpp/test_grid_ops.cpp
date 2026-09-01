#include <gtest/gtest.h>
#include "maptitude/Grid.h"
#include "maptitude/GridOps.h"
#include "maptitude/Error.h"

#include <oechem.h>
#include <oegrid.h>
#include <oesystem.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
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
    // Wrapping only lands on the sampled lattice when the cell is a whole number
    // of node intervals, and that number is one a period could span: the node
    // count, or one less on a grid whose closing plane duplicates its first. A
    // cell of twelve intervals on a ten-node grid is neither of those, and no
    // wrap can recover density that was never sampled.
    auto grid = MakeTestGrid();  // 10 nodes at spacing 1.0, so 10 or 9 intervals.
    EXPECT_THROW(interpolate_density_periodic(grid, 2.0, 3.0, 4.0, 12.0, 10.0, 10.0), CellError);
    EXPECT_THROW(interpolate_density_periodic(grid, 2.0, 3.0, 4.0, 10.0, 12.0, 10.0), CellError);
    EXPECT_THROW(interpolate_density_periodic(grid, 2.0, 3.0, 4.0, 10.0, 10.0, 12.0), CellError);

    const std::vector<double> points = {2.0, 3.0, 4.0};
    EXPECT_THROW(interpolate_density_periodic_batch(grid, points, 1, 12.0, 10.0, 10.0), CellError);
    EXPECT_NO_THROW(interpolate_density_periodic_batch(grid, points, 1, 10.0, 10.0, 10.0));
}

TEST(GridOpsTest, InterpolateDensityPeriodicTakesThePeriodFromTheCellNotTheNodeCount) {
    // Two interval counts are accepted on any grid and only one of them describes
    // it, so the cell is what chooses between them. MakeTestGrid holds ten
    // distinct nodes and no duplicated closing plane; a caller passing 9.0 is
    // asserting otherwise, and the guard takes that assertion rather than
    // inspecting the values to overrule it. Comparing the end planes would cost
    // O(n^2) per call on the hot path, and the batch entry point would pay it for
    // a property of the grid rather than of the batch.
    auto grid = MakeTestGrid();  // 10 nodes at spacing 1.0, values x + 10y + 100z.

    EXPECT_NO_THROW(interpolate_density_periodic(grid, 2.0, 3.0, 4.0, 10.0, 10.0, 10.0));
    EXPECT_NO_THROW(interpolate_density_periodic(grid, 2.0, 3.0, 4.0, 9.0, 9.0, 9.0));

    // One interval to either side of the pair is a lattice the grid was never
    // sampled on, and half an interval off rounds to an accepted count and then
    // misses it by far more than the float allowance.
    EXPECT_THROW(interpolate_density_periodic(grid, 2.0, 3.0, 4.0, 11.0, 10.0, 10.0), CellError);
    EXPECT_THROW(interpolate_density_periodic(grid, 2.0, 3.0, 4.0, 10.0, 8.0, 10.0), CellError);
    EXPECT_THROW(interpolate_density_periodic(grid, 2.0, 3.0, 4.0, 10.0, 10.0, 9.5), CellError);

    // The period the cell selects is the one the wrap runs on, so the two accepted
    // cells give different answers at the same point rather than the check being
    // cosmetic: at y = 3, z = 4 node 9 holds 439 and node 0 holds 430, and x = 9.0
    // is node 9 under a ten-interval period and node 0's image under a nine.
    EXPECT_NEAR(interpolate_density_periodic(grid, 9.0, 3.0, 4.0, 10.0, 10.0, 10.0),
                439.0, 1e-4);
    EXPECT_NEAR(interpolate_density_periodic(grid, 9.0, 3.0, 4.0, 9.0, 10.0, 10.0),
                430.0, 1e-4);
}

// --- Periodic interpolation on grids that came off a reader ---
//
// These are the suite's only periodic coverage of a grid that came out of a
// reader rather than a constructor: no other file in tests/cpp calls OEReadGrid
// and a periodic entry point. A grid MakeCubicGrid builds holds n distinct nodes
// and is passed a cell of n * spacing, and that pairing cannot show the closing
// plane OpenEye's reader appends, which is why the defect these cases pin went
// unseen.

namespace {

const char* const MAP_ASSETS[] = {
    "1d26_2fofc.ccp4", "340d_2fofc.ccp4", "390_emd_30342_A_z4.mrc", "3q9g_2fofc.ccp4"};

OESystem::OESkewGrid ReadMapAsset(const char* const name) {
    OESystem::OESkewGrid grid;
    const std::string path = std::string(MAPTITUDE_TEST_ASSET_DIR) + "/" + name;
    EXPECT_TRUE(OESystem::OEReadGrid(path, grid)) << "failed to read " << path;
    return grid;
}

}  // namespace

TEST(GridOpsTest, InterpolateDensityPeriodicAcceptsAMapReadGridsOwnCell) {
    // Each of these four stores exactly its MX,MY,MZ sections at NxSTART = 0 and
    // comes back with M + 1 nodes per axis, so the header cell each carries is
    // (n - 1) * spacing on all three axes. Measured against the n * spacing rule
    // this fix replaces, all four were refused their own cell.
    for (const char* const name : MAP_ASSETS) {
        SCOPED_TRACE(name);
        OESystem::OESkewGrid grid = ReadMapAsset(name);
        const GridParams gp = get_grid_params(grid);
        const UnitCellParams uc = get_unit_cell(grid);

        const double cell[3] = {uc.a, uc.b, uc.c};
        const unsigned int n[3] = {gp.x_dim, gp.y_dim, gp.z_dim};
        const double spacing[3] = {gp.x_spacing, gp.y_spacing, gp.z_spacing};

        // The premise. Without it an asset whose cell happened to be n * spacing
        // would pass here on the arm the old rule already accepted, and the case
        // would pin nothing.
        for (int i = 0; i < 3; ++i) {
            SCOPED_TRACE(i);
            EXPECT_LT(std::abs(cell[i] - (n[i] - 1u) * spacing[i]), 1e-4)
                << "this asset no longer reports the short cell";
            EXPECT_GT(std::abs(cell[i] - n[i] * spacing[i]), 0.5 * spacing[i])
                << "this asset's cell is now within rounding of n * spacing, which the "
                   "old rule already accepted";
        }

        EXPECT_NO_THROW(interpolate_density_periodic(
            grid, gp.x_origin, gp.y_origin, gp.z_origin, cell[0], cell[1], cell[2]));
        const std::vector<double> points = {gp.x_origin, gp.y_origin, gp.z_origin};
        EXPECT_NO_THROW(interpolate_density_periodic_batch(
            grid, points, 1, cell[0], cell[1], cell[2]));
    }
}

TEST(GridOpsTest, InterpolateDensityPeriodicDoesNotFlattenAMapReadGridsFinalInterval) {
    // Correcting the guard alone would not have been enough: the period is what
    // the wrap runs on. With a period of n_i the final interval blends node
    // n_i - 1 into node 0, and on a grid whose closing plane duplicates plane 0
    // those are the same values, so the interval comes back constant. Both cells
    // are accepted on this grid, which is what makes the wrong one callable here
    // as the control.
    OESystem::OESkewGrid grid = ReadMapAsset("1d26_2fofc.ccp4");
    const GridParams gp = get_grid_params(grid);
    const UnitCellParams uc = get_unit_cell(grid);
    const double map_cell[3] = {uc.a, uc.b, uc.c};
    const unsigned int n[3] = {gp.x_dim, gp.y_dim, gp.z_dim};
    const double spacing[3] = {gp.x_spacing, gp.y_spacing, gp.z_spacing};
    const double origin[3] = {gp.x_origin, gp.y_origin, gp.z_origin};
    const double full_cell[3] = {
        n[0] * spacing[0], n[1] * spacing[1], n[2] * spacing[2]};

    constexpr int SAMPLES = 21;
    for (int axis = 0; axis < 3; ++axis) {
        SCOPED_TRACE(axis);
        double map_lo = std::numeric_limits<double>::max();
        double map_hi = std::numeric_limits<double>::lowest();
        double full_lo = map_lo;
        double full_hi = map_hi;

        for (int k = 0; k < SAMPLES; ++k) {
            const double frac = k / static_cast<double>(SAMPLES);  // [0, 1)
            // The off-axis coordinates sit mid-interval so the sweep reads varying
            // density rather than a node plane.
            double pt[3] = {origin[0] + 3.5 * spacing[0],
                            origin[1] + 4.5 * spacing[1],
                            origin[2] + 5.5 * spacing[2]};

            // Under the map's own cell the period is n - 1, so its final interval
            // starts at node n - 2; under n * spacing the period is n and it
            // starts at node n - 1.
            pt[axis] = origin[axis] + (n[axis] - 2u + frac) * spacing[axis];
            const double map_v = interpolate_density_periodic(
                grid, pt[0], pt[1], pt[2], map_cell[0], map_cell[1], map_cell[2]);
            map_lo = std::min(map_lo, map_v);
            map_hi = std::max(map_hi, map_v);

            pt[axis] = origin[axis] + (n[axis] - 1u + frac) * spacing[axis];
            const double full_v = interpolate_density_periodic(
                grid, pt[0], pt[1], pt[2], full_cell[0], full_cell[1], full_cell[2]);
            full_lo = std::min(full_lo, full_v);
            full_hi = std::max(full_hi, full_v);
        }

        EXPECT_GT(map_hi - map_lo, 1e-3)
            << "the final interval under the map's own cell is flat";
        EXPECT_LT(full_hi - full_lo, 1e-12)
            << "the n * spacing period no longer flattens the final interval, so this "
               "case cannot tell the two periods apart and pins nothing";
    }
}

TEST(GridOpsTest, InterpolateDensityPeriodicRejectsAGridCoveringPartOfItsCell) {
    // A grid can come off the reader with no period at all. test_map.ccp4 stores
    // 21 of its cell's 42 samples at NxSTART = -10, so it covers half the cell,
    // and its 21 A edge is 42 of the grid's 0.5 A intervals -- neither the 21
    // nodes nor the 20 intervals the rule admits. Widening the rule to admit
    // n - 1 must not have widened it to this.
    OESystem::OESkewGrid grid;
    const std::string path = std::string(MAPTITUDE_TEST_DATA_DIR) + "/test_map.ccp4";
    ASSERT_TRUE(OESystem::OEReadGrid(path, grid)) << "failed to read " << path;

    const GridParams gp = get_grid_params(grid);
    const UnitCellParams uc = get_unit_cell(grid);
    ASSERT_EQ(gp.x_dim, 21u);
    ASSERT_NEAR(gp.x_spacing, 0.5, 1e-9);
    ASSERT_NEAR(uc.a, 21.0, 1e-6);

    EXPECT_THROW(interpolate_density_periodic(grid, gp.x_origin, gp.y_origin, gp.z_origin,
                                              uc.a, uc.b, uc.c),
                 CellError);
    const std::vector<double> points = {gp.x_origin, gp.y_origin, gp.z_origin};
    EXPECT_THROW(interpolate_density_periodic_batch(grid, points, 1, uc.a, uc.b, uc.c),
                 CellError);
}

TEST(GridOpsTest, InterpolateDensityAcceptsTheNominalNodeOrigin) {
    // The node span is derived from ElementToSpatialCoord, and the error in it is
    // proportional to the magnitude of the floats the geometry is held in:
    // OpenEye holds the grid centre as a float, so a grid nominally starting at
    // 12.3 A reports its first node 1.9e-7 A away, eight orders further out than
    // the 1.5e-15 A a grid at the Cartesian origin shows from the cos(90 deg)
    // matrix term alone. Without a magnitude-aware tolerance the caller's own
    // construction coordinate falls outside the node span at most origins float
    // does not represent exactly, and float represents few origins exactly. Only
    // 0.0 was covered before.
    constexpr unsigned int N = 5u;
    const double NODE0[] = {0.0, 0.1, 3.7, 12.3, -37.45};

    for (const double node0 : NODE0) {
        SCOPED_TRACE(node0);
        OESystem::OESkewGrid grid = MakeCubicGrid(N, 1.0, node0);
        float* values = grid.GetValues();
        for (unsigned int i = 0; i < grid.GetSize(); ++i) values[i] = 42.0f;

        const GridParams gp = get_grid_params(grid);
        ASSERT_NE(gp.x_origin, node0)
            << "the derived origin now matches its nominal value exactly, so this case "
               "pins nothing";

        const double lo = node0;
        const double hi = node0 + (N - 1);
        EXPECT_DOUBLE_EQ(interpolate_density(grid, lo, lo, lo, -99.0), 42.0);
        EXPECT_DOUBLE_EQ(interpolate_density(grid, hi, hi, hi, -99.0), 42.0);
        EXPECT_DOUBLE_EQ(
            interpolate_density_periodic(grid, lo, lo, lo, N, N, N, -99.0), 42.0);
    }
}

TEST(GridOpsTest, InterpolateDensityPeriodicAcceptsAGridsOwnExtentFarFromTheOrigin) {
    // The commensurability check compares a caller-exact edge against a derived
    // n * spacing, so unlike same_grid_geometry -- which compares two derived
    // quantities, where the float noise is common-mode -- nothing cancels here.
    // The noise grows with the coordinates the grid sits at, so a fixed relative
    // tolerance starts refusing grids that tile their cell perfectly once the
    // origin is a few hundred spans from zero.
    constexpr unsigned int N = 20u;
    constexpr double SPACING = 0.9020833333333;  // 1d26's node interval
    // The last entry centres the grid on the Cartesian origin, which is the
    // hardest of these five: the endpoints are as close to zero as this span
    // allows while the cell edge -- itself a stored float -- is 2n / (n - 1)
    // times larger, which at twenty nodes is 2.1x. That ratio grows as n falls,
    // and the case that actually needs the cell edge in the scale is the small
    // one below; these cases cover the origin-distance axis instead.
    const double NODE0[] = {0.0, 100.0, 1000.0, 3000.0, -0.5 * (N - 1u) * SPACING};
    const double extent = N * SPACING;

    for (const double node0 : NODE0) {
        SCOPED_TRACE(node0);
        OESystem::OESkewGrid grid = MakeCubicGrid(N, SPACING, node0);
        float* values = grid.GetValues();
        for (unsigned int i = 0; i < grid.GetSize(); ++i) values[i] = 1.0f;

        EXPECT_NO_THROW(interpolate_density_periodic(grid, node0, node0, node0,
                                                     extent, extent, extent, -99.0));
    }

    // Discriminating power survives at the far end: one node interval too many is
    // still a different lattice, and is still refused.
    OESystem::OESkewGrid distant = MakeCubicGrid(N, SPACING, 3000.0);
    float* values = distant.GetValues();
    for (unsigned int i = 0; i < distant.GetSize(); ++i) values[i] = 1.0f;
    const double one_node_too_wide = (N + 1u) * SPACING;
    EXPECT_THROW(interpolate_density_periodic(distant, 3000.0, 3000.0, 3000.0,
                                              one_node_too_wide, extent, extent, -99.0),
                 CellError);
}

TEST(GridOpsTest, InterpolateDensityPeriodicKeepsItsAllowanceBelowANodeInterval) {
    // The count test rounds the edge-to-spacing ratio with llround, which leaves
    // |edge - p * spacing| <= spacing / 2 for the p it picks. An allowance at half
    // a node interval is therefore met by every edge that clears the count test,
    // and the commensurability comparison stops discriminating. Distance from the
    // Cartesian origin alone takes the allowance there, because AxisMagnitude
    // carries the origin: the case below would otherwise accept a cell edge that
    // misses every whole multiple of the grid's node spacing by at least four
    // tenths of an interval.
    constexpr unsigned int N = 12u;
    constexpr double SPACING = 2.0;
    constexpr double NODE0 = 3.0e6;

    OESystem::OESkewGrid grid = MakeCubicGrid(N, SPACING, NODE0);
    float* values = grid.GetValues();
    for (unsigned int i = 0; i < grid.GetSize(); ++i) values[i] = 1.0f;

    const GridParams gp = get_grid_params(grid);
    const double extent = gp.x_dim * gp.x_spacing;
    const double probe = gp.x_origin + 1.5 * gp.x_spacing;

    // Mirrors src/Grid.cpp, which keeps both file-local; see the note on
    // CELL_EXTENT_ROUNDINGS there.
    constexpr double FLOAT_HALF_ULP = 0x1p-24;
    constexpr double CELL_EXTENT_ROUNDINGS = 8.0;
    const double axis_magnitude =
        std::max(std::max(std::abs(gp.x_origin),
                          std::abs(gp.x_origin + (gp.x_dim - 1u) * gp.x_spacing)),
                 gp.x_dim * gp.x_spacing);

    // The premise, asserted rather than asserted-by-comment: this grid can only
    // show what the cap does while its uncapped allowance really does clear half
    // an interval. A later change to the rounding count that moves it back under
    // fails here instead of leaving an EXPECT_THROW that passes for the wrong
    // reason.
    ASSERT_GT(CELL_EXTENT_ROUNDINGS * FLOAT_HALF_ULP * axis_magnitude,
              0.5 * gp.x_spacing);

    // The cap narrows the allowance; it does not move what a commensurate edge
    // has to match. The grid's own extent still passes.
    EXPECT_NO_THROW(interpolate_density_periodic(grid, probe, probe, probe,
                                                 extent, extent, extent, -99.0));

    // Four tenths of an interval too wide still rounds to the same node count, so
    // the count test passes it through and the allowance is the only thing left
    // that can refuse it.
    const double off_lattice = extent + 0.4 * gp.x_spacing;
    EXPECT_THROW(interpolate_density_periodic(grid, probe, probe, probe,
                                              off_lattice, extent, extent, -99.0),
                 CellError);
    EXPECT_THROW(interpolate_density_periodic(grid, probe, probe, probe,
                                              extent, extent - 0.4 * gp.z_spacing,
                                              extent, -99.0),
                 CellError);
}

TEST(GridOpsTest, InterpolateDensityPeriodicAcceptsASmallCentredGridsOwnExtent) {
    // The allowance is scaled by the largest magnitude in the axis's geometry,
    // which includes the cell edge and not only the two endpoints. On a grid
    // centred on the Cartesian origin the endpoints are +/- (n - 1) * spacing / 2
    // while the edge is n * spacing, so the edge leads by 2n / (n - 1): 4x at two
    // nodes, 2.7x at four, and tending to 2x as n grows, which is why a small
    // centred grid is the shape to look at. Scaled by the endpoints alone, both
    // of these grids are refused a cell they tile exactly, where the twenty-node
    // case above is not.
    //
    // The ratio picks the shape but does not decide the outcome, because nothing
    // in it depends on the spacing. What also has to be large is the gap between
    // the derived extent and the nominal one, and that turns on how the geometry
    // rounds at this particular spacing: of 4000 spacings from 0.005 to 20 A,
    // only 46 put the gap past the endpoint-only allowance on the two-node
    // centred grid and 13 on the four-node one, and 5.5, 5.4 and 0.9 are all well
    // inside it. SPACING is load-bearing, not illustrative. The assertion below
    // is what says so: without it, rounding the constant to 5.5 would leave an
    // EXPECT_NO_THROW that passes under either scale.
    constexpr double SPACING = 5.45;
    const unsigned int DIMS[] = {2u, 4u};

    // Mirrors src/Grid.cpp, which keeps both file-local. Exporting them would
    // widen the public surface for a test's benefit, which is the worse trade.
    constexpr double FLOAT_HALF_ULP = 0x1p-24;
    constexpr double CELL_EXTENT_ROUNDINGS = 8.0;

    for (const unsigned int n : DIMS) {
        SCOPED_TRACE(n);
        OESystem::OESkewGrid grid = MakeCubicGrid(n, SPACING, -0.5 * (n - 1u) * SPACING);
        float* values = grid.GetValues();
        for (unsigned int i = 0; i < grid.GetSize(); ++i) values[i] = 1.0f;

        const double extent = n * SPACING;

        // The premise, in the unit the scale counts in: half-ulps of float times
        // the magnitude it is handed. Against the further endpoint -- the
        // magnitude the endpoint-only scale would use -- the derived extent is
        // 8.2 of those units off nominal at two nodes and 9.4 at four, past the
        // eight allowed, so that scale really would refuse both of these grids.
        const GridParams gp = get_grid_params(grid);
        const double far_endpoint = std::max(
            std::abs(gp.x_origin), std::abs(gp.x_origin + (gp.x_dim - 1u) * gp.x_spacing));
        const double deviation = std::abs(gp.x_dim * gp.x_spacing - extent);
        ASSERT_GT(deviation / (FLOAT_HALF_ULP * far_endpoint), CELL_EXTENT_ROUNDINGS)
            << "SPACING no longer puts the derived extent outside the endpoint-only "
               "allowance, so this grid would be accepted with the cell edge dropped "
               "from the scale and the case pins nothing";

        EXPECT_NO_THROW(
            interpolate_density_periodic(grid, 0.0, 0.0, 0.0, extent, extent, extent, -99.0));

        // One node interval too many is still a different lattice at this shape,
        // so widening the allowance to admit the grid's own extent has not cost
        // the check its purpose.
        EXPECT_THROW(interpolate_density_periodic(grid, 0.0, 0.0, 0.0,
                                                  (n + 1u) * SPACING, extent, extent, -99.0),
                     CellError);
    }
}

TEST(GridOpsTest, InterpolateDensityPeriodicReturnsTheDefaultForANonFiniteCoordinate) {
    // fmod of a non-finite fractional index is NaN, and floor(NaN) cast to
    // unsigned is undefined, so the periodic path's isfinite guard is what stands
    // between an infinite query and an arbitrary array index.
    auto grid = MakeTestGrid();
    const double nan_value = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    constexpr double DEFAULT = -99.0;

    EXPECT_DOUBLE_EQ(
        interpolate_density_periodic(grid, nan_value, 3.0, 4.0, 10.0, 10.0, 10.0, DEFAULT),
        DEFAULT);
    EXPECT_DOUBLE_EQ(
        interpolate_density_periodic(grid, 2.0, nan_value, 4.0, 10.0, 10.0, 10.0, DEFAULT),
        DEFAULT);
    EXPECT_DOUBLE_EQ(
        interpolate_density_periodic(grid, 2.0, 3.0, nan_value, 10.0, 10.0, 10.0, DEFAULT),
        DEFAULT);
    EXPECT_DOUBLE_EQ(
        interpolate_density_periodic(grid, inf, 3.0, 4.0, 10.0, 10.0, 10.0, DEFAULT), DEFAULT);
    EXPECT_DOUBLE_EQ(
        interpolate_density_periodic(grid, -inf, 3.0, 4.0, 10.0, 10.0, 10.0, DEFAULT), DEFAULT);
    EXPECT_DOUBLE_EQ(
        interpolate_density_periodic(grid, 2.0, 3.0, inf, 10.0, 10.0, 10.0, DEFAULT), DEFAULT);
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

TEST(GridOpsTest, LeavesAShiftedMoleculeAtALatticeTranslateWhenTheSizingThrows) {
    // The sizing errors are raised after the in-place shift, and the case above uses a
    // molecule the shift leaves alone, so nothing pinned where a molecule that is moved
    // ends up. The header promises the original position plus whole multiples of the
    // cell vectors -- a crystallographically equivalent place rather than a corrupted
    // one -- and that promise is what makes not rolling the shift back defensible.
    //
    // The atom sits one cell along x from the 9.4 A the case above uses, so
    // round((4.5 - 19.4) / 10) is -1 and the shift moves it by exactly one cell vector,
    // landing it where that case starts: past the node span, where a zero padding
    // cannot give the padded x axis a second node.
    constexpr double CELL = 10.0;
    OESystem::OESkewGrid grid = MakeTestGrid();  // spacing 1.0, node span [0, 9], centre 4.5
    auto mol = MakeTestMol(9.4 + CELL, 4.5, 4.5);

    OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms();
    float before[3];
    mol.GetCoords(&(*atom), before);

    EXPECT_THROW(wrap_and_pad_grid(grid, mol, CELL, CELL, CELL, 0.0), GridError);

    float after[3];
    mol.GetCoords(&(*atom), after);

    // Without a shift the assertion below would hold vacuously, on the very path the
    // case above already covers.
    ASSERT_NE(after[0], before[0]) << "no shift ran, so this case pins nothing";

    for (int i = 0; i < 3; ++i) {
        SCOPED_TRACE(i);
        const double displacement =
            static_cast<double>(after[i]) - static_cast<double>(before[i]);
        const double cells = displacement / CELL;
        EXPECT_NEAR(cells, std::round(cells), 1e-6)
            << "the molecule moved " << displacement << " A on this axis, which is not a"
            << " whole number of " << CELL << " A cell vectors";
    }
}

TEST(GridOpsTest, WrapAndPadGridRejectsACellTheGridDoesNotTile) {
    // The padded grid is filled by periodic sampling, so it inherits the periodic
    // path's precondition: the cell has to be n or n - 1 of the grid's own node
    // intervals. Twelve is neither, on a grid that samples ten.
    OESystem::OESkewGrid grid = MakeTestGrid();  // 10 nodes at spacing 1.0
    auto mol = MakeTestMol(4.5, 4.5, 4.5);
    EXPECT_THROW(wrap_and_pad_grid(grid, mol, 12.0, 10.0, 10.0, 4.75), CellError);
}

TEST(GridOpsTest, WrapAndPadGridRejectsAnIncommensurateCellBeforeMovingTheMolecule) {
    // The cell's first use is the centroid shift, so that is what it has to be
    // valid for. Checked only on the sampling path, a molecule that fits after the
    // shift took the nullptr shortcut before any check ran: the caller got their
    // coordinates translated by a vector that is not a lattice vector of the map,
    // no padded grid, and no error. Here the wrong cell would move the atom to
    // x = 25 + round((4.5 - 25)/12) * 12 = 25 + round(-1.71) * 12 = 1.0 -- inside
    // the span with room for a 0.5 A padding -- where the correct cell puts it
    // at 5.0.
    OESystem::OESkewGrid grid = MakeTestGrid();  // 10 nodes at spacing 1.0, extent 10
    auto mol = MakeTestMol(25.0, 4.5, 4.5);

    EXPECT_THROW(wrap_and_pad_grid(grid, mol, 12.0, 10.0, 10.0, 0.5), CellError);

    float coords[3];
    OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms();
    mol.GetCoords(&(*atom), coords);
    EXPECT_FLOAT_EQ(coords[0], 25.0f) << "the molecule was moved by a rejected cell";
    EXPECT_FLOAT_EQ(coords[1], 4.5f);
    EXPECT_FLOAT_EQ(coords[2], 4.5f);
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

// The cell is now validated against the grid before the molecule is inspected, so
// these three pass MakeTestGrid's own sampled extent of 10 A. An arbitrary cell
// would reach CellError first and stop testing the heavy-atom contract.
TEST(GridOpsTest, WrapAndPadThrowsWhenTheMoleculeHasNoHeavyAtoms) {
    OEChem::OEGraphMol mol;  // empty
    OESystem::OESkewGrid grid = MakeTestGrid();
    EXPECT_THROW(wrap_and_pad_grid(grid, mol, 10.0, 10.0, 10.0), StructureError);
}

TEST(GridOpsTest, WrapAndPadThrowsForAMoleculeOfOnlyDummyAtoms) {
    OEChem::OEGraphMol mol;
    OEChem::OEAtomBase* atom = mol.NewAtom(0);  // Dummy atom (Z=0)
    const float coords[3] = {4.5f, 4.5f, 4.5f};
    mol.SetCoords(atom, coords);

    OESystem::OESkewGrid grid = MakeTestGrid();
    EXPECT_THROW(wrap_and_pad_grid(grid, mol, 10.0, 10.0, 10.0), StructureError);
}

TEST(GridOpsTest, WrapAndPadThrowsForAnAllHydrogenMolecule) {
    OEChem::OEGraphMol mol;
    OEChem::OEAtomBase* atom = mol.NewAtom(1);  // Hydrogen (Z=1)
    const float coords[3] = {4.5f, 4.5f, 4.5f};
    mol.SetCoords(atom, coords);

    OESystem::OESkewGrid grid = MakeTestGrid();
    EXPECT_THROW(wrap_and_pad_grid(grid, mol, 10.0, 10.0, 10.0), StructureError);
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
        wrap_and_pad_grid(grid, mol, 10.0, 10.0, 10.0));
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
        //
        // The slack is the sizing code's own PAD_INTERVAL_COUNT_TOL against the
        // extent, not a tighter figure of this test's choosing. Snapping a
        // near-whole interval count deliberately accepts a shortfall of up to
        // that much -- a hundred-thousandth of an Angstrom on a ten-Angstrom
        // extent, far below any padding a caller would ask for -- and a test that
        // demanded more would be asserting a guarantee the code does not make.
        const double origin[3] = {gp.x_origin, gp.y_origin, gp.z_origin};
        const double spacing[3] = {gp.x_spacing, gp.y_spacing, gp.z_spacing};
        const unsigned int dim[3] = {gp.x_dim, gp.y_dim, gp.z_dim};
        for (int i = 0; i < 3; ++i) {
            const double lo = std::min(c.atom_lo[i], c.atom_hi[i]) - c.padding;
            const double hi = std::max(c.atom_lo[i], c.atom_hi[i]) + c.padding;
            const double slack = PAD_INTERVAL_COUNT_TOL * (hi - lo);
            EXPECT_LE(origin[i], lo + slack)
                << "axis " << i << " node span starts inside the required extent";
            EXPECT_GE(origin[i] + (dim[i] - 1) * spacing[i], hi - slack)
                << "axis " << i << " node span ends inside the required extent";
        }
    }
}
