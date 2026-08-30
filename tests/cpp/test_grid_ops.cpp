#include <gtest/gtest.h>
#include "maptitude/Grid.h"
#include "maptitude/GridOps.h"
#include "maptitude/Error.h"

#include <oechem.h>
#include <oegrid.h>

#include <chrono>
#include <cmath>
#include <iostream>
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

    OESystem::OEScalarGrid* result = wrap_and_pad_grid(
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

    OESystem::OEScalarGrid* result = wrap_and_pad_grid(
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

    OESystem::OEScalarGrid* result = wrap_and_pad_grid(
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
    EXPECT_NEAR((*result)[0], static_cast<float>(val), 0.1);

    delete result;
}

static OESystem::OEScalarGrid MakeShiftedGrid(double shift) {
    double minmax[6] = {shift, shift, shift, 9.0 + shift, 9.0 + shift, 9.0 + shift};
    OESystem::OEScalarGrid grid(minmax, 1.0);
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        grid[i] = 1.0f;
    }
    return grid;
}

TEST(GridOpsTest, CombineRejectsGridsWithDifferentOrigins) {
    // Same dims, same spacing, different origin. Element-wise combination would
    // mix densities from different points in space.
    OESystem::OEScalarGrid a = MakeShiftedGrid(0.0);
    OESystem::OEScalarGrid b = MakeShiftedGrid(5.0);
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
    OESystem::OEScalarGrid a(19, 19, 19, 4.5, 4.5, 4.5, 0.5);
    OESystem::OEScalarGrid b(19, 19, 19, 4.5, 4.5, 4.5, static_cast<float>(0.5 + 1e-7));
    // Spacing is the only difference among the quantities OEGridSameGeometry compares --
    // dimensions, midpoints, and spacing -- which is what makes the rejection below the
    // spacing's. It is not the only difference between the two grids: OEScalarGrid stores
    // a midpoint and derives the origin from it, so the spacing delta moves the origin as
    // well, by 9.54e-7 A here. An earlier version of this comment claimed spacing was the
    // only difference outright.
    ASSERT_EQ(a.GetXDim(), b.GetXDim()) << "the dimensions must match or this pins nothing";
    ASSERT_EQ(a.GetXMid(), b.GetXMid()) << "the midpoints must match or this pins nothing";
    ASSERT_NEAR(std::fabs(a.GetXMin() - b.GetXMin()), 9.5367431640625e-7, 1e-13)
        << "the origins were expected to move with the spacing";
    ASSERT_NE(a.GetSpacing(), b.GetSpacing()) << "the two spacings collapsed to one float";
    ASSERT_LT(std::fabs(a.GetSpacing() - b.GetSpacing()), 1e-6)
        << "the delta must sit inside the old tolerance or this pins nothing";
    EXPECT_THROW(combine_maps(a, b, MapOp::ADD), GridError);
}

TEST(GridOpsTest, CombineStillAcceptsIdenticalGeometry) {
    OESystem::OEScalarGrid a = MakeTestGrid();
    OESystem::OEScalarGrid b = MakeTestGrid();
    EXPECT_NO_THROW({
        std::unique_ptr<OESystem::OEScalarGrid> result(combine_maps(a, b, MapOp::ADD));
        ASSERT_NE(result, nullptr);
    });
}

TEST(GridOpsTest, DiffToCalcRejectsMismatchedGeometry) {
    OESystem::OEScalarGrid obs = MakeShiftedGrid(0.0);
    OESystem::OEScalarGrid diff = MakeShiftedGrid(5.0);
    EXPECT_THROW(diff_to_calc(obs, diff), GridError);
}

TEST(GridOpsTest, WrapAndPadThrowsWhenTheMoleculeHasNoHeavyAtoms) {
    OEChem::OEGraphMol mol;  // empty
    OESystem::OEScalarGrid grid = MakeTestGrid();
    EXPECT_THROW(wrap_and_pad_grid(grid, mol, 20.0, 25.0, 30.0), StructureError);
}

TEST(GridOpsTest, WrapAndPadThrowsForAMoleculeOfOnlyDummyAtoms) {
    OEChem::OEGraphMol mol;
    OEChem::OEAtomBase* atom = mol.NewAtom(0);  // Dummy atom (Z=0)
    const float coords[3] = {4.5f, 4.5f, 4.5f};
    mol.SetCoords(atom, coords);

    OESystem::OEScalarGrid grid = MakeTestGrid();
    EXPECT_THROW(wrap_and_pad_grid(grid, mol, 20.0, 25.0, 30.0), StructureError);
}

TEST(GridOpsTest, WrapAndPadThrowsForAnAllHydrogenMolecule) {
    OEChem::OEGraphMol mol;
    OEChem::OEAtomBase* atom = mol.NewAtom(1);  // Hydrogen (Z=1)
    const float coords[3] = {4.5f, 4.5f, 4.5f};
    mol.SetCoords(atom, coords);

    OESystem::OEScalarGrid grid = MakeTestGrid();
    EXPECT_THROW(wrap_and_pad_grid(grid, mol, 20.0, 25.0, 30.0), StructureError);
}

TEST(GridOpsTest, WrapAndPadReturnsNullptrOnlyWhenNoPaddingIsNeeded) {
    // A molecule already well inside the grid needs no padding: nullptr means
    // "unchanged", and nothing else.
    OEChem::OEGraphMol mol;
    OEChem::OEAtomBase* atom = mol.NewAtom(6);
    const double coords[3] = {4.5, 4.5, 4.5};
    mol.SetCoords(atom, coords);

    OESystem::OEScalarGrid grid = MakeTestGrid();
    std::unique_ptr<OESystem::OEScalarGrid> result(
        wrap_and_pad_grid(grid, mol, 20.0, 25.0, 30.0));
    EXPECT_EQ(result, nullptr);
}

// Disabled by default: it reads a 222 KB asset and runs millions of samples.
// Run with:
//   ./build-debug/tests/cpp/maptitude_tests \
//     --gtest_also_run_disabled_tests --gtest_filter=Interpolation*Timing*
TEST(InterpolationTiming, DISABLED_OpenEyeVersusMaptitudeOn1d26) {
    const std::string path = std::string(MAPTITUDE_TEST_ASSET_DIR) + "/1d26_2fofc.ccp4";

    // The OpenEye arm of the timing comparison needs the scalar overload of
    // OEFloatGridLinearInterpolate, which is what this test measures against.
    OESystem::OEScalarGrid scalar;  // OE-SCALARGRID-OK: the OpenEye arm of the comparison
    ASSERT_TRUE(OESystem::OEReadGrid(path, scalar)) << "failed to read " << path;
    const OESystem::OESkewGrid skew(scalar);

    const GridParams gp = get_grid_params(skew);
    const float* values = skew.GetValues();
    ASSERT_NE(values, nullptr);

    // A deterministic sweep of the interior, avoiding the faces so both arms
    // sample the same domain.
    std::vector<double> pts;
    for (unsigned int iz = 1; iz + 1 < gp.z_dim; ++iz) {
        for (unsigned int iy = 1; iy + 1 < gp.y_dim; ++iy) {
            for (unsigned int ix = 1; ix + 1 < gp.x_dim; ++ix) {
                pts.push_back(gp.x_origin + (ix + 0.37) * gp.x_spacing);
                pts.push_back(gp.y_origin + (iy + 0.61) * gp.y_spacing);
                pts.push_back(gp.z_origin + (iz + 0.19) * gp.z_spacing);
            }
        }
    }
    const size_t n = pts.size() / 3;
    constexpr int REPEATS = 20;

    double sink = 0.0;
    auto t0 = std::chrono::steady_clock::now();
    for (int r = 0; r < REPEATS; ++r) {
        for (size_t i = 0; i < n; ++i) {
            sink += OESystem::OEFloatGridLinearInterpolate(
                scalar, static_cast<float>(pts[i * 3]),
                static_cast<float>(pts[i * 3 + 1]),
                static_cast<float>(pts[i * 3 + 2]), 0.0f);
        }
    }
    auto t1 = std::chrono::steady_clock::now();
    for (int r = 0; r < REPEATS; ++r) {
        for (size_t i = 0; i < n; ++i) {
            sink += interpolate_density_at(gp, values, pts[i * 3], pts[i * 3 + 1],
                                           pts[i * 3 + 2], 0.0);
        }
    }
    auto t2 = std::chrono::steady_clock::now();

    const double oe_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    const double mt_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
    std::cout << "points=" << n << " repeats=" << REPEATS
              << " OEFloatGridLinearInterpolate=" << oe_ms << " ms"
              << " interpolate_density_at=" << mt_ms << " ms"
              << " ratio=" << (oe_ms / mt_ms) << "\n";
    EXPECT_GT(sink, 0.0) << "sink keeps both loops from being optimized away";
}

// ---- Interpolation equivalence guard (spec §3.4) ----
//
// Fixture: a 4x4x4 isotropic grid with unit spacing whose node origin is the
// Cartesian origin, so a coordinate is its own fractional index.
//
// The value field is deliberately NOT multilinear. A trilinear blend
// reproduces a multilinear field exactly, so a guard built on one would agree
// with any index or weight error that happens to be self-consistent.
namespace {

float GuardFieldValue(const double x, const double y, const double z) {
    return static_cast<float>(x * x + 2.0 * y * y * y - 0.5 * z * z +
                              3.0 * x * y - y * z);
}

// The guard's whole purpose is to compare against OpenEye's scalar-grid
// interpolation, so this arm cannot move to the skew carrier.
OESystem::OEScalarGrid MakeGuardScalarGrid() {  // OE-SCALARGRID-OK: the guard's OpenEye arm
    double minmax[6] = {0.0, 0.0, 0.0, 3.0, 3.0, 3.0};
    OESystem::OEScalarGrid grid(minmax, 1.0);  // OE-SCALARGRID-OK: the guard's OpenEye arm
    for (unsigned int iz = 0; iz < 4u; ++iz) {
        for (unsigned int iy = 0; iy < 4u; ++iy) {
            for (unsigned int ix = 0; ix < 4u; ++ix) {
                grid[iz * 16u + iy * 4u + ix] = GuardFieldValue(ix, iy, iz);
            }
        }
    }
    return grid;
}

// Fractional index, then OEFloatGridLinearInterpolate's return at it.
// Captured in Task 2; the guard's reference after Task 3 removes the call.
struct GuardSample { double f[3]; double expected; };
static constexpr GuardSample GUARD_SAMPLES[] = {
    // Filled in by Step 9 below. Every entry is strictly interior:
    // 0 < f_i < 3 on all three axes, so the domain split does not apply.
    {{0.13, 0.13, 0.13}, 0.35879999399185181},
    {{0.13, 0.13, 0.5}, 0.12569998204708099},
    {{0.13, 0.13, 1}, -0.18930001556873322},
    {{0.13, 0.13, 1.47}, -0.95540004968643188},
    {{0.13, 0.13, 2}, -1.8193000555038452},
    {{0.13, 0.13, 2.8599999999999999}, -4.0810995101928711},
    {{0.13, 0.5, 0.13}, 1.1949999332427979},
    {{0.13, 0.5, 0.5}, 0.82499998807907104},
    {{0.13, 0.5, 1}, 0.32499998807907104},
    {{0.13, 0.5, 1.47}, -0.61500006914138794},
    {{0.13, 0.5, 2}, -1.6749999523162842},
    {{0.13, 0.5, 2.8599999999999999}, -4.2549996376037598},
    {{0.13, 1, 0.13}, 2.3250000476837158},
    {{0.13, 1, 0.5}, 1.7699999809265137},
    {{0.13, 1, 1}, 1.0199999809265137},
    {{0.13, 1, 1.47}, -0.15500009059906006},
    {{0.13, 1, 2}, -1.4800000190734863},
    {{0.13, 1, 2.8599999999999999}, -4.4899997711181641},
    {{0.13, 1.47, 0.13}, 9.0272006988525391},
    {{0.13, 1.47, 0.5}, 8.2983007431030273},
    {{0.13, 1.47, 1}, 7.3133001327514648},
    {{0.13, 1.47, 1.47}, 5.9174003601074219},
    {{0.13, 1.47, 2}, 4.3433003425598145},
    {{0.13, 1.47, 2.8599999999999999}, 0.92910069227218628},
    {{0.13, 2, 0.13}, 16.584999084472656},
    {{0.13, 2, 0.5}, 15.659999847412109},
    {{0.13, 2, 1}, 14.409999847412109},
    {{0.13, 2, 1.47}, 12.764999389648438},
    {{0.13, 2, 2}, 10.909999847412109},
    {{0.13, 2, 2.8599999999999999}, 7.0400004386901855},
    {{0.13, 2.8599999999999999, 0.13}, 49.488594055175781},
    {{0.13, 2.8599999999999999, 0.5}, 48.245395660400391},
    {{0.13, 2.8599999999999999, 1}, 46.565395355224609},
    {{0.13, 2.8599999999999999, 1.47}, 44.516197204589844},
    {{0.13, 2.8599999999999999, 2}, 42.205394744873047},
    {{0.13, 2.8599999999999999, 2.8599999999999999}, 37.595798492431641},
    {{0.5, 0.13, 0.13}, 0.87309998273849487},
    {{0.5, 0.13, 0.5}, 0.63999998569488525},
    {{0.5, 0.13, 1}, 0.32499998807907104},
    {{0.5, 0.13, 1.47}, -0.44110006093978882},
    {{0.5, 0.13, 2}, -1.3050000667572021},
    {{0.5, 0.13, 2.8599999999999999}, -3.5667996406555176},
    {{0.5, 0.5, 0.13}, 2.119999885559082},
    {{0.5, 0.5, 0.5}, 1.75},
    {{0.5, 0.5, 1}, 1.25},
    {{0.5, 0.5, 1.47}, 0.30999994277954102},
    {{0.5, 0.5, 2}, -0.75},
    {{0.5, 0.5, 2.8599999999999999}, -3.3299996852874756},
    {{0.5, 1, 0.13}, 3.8050000667572021},
    {{0.5, 1, 0.5}, 3.25},
    {{0.5, 1, 1}, 2.5},
    {{0.5, 1, 1.47}, 1.3249999284744263},
    {{0.5, 1, 2}, 0},
    {{0.5, 1, 2.8599999999999999}, -3.0099997520446777},
    {{0.5, 1.47, 0.13}, 11.028900146484375},
    {{0.5, 1.47, 0.5}, 10.300000190734863},
    {{0.5, 1.47, 1}, 9.3150005340576172},
    {{0.5, 1.47, 1.47}, 7.919100284576416},
    {{0.5, 1.47, 2}, 6.3450002670288086},
    {{0.5, 1.47, 2.8599999999999999}, 2.9308006763458252},
    {{0.5, 2, 0.13}, 19.174999237060547},
    {{0.5, 2, 0.5}, 18.25},
    {{0.5, 2, 1}, 17},
    {{0.5, 2, 1.47}, 15.354999542236328},
    {{0.5, 2, 2}, 13.5},
    {{0.5, 2, 2.8599999999999999}, 9.630000114440918},
    {{0.5, 2.8599999999999999, 0.13}, 53.033195495605469},
    {{0.5, 2.8599999999999999, 0.5}, 51.789997100830078},
    {{0.5, 2.8599999999999999, 1}, 50.109996795654297},
    {{0.5, 2.8599999999999999, 1.47}, 48.060794830322266},
    {{0.5, 2.8599999999999999, 2}, 45.749996185302734},
    {{0.5, 2.8599999999999999, 2.8599999999999999}, 41.140396118164062},
    {{1, 0.13, 0.13}, 1.5680999755859375},
    {{1, 0.13, 0.5}, 1.3350000381469727},
    {{1, 0.13, 1}, 1.0199999809265137},
    {{1, 0.13, 1.47}, 0.25389993190765381},
    {{1, 0.13, 2}, -0.61000001430511475},
    {{1, 0.13, 2.8599999999999999}, -2.8717997074127197},
    {{1, 0.5, 0.13}, 3.369999885559082},
    {{1, 0.5, 0.5}, 3},
    {{1, 0.5, 1}, 2.5},
    {{1, 0.5, 1.47}, 1.559999942779541},
    {{1, 0.5, 2}, 0.5},
    {{1, 0.5, 2.8599999999999999}, -2.0799996852874756},
    {{1, 1, 0.13}, 5.804999828338623},
    {{1, 1, 0.5}, 5.25},
    {{1, 1, 1}, 4.5},
    {{1, 1, 1.47}, 3.3249998092651367},
    {{1, 1, 2}, 2},
    {{1, 1, 2.8599999999999999}, -1.0099996328353882},
    {{1, 1.47, 0.13}, 13.73390007019043},
    {{1, 1.47, 0.5}, 13.005000114440918},
    {{1, 1.47, 1}, 12.020000457763672},
    {{1, 1.47, 1.47}, 10.624100685119629},
    {{1, 1.47, 2}, 9.0500001907348633},
    {{1, 1.47, 2.8599999999999999}, 5.635800838470459},
    {{1, 2, 0.13}, 22.674999237060547},
    {{1, 2, 0.5}, 21.75},
    {{1, 2, 1}, 20.5},
    {{1, 2, 1.47}, 18.854999542236328},
    {{1, 2, 2}, 17},
    {{1, 2, 2.8599999999999999}, 13.130000114440918},
    {{1, 2.8599999999999999, 0.13}, 57.823196411132812},
    {{1, 2.8599999999999999, 0.5}, 56.579994201660156},
    {{1, 2.8599999999999999, 1}, 54.899993896484375},
    {{1, 2.8599999999999999, 1.47}, 52.850795745849609},
    {{1, 2.8599999999999999, 2}, 50.539997100830078},
    {{1, 2.8599999999999999, 2.8599999999999999}, 45.930397033691406},
    {{1.47, 0.13, 0.13}, 3.1614000797271729},
    {{1.47, 0.13, 0.5}, 2.928300142288208},
    {{1.47, 0.13, 1}, 2.613300085067749},
    {{1.47, 0.13, 1.47}, 1.8472000360488892},
    {{1.47, 0.13, 2}, 0.98330008983612061},
    {{1.47, 0.13, 2.8599999999999999}, -1.2784996032714844},
    {{1.47, 0.5, 0.13}, 5.4850001335144043},
    {{1.47, 0.5, 0.5}, 5.1150002479553223},
    {{1.47, 0.5, 1}, 4.6150002479553223},
    {{1.47, 0.5, 1.47}, 3.6750001907348633},
    {{1.47, 0.5, 2}, 2.6150002479553223},
    {{1.47, 0.5, 2.8599999999999999}, 0.035000443458557129},
    {{1.47, 1, 0.13}, 8.625},
    {{1.47, 1, 0.5}, 8.0699996948242188},
    {{1.47, 1, 1}, 7.320000171661377},
    {{1.47, 1, 1.47}, 6.1449999809265137},
    {{1.47, 1, 2}, 4.820000171661377},
    {{1.47, 1, 2.8599999999999999}, 1.8100005388259888},
    {{1.47, 1.47, 0.13}, 17.21660041809082},
    {{1.47, 1.47, 0.5}, 16.487701416015625},
    {{1.47, 1.47, 1}, 15.502700805664062},
    {{1.47, 1.47, 1.47}, 14.10680103302002},
    {{1.47, 1.47, 2}, 12.532700538635254},
    {{1.47, 1.47, 2.8599999999999999}, 9.1185007095336914},
    {{1.47, 2, 0.13}, 26.905000686645508},
    {{1.47, 2, 0.5}, 25.979999542236328},
    {{1.47, 2, 1}, 24.729999542236328},
    {{1.47, 2, 1.47}, 23.085000991821289},
    {{1.47, 2, 2}, 21.229999542236328},
    {{1.47, 2, 2.8599999999999999}, 17.360000610351562},
    {{1.47, 2.8599999999999999, 0.13}, 63.265796661376953},
    {{1.47, 2.8599999999999999, 0.5}, 62.022594451904297},
    {{1.47, 2.8599999999999999, 1}, 60.342594146728516},
    {{1.47, 2.8599999999999999, 1.47}, 58.29339599609375},
    {{1.47, 2.8599999999999999, 2}, 55.982597351074219},
    {{1.47, 2.8599999999999999, 2.8599999999999999}, 51.372997283935547},
    {{2, 0.13, 0.13}, 4.9580998420715332},
    {{2, 0.13, 0.5}, 4.7249999046325684},
    {{2, 0.13, 1}, 4.4099998474121094},
    {{2, 0.13, 1.47}, 3.6438999176025391},
    {{2, 0.13, 2}, 2.7799999713897705},
    {{2, 0.13, 2.8599999999999999}, 0.51820027828216553},
    {{2, 0.5, 0.13}, 7.869999885559082},
    {{2, 0.5, 0.5}, 7.5},
    {{2, 0.5, 1}, 7},
    {{2, 0.5, 1.47}, 6.059999942779541},
    {{2, 0.5, 2}, 5},
    {{2, 0.5, 2.8599999999999999}, 2.4200003147125244},
    {{2, 1, 0.13}, 11.805000305175781},
    {{2, 1, 0.5}, 11.25},
    {{2, 1, 1}, 10.5},
    {{2, 1, 1.47}, 9.3249998092651367},
    {{2, 1, 2}, 8},
    {{2, 1, 2.8599999999999999}, 4.9900002479553223},
    {{2, 1.47, 0.13}, 21.143899917602539},
    {{2, 1.47, 0.5}, 20.415000915527344},
    {{2, 1.47, 1}, 19.430000305175781},
    {{2, 1.47, 1.47}, 18.034099578857422},
    {{2, 1.47, 2}, 16.460000991821289},
    {{2, 1.47, 2.8599999999999999}, 13.045801162719727},
    {{2, 2, 0.13}, 31.674999237060547},
    {{2, 2, 0.5}, 30.75},
    {{2, 2, 1}, 29.5},
    {{2, 2, 1.47}, 27.854999542236328},
    {{2, 2, 2}, 26},
    {{2, 2, 2.8599999999999999}, 22.130001068115234},
    {{2, 2.8599999999999999, 0.13}, 69.4031982421875},
    {{2, 2.8599999999999999, 0.5}, 68.159996032714844},
    {{2, 2.8599999999999999, 1}, 66.479995727539062},
    {{2, 2.8599999999999999, 1.47}, 64.430793762207031},
    {{2, 2.8599999999999999, 2}, 62.1199951171875},
    {{2, 2.8599999999999999, 2.8599999999999999}, 57.510395050048828},
    {{2.8599999999999999, 0.13, 0.13}, 9.5934991836547852},
    {{2.8599999999999999, 0.13, 0.5}, 9.3603992462158203},
    {{2.8599999999999999, 0.13, 1}, 9.0453996658325195},
    {{2.8599999999999999, 0.13, 1.47}, 8.2792997360229492},
    {{2.8599999999999999, 0.13, 2}, 7.4153995513916016},
    {{2.8599999999999999, 0.13, 2.8599999999999999}, 5.153599739074707},
    {{2.8599999999999999, 0.5, 0.13}, 13.459999084472656},
    {{2.8599999999999999, 0.5, 0.5}, 13.089999198913574},
    {{2.8599999999999999, 0.5, 1}, 12.589999198913574},
    {{2.8599999999999999, 0.5, 1.47}, 11.649999618530273},
    {{2.8599999999999999, 0.5, 2}, 10.589999198913574},
    {{2.8599999999999999, 0.5, 2.8599999999999999}, 8.0099992752075195},
    {{2.8599999999999999, 1, 0.13}, 18.684999465942383},
    {{2.8599999999999999, 1, 0.5}, 18.129999160766602},
    {{2.8599999999999999, 1, 1}, 17.379999160766602},
    {{2.8599999999999999, 1, 1.47}, 16.204999923706055},
    {{2.8599999999999999, 1, 2}, 14.879999160766602},
    {{2.8599999999999999, 1, 2.8599999999999999}, 11.869999885559082},
    {{2.8599999999999999, 1.47, 0.13}, 29.236499786376953},
    {{2.8599999999999999, 1.47, 0.5}, 28.507598876953125},
    {{2.8599999999999999, 1.47, 1}, 27.522600173950195},
    {{2.8599999999999999, 1.47, 1.47}, 26.126699447631836},
    {{2.8599999999999999, 1.47, 2}, 24.55259895324707},
    {{2.8599999999999999, 1.47, 2.8599999999999999}, 21.138399124145508},
    {{2.8599999999999999, 2, 0.13}, 41.134998321533203},
    {{2.8599999999999999, 2, 0.5}, 40.209999084472656},
    {{2.8599999999999999, 2, 1}, 38.959999084472656},
    {{2.8599999999999999, 2, 1.47}, 37.314998626708984},
    {{2.8599999999999999, 2, 2}, 35.459999084472656},
    {{2.8599999999999999, 2, 2.8599999999999999}, 31.590000152587891},
    {{2.8599999999999999, 2.8599999999999999, 0.13}, 81.081993103027344},
    {{2.8599999999999999, 2.8599999999999999, 0.5}, 79.838790893554688},
    {{2.8599999999999999, 2.8599999999999999, 1}, 78.158790588378906},
    {{2.8599999999999999, 2.8599999999999999, 1.47}, 76.109596252441406},
    {{2.8599999999999999, 2.8599999999999999, 2}, 73.798797607421875},
    {{2.8599999999999999, 2.8599999999999999, 2.8599999999999999}, 69.189193725585938},
};

}  // namespace

TEST(InterpolationEquivalence, MatchesOpenEyeOnTheSharedInteriorDomain) {
    const OESystem::OEScalarGrid scalar = MakeGuardScalarGrid();  // OE-SCALARGRID-OK: the guard's OpenEye arm
    const OESystem::OESkewGrid skew(scalar);
    const GridParams gp = get_grid_params(skew);
    const float* values = skew.GetValues();
    ASSERT_NE(values, nullptr);

    for (const GuardSample& s : GUARD_SAMPLES) {
        const float oe = OESystem::OEFloatGridLinearInterpolate(
            scalar, static_cast<float>(s.f[0]), static_cast<float>(s.f[1]),
            static_cast<float>(s.f[2]), 0.0f);
        EXPECT_NEAR(oe, s.expected, 1e-5)
            << "pinned reference drifted at (" << s.f[0] << ", " << s.f[1]
            << ", " << s.f[2] << ")";
        EXPECT_NEAR(interpolate_density_at(gp, values, s.f[0], s.f[1], s.f[2], 0.0),
                    s.expected, 1e-5)
            << "maptitude disagrees at (" << s.f[0] << ", " << s.f[1]
            << ", " << s.f[2] << ")";
    }
}

// ---- Far-face divergence (spec §3.4) ----
//
// OEFloatGridLinearInterpolate treats the far face as outside the grid and
// returns the caller's default. maptitude closes the interval and blends from
// the nodes that are there. Each test pins both halves; Task 3 keeps only the
// maptitude half, because by then there is no OpenEye call left to make.
//
// -99.0f is the sentinel default: no value in this field is near it, so a
// "returned the default" result cannot be confused with a blend.

TEST(InterpolationFarFace, DivergesFromOpenEyeAtTheFarCorner) {
    const OESystem::OEScalarGrid scalar = MakeGuardScalarGrid();  // OE-SCALARGRID-OK: the divergence test's OpenEye arm
    const OESystem::OESkewGrid skew(scalar);
    const GridParams gp = get_grid_params(skew);
    const float* values = skew.GetValues();
    ASSERT_NE(values, nullptr);

    EXPECT_FLOAT_EQ(
        OESystem::OEFloatGridLinearInterpolate(scalar, 3.0f, 3.0f, 3.0f, -99.0f),
        -99.0f);
    // All three fractional indices at n - 1: no blend, the corner node itself.
    EXPECT_NEAR(interpolate_density_at(gp, values, 3.0, 3.0, 3.0, -99.0),
                GuardFieldValue(3.0, 3.0, 3.0), 1e-4);
}

TEST(InterpolationFarFace, DivergesFromOpenEyeOnAFarFace) {
    const OESystem::OEScalarGrid scalar = MakeGuardScalarGrid();  // OE-SCALARGRID-OK: the divergence test's OpenEye arm
    const OESystem::OESkewGrid skew(scalar);
    const GridParams gp = get_grid_params(skew);
    const float* values = skew.GetValues();
    ASSERT_NE(values, nullptr);

    EXPECT_FLOAT_EQ(
        OESystem::OEFloatGridLinearInterpolate(scalar, 3.0f, 1.5f, 1.5f, -99.0f),
        -99.0f);
    // One index pinned at n - 1, two interior: a bilinear blend of four nodes.
    const double expected =
        0.25 * (GuardFieldValue(3.0, 1.0, 1.0) + GuardFieldValue(3.0, 2.0, 1.0) +
                GuardFieldValue(3.0, 1.0, 2.0) + GuardFieldValue(3.0, 2.0, 2.0));
    EXPECT_NEAR(interpolate_density_at(gp, values, 3.0, 1.5, 1.5, -99.0),
                expected, 1e-4);
}

TEST(InterpolationFarFace, DivergesFromOpenEyeOnAFarEdge) {
    const OESystem::OEScalarGrid scalar = MakeGuardScalarGrid();  // OE-SCALARGRID-OK: the divergence test's OpenEye arm
    const OESystem::OESkewGrid skew(scalar);
    const GridParams gp = get_grid_params(skew);
    const float* values = skew.GetValues();
    ASSERT_NE(values, nullptr);

    EXPECT_FLOAT_EQ(
        OESystem::OEFloatGridLinearInterpolate(scalar, 3.0f, 3.0f, 1.5f, -99.0f),
        -99.0f);
    // Two indices pinned at n - 1, one interior: a linear blend of two nodes.
    const double expected =
        0.5 * (GuardFieldValue(3.0, 3.0, 1.0) + GuardFieldValue(3.0, 3.0, 2.0));
    EXPECT_NEAR(interpolate_density_at(gp, values, 3.0, 3.0, 1.5, -99.0),
                expected, 1e-4);
}
