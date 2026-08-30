/// Tier 3: exercise the real CCP4 file path.
///
/// The synthetic fixtures build grids in memory and never touch a reader. This
/// test is the only coverage of OEReadGrid and of scoring against a map that
/// came off disk.
#include <gtest/gtest.h>

#include <cmath>
#include <string>

#include <oegrid.h>
#include <oesystem.h>

#include "maptitude/Grid.h"
#include "maptitude/Metric.h"

#include "fixtures.h"

using namespace Maptitude;
using namespace MaptitudeTest;

namespace {
std::string TestMapPath() {
    return std::string(MAPTITUDE_TEST_DATA_DIR) + "/test_map.ccp4";
}
}  // namespace

TEST(GridIoTest, ReadsCommittedCcp4Map) {
    OESystem::OESkewGrid grid;
    ASSERT_TRUE(OESystem::OEReadGrid(TestMapPath(), grid)) << "failed to read " << TestMapPath();

    // The fixture is committed, so its geometry is an exact contract rather
    // than a range to bound.
    EXPECT_EQ(grid.GetSize(), 9261u);
    EXPECT_EQ(grid.GetXDim(), 21u);
    EXPECT_EQ(grid.GetYDim(), 21u);
    EXPECT_EQ(grid.GetZDim(), 21u);
    // The fixture is isotropic, so one spacing held for it before; asserting all
    // three is what would now catch a reader that lost an axis's scale.
    const GridParams gp = get_grid_params(grid);
    EXPECT_DOUBLE_EQ(gp.x_spacing, 0.5);
    EXPECT_DOUBLE_EQ(gp.y_spacing, 0.5);
    EXPECT_DOUBLE_EQ(gp.z_spacing, 0.5);

    // Element index to coordinate, probed one step along each stride. This is
    // the axis-order guard: x varies fastest at stride 1, y at 21, z at 441, so
    // a reader that transposed the axes moves these three coordinates even
    // though the cubic isotropic payload itself is transposition-invariant.
    float x = 0.0f, y = 0.0f, z = 0.0f;
    ASSERT_TRUE(grid.ElementToSpatialCoord(0u, x, y, z));
    EXPECT_FLOAT_EQ(x, -5.0f);
    EXPECT_FLOAT_EQ(y, -5.0f);
    EXPECT_FLOAT_EQ(z, -5.0f);
    ASSERT_TRUE(grid.ElementToSpatialCoord(1u, x, y, z));
    EXPECT_FLOAT_EQ(x, -4.5f);
    EXPECT_FLOAT_EQ(y, -5.0f);
    EXPECT_FLOAT_EQ(z, -5.0f);
    ASSERT_TRUE(grid.ElementToSpatialCoord(21u, x, y, z));
    EXPECT_FLOAT_EQ(x, -5.0f);
    EXPECT_FLOAT_EQ(y, -4.5f);
    EXPECT_FLOAT_EQ(z, -5.0f);
    ASSERT_TRUE(grid.ElementToSpatialCoord(441u, x, y, z));
    EXPECT_FLOAT_EQ(x, -5.0f);
    EXPECT_FLOAT_EQ(y, -5.0f);
    EXPECT_FLOAT_EQ(z, -4.5f);

    // Payload: the Gaussian peaks at the centre and decays to a tiny but
    // strictly positive corner. Counting the nonzero voxels is what rules out a
    // read that populated only part of the map.
    const float* values = grid.GetValues();
    ASSERT_NE(values, nullptr);
    EXPECT_FLOAT_EQ(values[4630], 1.0f);
    EXPECT_GT(values[0], 0.0f);
    EXPECT_LT(values[0], 1.0e-9f);
    unsigned int positive = 0u;
    for (unsigned int i = 0u; i < grid.GetSize(); ++i) {
        if (values[i] > 0.0f) {
            ++positive;
        }
    }
    EXPECT_EQ(positive, grid.GetSize());
}

TEST(GridIoTest, ScoresAgainstAMapReadFromDisk) {
    OESystem::OEScalarGrid grid;
    ASSERT_TRUE(OESystem::OEReadGrid(TestMapPath(), grid));

    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    DensityScoreResult result = qscore(mol, grid, 2.0);

    // A smoke test, not a pin: assert only that a real map produces a finite,
    // in-range score through the full scoring path.
    EXPECT_TRUE(std::isfinite(result.overall));
    EXPECT_GE(result.overall, -1.0);
    EXPECT_LE(result.overall, 1.0);
    // Measured 0.94143 for a carbon on the Gaussian peak. The 0.8 floor is a
    // wide-margin guard against a read that lands on the wrong region and still
    // returns a small positive score -- not a pinned value.
    EXPECT_GT(result.overall, 0.8);
    ASSERT_EQ(result.by_atom.size(), 1u);
    EXPECT_DOUBLE_EQ(result.by_atom.at(0), result.overall);
}
