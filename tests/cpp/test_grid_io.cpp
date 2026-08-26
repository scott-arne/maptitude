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
    OESystem::OEScalarGrid grid;
    ASSERT_TRUE(OESystem::OEReadGrid(TestMapPath(), grid)) << "failed to read " << TestMapPath();

    EXPECT_GT(grid.GetSize(), 0u);
    EXPECT_GT(grid.GetXDim(), 1u);
    EXPECT_GT(grid.GetYDim(), 1u);
    EXPECT_GT(grid.GetZDim(), 1u);
    EXPECT_GT(grid.GetSpacing(), 0.0);
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
    // A carbon on the Gaussian's peak scores ~0.94. Asserting strict positivity
    // is what distinguishes this from a test that an always-zero regression in
    // the disk-read scoring path would still pass.
    EXPECT_GT(result.overall, 0.0);
    EXPECT_FALSE(result.by_atom.empty());
}
