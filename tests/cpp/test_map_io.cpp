/**
 * @file test_map_io.cpp
 * @brief Read and write pins for the CCP4/MRC map I/O path.
 *
 * The values pinned here come from the measurement table in section 2.1 of the
 * map I/O design, taken against OpenEye 2026.1.0. They are pins in the Tier 2
 * sense: they record what the toolkit and this module actually do on the
 * committed fixtures, so a change in either shows up as a failure here rather
 * than as a silently relocated map.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

#include <oegrid.h>
#include <oesystem.h>

#include "maptitude/Error.h"
#include "maptitude/Grid.h"
#include "maptitude/MapIO.h"
#include "maptitude/SymOp.h"

#include "fixtures.h"

using namespace Maptitude;
using namespace MaptitudeTest;

namespace {

std::string AssetPath(const std::string& name) {
    return std::string(MAPTITUDE_TEST_ASSET_DIR) + "/" + name;
}

std::string DataPath(const std::string& name) {
    return std::string(MAPTITUDE_TEST_DATA_DIR) + "/" + name;
}

/// One row of the section 2.1 read table.
///
/// `name` is resolved against a directory at use, not stored as a pointer into
/// a temporary: `AssetPath(...).c_str()` in an aggregate initializer would
/// dangle the moment the initialization's full expression ended.
struct ReadPin {
    const char* name;         ///< Fixture file name.
    bool in_data_dir;         ///< True for tests/data, false for the asset dir.
    unsigned int dim[3];      ///< Node count per axis.
    double cell[3];           ///< Declared cell edges (Angstroms).
    double spacing[3];        ///< Node interval per axis (Angstroms).
    unsigned int space_group; ///< ISPG as the carrier reports it.
    std::size_t symop_count;  ///< Number of 80-byte records in the block.
    const char* symop_first;  ///< First record, stripped.
    const char* symop_second; ///< Second record, stripped; "" when there is one.
    double node0[3];          ///< Position of element 0 after read_map.
};

std::vector<std::string> SplitLines(const std::string& text) {
    std::vector<std::string> lines;
    std::size_t start = 0;
    while (start <= text.size() && !text.empty()) {
        const std::size_t end = text.find('\n', start);
        if (end == std::string::npos) {
            lines.push_back(text.substr(start));
            break;
        }
        lines.push_back(text.substr(start, end - start));
        start = end + 1;
    }
    return lines;
}

}  // namespace

TEST(MapIoReadTest, PinsTheFiveFixtures) {
    const ReadPin pins[] = {
        {"1d26_2fofc.ccp4", false,
         {49u, 49u, 25u},
         {43.3, 43.3, 24.52},
         {0.902083, 0.902083, 1.021667},
         96u, 8u, "x,y,z", "-y+1/2,x+1/2,z+3/4",
         {0.0, 0.0, 0.0}},
        {"340d_2fofc.ccp4", false,
         {61u, 61u, 37u},
         {43.18, 43.18, 25.12},
         {0.719667, 0.719667, 0.697778},
         96u, 8u, "x,y,z", "-y+1/2,x+1/2,z+3/4",
         {0.0, 0.0, 0.0}},
        {"3q9g_2fofc.ccp4", false,
         {37u, 37u, 61u},
         {32.867, 32.867, 55.413},
         {0.912972, 0.912972, 0.92355},
         98u, 16u, "x,y,z", "-y,x+1/2,z+1/4",
         {0.0, 0.0, 0.0}},
        {"390_emd_30342_A_z4.mrc", false,
         {90u, 82u, 64u},
         {99.057, 90.153, 70.119},
         {1.113, 1.113, 1.113},
         0u, 0u, "", "",
         {145.825, 112.825, 120.517}},
        // test_map.ccp4 is load-bearing: the only fixture whose NC differs from
        // its NX, and the only one with a nonzero NxSTART. It is also the only
        // one that can tell cell / NX apart from cell / (dim - 1).
        //
        // Its ISPG is 146 even though tests/data/make_test_map.py:50 calls
        // SetSpaceGroup(1), which section 9 measures to write ISPG 1: the
        // committed file predates the generator's current form. Pin the 146 the
        // file carries. Regenerating it to agree with its generator would move
        // the space-group pins underneath the tests that read them.
        {"test_map.ccp4", true,
         {21u, 21u, 21u},
         {21.0, 21.0, 21.0},
         {0.5, 0.5, 0.5},
         146u, 0u, "", "",
         {-5.0, -5.0, -5.0}},
    };

    for (const ReadPin& pin : pins) {
        SCOPED_TRACE(pin.name);
        const std::string path =
            pin.in_data_dir ? DataPath(pin.name) : AssetPath(pin.name);
        const MapFile map = read_map(path);
        ASSERT_TRUE(map.grid);

        EXPECT_EQ(map.grid->GetXDim(), pin.dim[0]);
        EXPECT_EQ(map.grid->GetYDim(), pin.dim[1]);
        EXPECT_EQ(map.grid->GetZDim(), pin.dim[2]);

        const UnitCellParams cell = get_unit_cell(*map.grid);
        EXPECT_NEAR(cell.a, pin.cell[0], 1e-3);
        EXPECT_NEAR(cell.b, pin.cell[1], 1e-3);
        EXPECT_NEAR(cell.c, pin.cell[2], 1e-3);
        EXPECT_NEAR(cell.alpha, 90.0, 1e-4);
        EXPECT_NEAR(cell.beta, 90.0, 1e-4);
        EXPECT_NEAR(cell.gamma, 90.0, 1e-4);

        const GridParams gp = get_grid_params(*map.grid);
        EXPECT_NEAR(gp.x_spacing, pin.spacing[0], 1e-5);
        EXPECT_NEAR(gp.y_spacing, pin.spacing[1], 1e-5);
        EXPECT_NEAR(gp.z_spacing, pin.spacing[2], 1e-5);

        EXPECT_NEAR(gp.x_origin, pin.node0[0], 1e-3);
        EXPECT_NEAR(gp.y_origin, pin.node0[1], 1e-3);
        EXPECT_NEAR(gp.z_origin, pin.node0[2], 1e-3);

        EXPECT_EQ(map.grid->GetSpaceGroup(), pin.space_group);

        if (pin.symop_count == 0u) {
            EXPECT_TRUE(map.symops.empty());
        } else {
            const std::vector<std::string> lines = SplitLines(map.symops);
            ASSERT_EQ(lines.size(), pin.symop_count);
            EXPECT_EQ(lines[0], pin.symop_first);
            EXPECT_EQ(lines[1], pin.symop_second);
            EXPECT_EQ(SymOp::ParseAll(map.symops).size(), pin.symop_count);
        }
    }
}

TEST(MapIoReadTest, PlacesTheEmMapAtItsOriginRecord) {
    // The one fixture with a nonzero ORIGIN. OEReadGrid never consults that
    // record (section 2.1), so a bare read lands 145 A away; this is the whole
    // reason read_map exists.
    OESystem::OESkewGrid bare;
    ASSERT_TRUE(OESystem::OEReadGrid(AssetPath("390_emd_30342_A_z4.mrc"), bare));
    const GridParams bare_gp = get_grid_params(bare);
    EXPECT_NEAR(bare_gp.x_origin, 0.0, 1e-3);

    const MapFile map = read_map(AssetPath("390_emd_30342_A_z4.mrc"));
    const GridParams gp = get_grid_params(*map.grid);
    EXPECT_NEAR(gp.x_origin, 145.825, 1e-3);
    EXPECT_NEAR(gp.y_origin, 112.825, 1e-3);
    EXPECT_NEAR(gp.z_origin, 120.517, 1e-3);
}

TEST(MapIoReadTest, MovesRatherThanResamplesTheEmPayload) {
    // The ORIGIN correction adds a constant to GetMid, so every voxel value is
    // untouched. A resample would change them.
    OESystem::OESkewGrid bare;
    ASSERT_TRUE(OESystem::OEReadGrid(AssetPath("390_emd_30342_A_z4.mrc"), bare));
    const MapFile map = read_map(AssetPath("390_emd_30342_A_z4.mrc"));

    ASSERT_EQ(map.grid->GetSize(), bare.GetSize());
    const float* moved = map.grid->GetValues();
    const float* original = bare.GetValues();
    for (unsigned int i = 0; i < bare.GetSize(); ++i) {
        ASSERT_EQ(moved[i], original[i]) << "voxel " << i;
    }
}
