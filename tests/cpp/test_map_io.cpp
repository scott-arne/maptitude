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
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>
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

/// Copy a fixture into a scratch file with named header words rewritten.
///
/// The variants section 2.1 measured (permuted MAPC/MAPR/MAPS, ORIGIN and
/// NxSTART in each combination) are cheap to rebuild this way, so each
/// measurement read_map depends on becomes a pin rather than a number that
/// lives only in the design document.
class HeaderVariant {
public:
    explicit HeaderVariant(const std::string& source)
        : path_(::testing::TempDir() + "/maptitude_variant_" +
                std::to_string(++counter_) + ".ccp4") {
        std::ifstream in(source, std::ios::binary);
        bytes_.assign(std::istreambuf_iterator<char>(in),
                      std::istreambuf_iterator<char>());
    }

    ~HeaderVariant() { std::remove(path_.c_str()); }

    HeaderVariant& SetInt(const std::size_t word, const std::int32_t value) {
        std::memcpy(&bytes_[(word - 1u) * 4u], &value, 4);
        return *this;
    }

    HeaderVariant& SetFloat(const std::size_t word, const float value) {
        std::memcpy(&bytes_[(word - 1u) * 4u], &value, 4);
        return *this;
    }

    /// Write the variant out and return its path.
    const std::string& Write() {
        std::ofstream out(path_, std::ios::binary);
        out.write(bytes_.data(), static_cast<std::streamsize>(bytes_.size()));
        out.close();
        return path_;
    }

    /// Splice a symmetry block in at offset 1024 and update NSYMBT to match.
    ///
    /// The payload moves rather than being overwritten, which is what the
    /// on-disk layout actually does. Truncating it instead would make
    /// OEReadGrid fail first and the malformed-symop cases below would report
    /// GridError, never reaching the parser they exist to exercise.
    HeaderVariant& SetSymopBlock(const std::string& block) {
        std::int32_t existing = 0;
        std::memcpy(&existing, &bytes_[(24 - 1) * 4], 4);
        bytes_.replace(1024u, static_cast<std::size_t>(existing), block);
        const std::int32_t updated = static_cast<std::int32_t>(block.size());
        std::memcpy(&bytes_[(24 - 1) * 4], &updated, 4);
        return *this;
    }

    /// Cut the file short, for the case where NSYMBT outruns the file.
    HeaderVariant& TruncateTo(const std::size_t size) {
        bytes_.resize(size);
        return *this;
    }

private:
    static int counter_;
    std::string path_;
    std::string bytes_;
};

int HeaderVariant::counter_ = 0;

/// Build a big-endian copy of a little-endian CCP4 file.
///
/// Words 1-52 and 55-56 are numeric and swap; word 53 is the "MAP " string and
/// words 57-256 are the ten 80-character labels, both text, so neither does.
/// Word 54 is MACHST, which is set to the big-endian stamp rather than swapped.
/// The symmetry block is text. The payload is float32 and swaps.
std::string MakeBigEndianCopy(const std::string& source) {
    std::ifstream in(source, std::ios::binary);
    std::string bytes((std::istreambuf_iterator<char>(in)),
                      std::istreambuf_iterator<char>());

    auto swap_word = [&bytes](const std::size_t word) {
        std::uint32_t value = 0u;
        std::memcpy(&value, &bytes[(word - 1u) * 4u], 4);
        value = ((value & 0x000000FFu) << 24) | ((value & 0x0000FF00u) << 8) |
                ((value & 0x00FF0000u) >> 8) | ((value & 0xFF000000u) >> 24);
        std::memcpy(&bytes[(word - 1u) * 4u], &value, 4);
    };

    std::int32_t nsymbt = 0;
    std::memcpy(&nsymbt, &bytes[(24 - 1) * 4], 4);

    for (std::size_t word = 1; word <= 52; ++word) {
        swap_word(word);
    }
    for (std::size_t word = 55; word <= 56; ++word) {
        swap_word(word);
    }
    bytes[(54 - 1) * 4 + 0] = static_cast<char>(0x11);
    bytes[(54 - 1) * 4 + 1] = static_cast<char>(0x11);
    bytes[(54 - 1) * 4 + 2] = static_cast<char>(0x00);
    bytes[(54 - 1) * 4 + 3] = static_cast<char>(0x00);

    const std::size_t payload_start = 1024u + static_cast<std::size_t>(nsymbt);
    for (std::size_t offset = payload_start; offset + 4u <= bytes.size();
         offset += 4u) {
        std::uint32_t value = 0u;
        std::memcpy(&value, &bytes[offset], 4);
        value = ((value & 0x000000FFu) << 24) | ((value & 0x0000FF00u) << 8) |
                ((value & 0x00FF0000u) >> 8) | ((value & 0xFF000000u) >> 24);
        std::memcpy(&bytes[offset], &value, 4);
    }

    const std::string path =
        ::testing::TempDir() + "/maptitude_bigendian_1d26.ccp4";
    std::ofstream out(path, std::ios::binary);
    out.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    return path;
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

TEST(MapIoReadTest, RejectsNsymbtExceedingRecordCap) {
    // Regression test for the symop record count cap. A header declaring more
    // than MAX_SYMOP_RECORDS (4096) must not drive an unbounded allocation.
    //
    // NSYMBT = 409600 is 5120 records, exceeding the cap. Measurement showed
    // OEReadGrid accepts this (7 ms) and the cap check fires. Without the cap,
    // read_map would allocate 400 KB; with larger values a 2 GB file declaring
    // 2 GB of symops would allocate that. The message check ensures this test
    // fails if the cap is removed.
    HeaderVariant variant(DataPath("test_map.ccp4"));
    variant.SetInt(24, 409600);

    try {
        MapFile map = read_map(variant.Write());
        FAIL() << "Expected GridError for NSYMBT exceeding record cap";
    } catch (const GridError& e) {
        const std::string msg(e.what());
        EXPECT_NE(msg.find("exceeds cap"), std::string::npos)
            << "Expected cap-check message, got: " << msg;
        EXPECT_NE(msg.find("5120 symmetry records"), std::string::npos)
            << "Expected cap-check message, got: " << msg;
        EXPECT_NE(msg.find("4096 records"), std::string::npos)
            << "Expected cap-check message, got: " << msg;
    }
}

TEST(MapIoReadTest, RejectsNsymbtLargerThanFile) {
    // Regression test for the NSYMBT file-size bound. A header declaring
    // NSYMBT larger than the file can contain must not drive an allocation off
    // the declared size; ReadSymopBlock must check the file's actual size first.
    //
    // NSYMBT = 40000 is 500 records (under the cap) but test_map.ccp4 is only
    // 38068 bytes, so the file ends before the symmetry block would. Measurement
    // showed OEReadGrid accepts this (9 ms) and the file-size check fires.
    // Without the bound, read_map would allocate 40 KB and throw after a short
    // read; with it, read_map throws before allocating, naming both sizes. The
    // message check ensures this test fails if the bound is removed.
    HeaderVariant variant(DataPath("test_map.ccp4"));
    variant.SetInt(24, 40000);

    try {
        MapFile map = read_map(variant.Write());
        FAIL() << "Expected GridError for NSYMBT exceeding file size";
    } catch (const GridError& e) {
        const std::string msg(e.what());
        EXPECT_NE(msg.find("file size is"), std::string::npos)
            << "Expected size-check message, got: " << msg;
        EXPECT_NE(msg.find("need at least"), std::string::npos)
            << "Expected size-check message, got: " << msg;
    }
}

TEST(MapIoHeaderVariantTest, HonoursAPermutedAxisOrder) {
    // MAPC/MAPR/MAPS at words 17-19. All five fixtures ship (1, 2, 3), so a
    // synthetic variant is the only evidence the reader honours a permuted
    // order.
    //
    // The committed NxSTART is (-10, -10, -10) and the grid is cubic, which
    // makes a permutation the identity on node 0: measured, the permuted and
    // unpermuted files both read back at (-5, -5, -5). So the starts have to
    // be made distinct first, or this variant tests nothing.
    //
    // Measured with NxSTART (-10, -4, 6) at spacing 0.5:
    //   MAPC/MAPR/MAPS  node 0
    //   (1, 2, 3)       (-5, -2, +3)   -- unpermuted
    //   (3, 1, 2)       (-2, +3, -5)
    //   (2, 3, 1)       (+3, -5, -2)
    //   (1, 3, 2)       (-5, +3, -2)
    // The rule the numbers show: word 17 names the crystal axis the file's
    // fast axis belongs to, so NCSTART's contribution lands on that axis.
    // Pin the placement, not an inequality.
    HeaderVariant variant(DataPath("test_map.ccp4"));
    variant.SetInt(5, -10).SetInt(6, -4).SetInt(7, 6);
    variant.SetInt(17, 3).SetInt(18, 1).SetInt(19, 2);
    const MapFile map = read_map(variant.Write());
    const GridParams gp = get_grid_params(*map.grid);

    EXPECT_NEAR(gp.x_origin, -2.0, 1e-9);
    EXPECT_NEAR(gp.y_origin, 3.0, 1e-9);
    EXPECT_NEAR(gp.z_origin, -5.0, 1e-9);

    // The fixture is cubic, so dims and spacing are unchanged by the
    // permutation and cannot testify either way -- assert them as the
    // invariants they are, not as evidence of the permutation.
    EXPECT_EQ(map.grid->GetXDim(), 21u);
    EXPECT_NEAR(gp.x_spacing, 0.5, 1e-9);
}

TEST(MapIoHeaderVariantTest, DistinctStartsAreWhatMakeThePermutationVisible) {
    // The premise of the test above, pinned so a later change to the fixture
    // cannot silently turn it into a tautology. On the committed equal starts
    // the permutation is unobservable in node 0; the test above is only
    // meaningful because it sets distinct ones.
    HeaderVariant permuted(DataPath("test_map.ccp4"));
    permuted.SetInt(17, 3).SetInt(18, 1).SetInt(19, 2);
    const GridParams permuted_gp =
        get_grid_params(*read_map(permuted.Write()).grid);
    const GridParams plain_gp =
        get_grid_params(*read_map(DataPath("test_map.ccp4")).grid);

    EXPECT_EQ(permuted_gp.x_origin, plain_gp.x_origin);
    EXPECT_EQ(permuted_gp.y_origin, plain_gp.y_origin);
    EXPECT_EQ(permuted_gp.z_origin, plain_gp.z_origin);
}

TEST(MapIoHeaderVariantTest, OriginRecordWinsWhenNxStartIsAlsoSet) {
    // test_map.ccp4 already carries NxSTART (-10, -10, -10), which the reader
    // applies as node 0 = (-5, -5, -5). Adding a nonzero ORIGIN makes both
    // records live, which is exactly the tiebreak case.
    HeaderVariant variant(DataPath("test_map.ccp4"));
    variant.SetFloat(50, 7.0f).SetFloat(51, 8.0f).SetFloat(52, 9.0f);
    const std::string path = variant.Write();

    const MapFile by_origin = read_map(path, OriginSource::ORIGIN_RECORD);
    const GridParams origin_gp = get_grid_params(*by_origin.grid);
    EXPECT_NEAR(origin_gp.x_origin, 7.0, 1e-6);
    EXPECT_NEAR(origin_gp.y_origin, 8.0, 1e-6);
    EXPECT_NEAR(origin_gp.z_origin, 9.0, 1e-6);

    const MapFile by_nxstart = read_map(path, OriginSource::NXSTART);
    const GridParams nxstart_gp = get_grid_params(*by_nxstart.grid);
    EXPECT_NEAR(nxstart_gp.x_origin, -5.0, 1e-6);
    EXPECT_NEAR(nxstart_gp.y_origin, -5.0, 1e-6);
    EXPECT_NEAR(nxstart_gp.z_origin, -5.0, 1e-6);
}

TEST(MapIoHeaderVariantTest, TiebreakIsIgnoredWhenOnlyOneRecordIsSet) {
    // NxSTART alone: both tiebreaks leave the reader's placement alone.
    const MapFile a = read_map(DataPath("test_map.ccp4"),
                               OriginSource::ORIGIN_RECORD);
    const MapFile b = read_map(DataPath("test_map.ccp4"),
                               OriginSource::NXSTART);
    EXPECT_NEAR(get_grid_params(*a.grid).x_origin,
                get_grid_params(*b.grid).x_origin, 1e-9);
    EXPECT_NEAR(get_grid_params(*a.grid).x_origin, -5.0, 1e-6);

    // ORIGIN alone: both tiebreaks apply it.
    HeaderVariant variant(DataPath("test_map.ccp4"));
    variant.SetInt(5, 0).SetInt(6, 0).SetInt(7, 0);
    variant.SetFloat(50, 3.5f).SetFloat(51, 3.5f).SetFloat(52, 3.5f);
    const std::string path = variant.Write();
    EXPECT_NEAR(get_grid_params(*read_map(path, OriginSource::NXSTART).grid)
                    .x_origin,
                3.5, 1e-6);
    EXPECT_NEAR(
        get_grid_params(*read_map(path, OriginSource::ORIGIN_RECORD).grid)
            .x_origin,
        3.5, 1e-6);
}

TEST(MapIoHeaderVariantTest, ResolvesLittleEndianForEveryShippedFixture) {
    // Every committed fixture reads correctly on this little-endian host, so
    // none of them can be carrying a big-endian MACHST. Pinning that keeps a
    // future fixture from silently taking the swap path.
    const char* fixtures[] = {"1d26_2fofc.ccp4", "340d_2fofc.ccp4",
                              "3q9g_2fofc.ccp4", "390_emd_30342_A_z4.mrc"};
    for (const char* name : fixtures) {
        SCOPED_TRACE(name);
        std::ifstream in(AssetPath(name), std::ios::binary);
        char header[1024] = {0};
        in.read(header, 1024);
        EXPECT_NE(static_cast<unsigned char>(header[(54 - 1) * 4]), 0x11u);
    }
}

TEST(MapIoErrorTest, RaisesOnAMissingFile) {
    EXPECT_THROW(read_map(DataPath("no_such_map.ccp4")), GridError);
}

TEST(MapIoErrorTest, RaisesOnATruncatedFile) {
    const std::string path =
        ::testing::TempDir() + "/maptitude_truncated.ccp4";
    {
        std::ofstream out(path, std::ios::binary);
        out << "not a map";
    }
    EXPECT_THROW(read_map(path), GridError);
    std::remove(path.c_str());
}

TEST(MapIoErrorTest, RaisesWhenNsymbtIsNotAMultipleOfEighty) {
    HeaderVariant variant(DataPath("test_map.ccp4"));
    variant.SetInt(24, 37);
    EXPECT_THROW(read_map(variant.Write()), GridError);
}

TEST(MapIoErrorTest, RaisesWhenNsymbtIsNegative) {
    HeaderVariant variant(DataPath("test_map.ccp4"));
    variant.SetInt(24, -80);

    try {
        MapFile map = read_map(variant.Write());
        FAIL() << "Expected GridError for negative NSYMBT";
    } catch (const GridError& e) {
        const std::string msg(e.what());
        EXPECT_NE(msg.find("not a non-negative multiple of"), std::string::npos)
            << "Expected guard 1 message, got: " << msg;
    }
}

TEST(MapIoErrorTest, RaisesWhenTheFileEndsBeforeTheSymopBlockDoes) {
    HeaderVariant variant(DataPath("test_map.ccp4"));
    variant.SetInt(24, 160).TruncateTo(1024u + 80u);
    EXPECT_THROW(read_map(variant.Write()), GridError);
}

TEST(MapIoErrorTest, RaisesSymOpErrorOnEachMalformedRecordShape) {
    // One record per malformed shape section 2.1 measured. Each is padded to
    // the fixed 80 bytes the on-disk block uses. A valid first record precedes
    // the bad one so the case shows the parser rejecting a record rather than
    // rejecting the block wholesale.
    const char* malformed[] = {
        "not a symop at all",  // non-triplet
        "x,y",                 // two components
        "x,y,z,w",             // four components
        "x,y,q*z",             // bad coefficient
    };
    for (const char* record : malformed) {
        SCOPED_TRACE(record);
        std::string block("x,y,z");
        block.resize(80, ' ');
        std::string second(record);
        second.resize(80, ' ');
        block += second;

        HeaderVariant variant(DataPath("test_map.ccp4"));
        variant.SetSymopBlock(block);
        EXPECT_THROW(read_map(variant.Write()), SymOpError);
    }
}

TEST(MapIoErrorTest, AnEmptySymopBlockYieldsNoOperators) {
    HeaderVariant variant(DataPath("test_map.ccp4"));
    variant.SetSymopBlock(std::string(80, ' '));
    const MapFile map = read_map(variant.Write());
    EXPECT_TRUE(map.symops.empty());
    EXPECT_TRUE(SymOp::ParseAll(map.symops).empty());
}

TEST(MapIoErrorTest, ReadsAWellFormedSymopBlockSplicedIntoAFixtureThatHadNone) {
    // The positive control for SetSymopBlock. Without it, a splice that
    // corrupted the file would make every malformed case above pass for the
    // wrong reason -- they only assert that something was rejected.
    std::string block("x,y,z");
    block.resize(80, ' ');
    std::string second("-x,-y,z+1/2");
    second.resize(80, ' ');
    block += second;

    HeaderVariant variant(DataPath("test_map.ccp4"));
    variant.SetSymopBlock(block);
    const MapFile map = read_map(variant.Write());
    EXPECT_EQ(map.symops, "x,y,z\n-x,-y,z+1/2");
    EXPECT_EQ(SymOp::ParseAll(map.symops).size(), 2u);
    // The payload moved rather than being overwritten, so the grid still reads.
    EXPECT_EQ(map.grid->GetSize(), 9261u);
}

TEST(MapIoEndiannessTest, ReadsABigEndianFileAsItsLittleEndianOriginal) {
    const std::string path = MakeBigEndianCopy(AssetPath("1d26_2fofc.ccp4"));
    const MapFile big = read_map(path);
    const MapFile little = read_map(AssetPath("1d26_2fofc.ccp4"));

    EXPECT_EQ(big.grid->GetXDim(), little.grid->GetXDim());
    EXPECT_EQ(big.grid->GetYDim(), little.grid->GetYDim());
    EXPECT_EQ(big.grid->GetZDim(), little.grid->GetZDim());

    const GridParams big_gp = get_grid_params(*big.grid);
    const GridParams little_gp = get_grid_params(*little.grid);
    EXPECT_NEAR(big_gp.x_spacing, little_gp.x_spacing, 1e-6);
    EXPECT_NEAR(big_gp.y_spacing, little_gp.y_spacing, 1e-6);
    EXPECT_NEAR(big_gp.z_spacing, little_gp.z_spacing, 1e-6);
    EXPECT_NEAR(big_gp.x_origin, little_gp.x_origin, 1e-6);
    EXPECT_NEAR(big_gp.y_origin, little_gp.y_origin, 1e-6);
    EXPECT_NEAR(big_gp.z_origin, little_gp.z_origin, 1e-6);

    EXPECT_EQ(big.symops, little.symops);

    ASSERT_EQ(big.grid->GetSize(), little.grid->GetSize());
    const float* big_values = big.grid->GetValues();
    const float* little_values = little.grid->GetValues();
    for (unsigned int i = 0; i < little.grid->GetSize(); ++i) {
        ASSERT_EQ(big_values[i], little_values[i]) << "voxel " << i;
    }

    std::remove(path.c_str());
}
