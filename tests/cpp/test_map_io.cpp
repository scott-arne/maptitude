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
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <locale>
#include <string>
#include <system_error>
#include <vector>

#include <sys/stat.h>
#include <unistd.h>

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
                std::to_string(::getpid()) + "_" +
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

/// Build a big-endian copy of a little-endian CCP4 file with RAII cleanup.
///
/// Words 1-52 and 55-56 are numeric and swap; word 53 is the "MAP " string and
/// words 57-256 are the ten 80-character labels, both text, so neither does.
/// Word 54 is MACHST, which is set to the big-endian stamp rather than swapped.
/// The symmetry block is text. The payload is float32 and swaps.
class BigEndianCopy {
public:
    explicit BigEndianCopy(const std::string& source)
        : path_(::testing::TempDir() + "/maptitude_bigendian_" +
                std::to_string(::getpid()) + "_" +
                std::to_string(++counter_) + ".ccp4") {
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

        std::ofstream out(path_, std::ios::binary);
        out.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    }

    ~BigEndianCopy() { std::remove(path_.c_str()); }

    const std::string& Path() const { return path_; }

private:
    static int counter_;
    std::string path_;
};

int BigEndianCopy::counter_ = 0;

/// RAII wrapper for a truncated test file.
class TruncatedFile {
public:
    TruncatedFile()
        : path_(::testing::TempDir() + "/maptitude_truncated_" +
                std::to_string(::getpid()) + "_" +
                std::to_string(++counter_) + ".ccp4") {
        std::ofstream out(path_, std::ios::binary);
        out << "not a map";
    }

    ~TruncatedFile() { std::remove(path_.c_str()); }

    const std::string& Path() const { return path_; }

private:
    static int counter_;
    std::string path_;
};

int TruncatedFile::counter_ = 0;

/// A scratch destination that removes itself, so a failing assertion cannot
/// leave a stale map behind for the next run to read.
///
/// The process id is part of the name because gtest_discover_tests gives every
/// case its own ctest process, where a static counter alone repeats.
class ScratchPath {
public:
    explicit ScratchPath(const std::string& extension)
        : path_(::testing::TempDir() + "/maptitude_write_" +
                std::to_string(::getpid()) + "_" +
                std::to_string(++counter_) + extension) {
        std::remove(path_.c_str());
    }
    ~ScratchPath() { std::remove(path_.c_str()); }
    const std::string& Str() const { return path_; }

private:
    static int counter_;
    std::string path_;
};

int ScratchPath::counter_ = 0;

/// Assert that a written file reads back as the grid it came from.
void ExpectRoundTrip(const OESystem::OESkewGrid& original,
                     const std::string& path,
                     const std::string& symops) {
    const MapFile back = read_map(path);
    ASSERT_TRUE(back.grid);
    EXPECT_EQ(compare_written_map(original, *back.grid), "");
    EXPECT_EQ(back.symops, symops);
}

/// Read @p count bytes of @p path starting at @p offset.
///
/// The symmetry-block test has to read the file as bytes. read_map strips and
/// rejoins the records, so comparing through it passes on a block whose NSYMBT
/// or record width is wrong but which strips to the same triplets, and those
/// bytes are what a consumer that is not read_map will read.
std::string ReadBytesAt(const std::string& path, const std::size_t offset,
                        const std::size_t count) {
    std::ifstream in(path, std::ios::binary);
    EXPECT_TRUE(in.good()) << "cannot open " << path;
    in.seekg(static_cast<std::streamoff>(offset));
    std::string bytes(count, '\0');
    in.read(&bytes[0], static_cast<std::streamsize>(count));
    EXPECT_EQ(in.gcount(), static_cast<std::streamsize>(count))
        << path << " holds fewer than " << (offset + count) << " bytes";
    return bytes;
}

/// Count the files beside @p destination whose names have the shape write_map's
/// temporary takes in *this* process: ".maptitude-" + <pid> + "-" + <hex> +
/// extension, in the destination's own directory.
///
/// The hex comes from std::random_device, so the name cannot be predicted and
/// the shape has to be matched instead. Enumerating the parent and requiring it
/// empty would not do: TempDir() is shared with the other scratch helpers in
/// this file, and under `ctest -j 8` with other processes' files as well.
///
/// The pattern carries no destination stem, because write_map's temporary
/// carries none either, so this counts every temporary of this process and of
/// the matching extension in the directory rather than only the ones belonging
/// to @p destination. Those coincide here: this file makes no concurrent
/// write_map call, so no other destination's temporary of this process is live
/// when a count runs.
///
/// The pid in the pattern is what keeps another process's temporary out. A
/// leftover from a run that was killed or that crashed does not disappear when
/// its process does, so without the pid it would be counted by every later
/// assertion whose destination carried the same extension, in a serial run as
/// readily as under `ctest -j 8`.
std::size_t CountTemporarySiblings(const std::string& destination) {
    const std::filesystem::path dest(destination);
    const std::string prefix =
        ".maptitude-" + std::to_string(::getpid()) + "-";
    const std::string suffix = dest.extension().string();

    std::error_code error;
    std::filesystem::directory_iterator entries(dest.parent_path(), error);
    EXPECT_FALSE(error) << "cannot enumerate " << dest.parent_path().string()
                        << ": " << error.message();

    std::size_t found = 0;
    for (const std::filesystem::directory_entry& entry : entries) {
        const std::string name = entry.path().filename().string();
        if (name.size() > prefix.size() + suffix.size() &&
            name.compare(0, prefix.size(), prefix) == 0 &&
            name.compare(name.size() - suffix.size(), suffix.size(), suffix) ==
                0) {
            ++found;
        }
    }
    return found;
}

/// A file with the name shape CountTemporarySiblings looks for, planted beside
/// a destination so a test can show that counter is capable of seeing one.
///
/// The pid sits where write_map's temporary carries it, so the decoy is of the
/// shape this process's own temporary takes. That also keeps two processes
/// planting a decoy in a shared TempDir from removing it from under each other.
class DecoySibling {
public:
    explicit DecoySibling(const std::string& destination) {
        const std::filesystem::path dest(destination);
        path_ = (dest.parent_path() /
                 (".maptitude-" + std::to_string(::getpid()) + "-decoy" +
                  dest.extension().string()))
                    .string();
        std::ofstream out(path_, std::ios::binary);
        out << "decoy";
    }
    ~DecoySibling() { std::remove(path_.c_str()); }
    const std::string& Str() const { return path_; }

private:
    std::string path_;
};

/// A file of the shape write_map's temporary takes in a *different* process,
/// planted beside a destination so a test can show this process's leak count
/// does not attribute it here.
///
/// The pid written into the name is 0, which no test process holds. That is
/// what makes the file safe to plant in a TempDir shared under `ctest -j`: a
/// concurrent test process scoping its own count to its own pid cannot match
/// it, so this file's presence cannot fail somebody else's leak assertion.
class StaleSibling {
public:
    explicit StaleSibling(const std::string& destination) {
        const std::filesystem::path dest(destination);
        path_ = (dest.parent_path() /
                 (".maptitude-0-stale" + dest.extension().string()))
                    .string();
        std::ofstream out(path_, std::ios::binary);
        out << "stale";
    }
    ~StaleSibling() { std::remove(path_.c_str()); }
    const std::string& Str() const { return path_; }

private:
    std::string path_;
};

/// A numpunct that groups every digit and separates the groups with '/'.
///
/// Grouping is what a global locale does to an integer written through an
/// ostream, and '/' is a separator no path component can carry. A name built
/// through a stream carrying this facet therefore resolves under a directory
/// that does not exist, which is what turns "the name was formatted through
/// the global locale" from an odd but usable filename into a create that
/// fails.
///
/// Built here rather than taken from the system, so the case using it does not
/// depend on which locales the runner has generated.
class SlashGrouping : public std::numpunct<char> {
protected:
    char do_thousands_sep() const override { return '/'; }
    std::string do_grouping() const override { return "\1"; }
};

/// Installs a global locale for a scope and puts the previous one back.
///
/// The restore runs from a destructor because a failing gtest assertion returns
/// out of the case, and a global locale left installed would follow every case
/// that runs after it in the same process.
class ScopedGlobalLocale {
public:
    explicit ScopedGlobalLocale(const std::locale& locale)
        : previous_(std::locale::global(locale)) {}
    ~ScopedGlobalLocale() { std::locale::global(previous_); }
    ScopedGlobalLocale(const ScopedGlobalLocale&) = delete;
    ScopedGlobalLocale& operator=(const ScopedGlobalLocale&) = delete;

private:
    std::locale previous_;
};

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
    TruncatedFile truncated;
    EXPECT_THROW(read_map(truncated.Path()), GridError);
}

TEST(MapIoErrorTest, RaisesOnAPathWithAnEmbeddedNul) {
    // read_map has no gate to escalate past the way write_map does, so this is
    // not the security case; it is the path contract. The prefix here names a
    // file that reads perfectly and the calls that open it stop at the NUL, so
    // this call returned that file under a name the caller never gave: with
    // the guard commented out it does not throw at all. The same split reached
    // the diagnostic. Against a prefix OpenEye refuses, the message named the
    // text up to the NUL and was cut off there, losing even its closing quote.
    const std::string path = DataPath("test_map.ccp4") + '\0' + ".ccp4";
    try {
        read_map(path);
        FAIL() << "expected GridError: the path carries an embedded NUL";
    } catch (const GridError& error) {
        const std::string message(error.what());
        EXPECT_NE(message.find("embedded NUL"), std::string::npos)
            << "refused by the wrong branch: " << message;
        EXPECT_EQ(message.find('\0'), std::string::npos)
            << "the message carries the NUL it is refusing, so whatever prints "
               "it will stop there";
    }
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
    BigEndianCopy big_endian(AssetPath("1d26_2fofc.ccp4"));
    const MapFile big = read_map(big_endian.Path());
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
}

TEST(MapIoEndiannessTest, SwapsTheOriginRecordOnABigEndianEmMap) {
    // 390_emd_30342_A_z4.mrc has nonzero ORIGIN (145.825, 112.825, 120.517)
    // and zero NCSTART, so this exercises read_map's own ORIGIN byte-swap at
    // src/MapIO.cpp:115 through WordAsFloat's Swap32. The 1d26 test above
    // compares (0,0,0) against (0,0,0) and cannot see a missing swap.
    BigEndianCopy big_endian(AssetPath("390_emd_30342_A_z4.mrc"));
    const MapFile little = read_map(AssetPath("390_emd_30342_A_z4.mrc"));
    const MapFile big = read_map(big_endian.Path());

    float big_x = 0.0f, big_y = 0.0f, big_z = 0.0f;
    ASSERT_TRUE(big.grid->ElementToSpatialCoord(0u, big_x, big_y, big_z));

    const GridParams little_gp = get_grid_params(*little.grid);
    EXPECT_NEAR(big_x, little_gp.x_origin, 1e-6);
    EXPECT_NEAR(big_y, little_gp.y_origin, 1e-6);
    EXPECT_NEAR(big_z, little_gp.z_origin, 1e-6);
    EXPECT_NEAR(big_x, 145.825, 1e-3);
    EXPECT_NEAR(big_y, 112.825, 1e-3);
    EXPECT_NEAR(big_z, 120.517, 1e-3);

    ASSERT_EQ(big.grid->GetSize(), little.grid->GetSize());
    const float* big_values = big.grid->GetValues();
    const float* little_values = little.grid->GetValues();
    for (unsigned int i = 0; i < little.grid->GetSize(); ++i) {
        ASSERT_EQ(big_values[i], little_values[i]) << "voxel " << i;
    }
}

TEST(MapIoWriteTest, RoundTripsEachOfTheFiveFixtures) {
    const char* fixtures[] = {"1d26_2fofc.ccp4", "340d_2fofc.ccp4",
                              "3q9g_2fofc.ccp4", "390_emd_30342_A_z4.mrc"};
    for (const char* name : fixtures) {
        SCOPED_TRACE(name);
        const MapFile source = read_map(AssetPath(name));
        ScratchPath out(".ccp4");
        write_map(out.Str(), *source.grid, source.symops);
        ExpectRoundTrip(*source.grid, out.Str(), source.symops);
    }

    SCOPED_TRACE("test_map.ccp4");
    const MapFile source = read_map(DataPath("test_map.ccp4"));
    ScratchPath out(".ccp4");
    write_map(out.Str(), *source.grid, source.symops);
    ExpectRoundTrip(*source.grid, out.Str(), source.symops);
}

TEST(MapIoWriteTest, WritesTheEmOriginIntoTheHeader) {
    // The one fixture whose node 0 exercises the ORIGIN patch at all.
    const MapFile source = read_map(AssetPath("390_emd_30342_A_z4.mrc"));
    ScratchPath out(".mrc");
    write_map(out.Str(), *source.grid, source.symops);

    std::ifstream in(out.Str(), std::ios::binary);
    char header[1024] = {0};
    in.read(header, 1024);
    float origin[3] = {0.0f, 0.0f, 0.0f};
    std::memcpy(origin, header + (50 - 1) * 4, 12);
    EXPECT_NEAR(origin[0], 145.825f, 1e-3f);
    EXPECT_NEAR(origin[1], 112.825f, 1e-3f);
    EXPECT_NEAR(origin[2], 120.517f, 1e-3f);
}

TEST(MapIoWriteTest, WritesTheSymmetryBlockAsFixedWidthRecordsOnDisk) {
    // The byte-level companion to RoundTripsEachOfTheFiveFixtures, and the test
    // side of write_map's raw symop check. The round-trip test compares
    // back.symops, which read_map has already stripped and rejoined: it passes
    // on a block whose NSYMBT is wrong, or whose records are newline-terminated
    // rather than 80-byte space-padded, because both strip to the same
    // triplets. A consumer that is not read_map reads these bytes.
    //
    // Measured on the two shipped fixtures that carry symmetry records
    // (390_emd_30342_A_z4.mrc and test_map.ccp4 both declare NSYMBT 0, so
    // neither can stand in here):
    //   1d26_2fofc.ccp4  NSYMBT  640   8 records
    //   3q9g_2fofc.ccp4  NSYMBT 1280  16 records
    // Both start "x,y,z" and end "y,x,-z", no record starts with a space, and
    // stripping each record and re-padding it to 80 with spaces reproduces the
    // on-disk block byte for byte. That last measurement is what lets this
    // demand equality with the source block rather than mere re-parsability.
    // Two fixtures, not one: a writer that hardcoded a record count would pass
    // on either alone.
    struct Fixture {
        const char* name;
        std::int32_t nsymbt;
        std::size_t records;
    };
    const Fixture fixtures[] = {{"1d26_2fofc.ccp4", 640, 8},
                                {"3q9g_2fofc.ccp4", 1280, 16}};

    for (const Fixture& fixture : fixtures) {
        SCOPED_TRACE(fixture.name);
        const std::string source_path = AssetPath(fixture.name);
        const MapFile source = read_map(source_path);
        ASSERT_FALSE(source.symops.empty())
            << "this fixture no longer carries symmetry records, so the "
               "comparisons below would hold vacuously";

        ScratchPath out(".ccp4");
        write_map(out.Str(), *source.grid, source.symops);

        const std::string header = ReadBytesAt(out.Str(), 0, 1024);
        std::int32_t nsymbt = 0;
        std::memcpy(&nsymbt, header.data() + (24 - 1) * 4, 4);
        EXPECT_EQ(nsymbt, fixture.nsymbt);

        const std::string written =
            ReadBytesAt(out.Str(), 1024, static_cast<std::size_t>(fixture.nsymbt));
        ASSERT_EQ(written.size() % 80u, 0u);
        EXPECT_EQ(written.size() / 80u, fixture.records);
        EXPECT_EQ(written.substr(0, 80),
                  std::string("x,y,z") + std::string(75, ' '));
        EXPECT_EQ(written.substr((fixture.records - 1u) * 80u, 80u),
                  std::string("y,x,-z") + std::string(74, ' '));

        // The shape assertions above name which property broke; this one is the
        // actual contract, and would otherwise fail as 640 bytes of diff.
        EXPECT_EQ(written,
                  ReadBytesAt(source_path, 1024,
                              static_cast<std::size_t>(fixture.nsymbt)));
    }
}

TEST(MapIoWriteTest, SplitsSemicolonSeparatedSymopsIntoOneRecordEach) {
    // write_map validates its symops through SymOp::ParseAll, which treats ';'
    // as an operator boundary as readily as '\n'. The canonicalization that
    // feeds the on-disk block has to agree, or the same two operators are one
    // 80-byte record in the ';' spelling and two in the '\n' spelling.
    //
    // This asserts on the raw bytes rather than through read_map, and on
    // hand-built symop text rather than on a fixture's. write_map's own NSYMBT
    // and block checks derive both sides from the canonicalization, so they
    // cannot see a canonicalization that is itself wrong; and read_map hands
    // the record back verbatim, so a one-record block with an embedded ';'
    // compares equal to the text that produced it. Fixture symops come from
    // read_map already split on newlines, which is the case that never fails.
    const std::string semicolons = "x,y,z;-x,y+1/2,-z";
    ASSERT_EQ(SymOp::ParseAll(semicolons).size(), 2u)
        << "the parser no longer reads this as two operators, so the split "
           "this test pins is not the one write_map validates against";

    const MapFile source = read_map(DataPath("test_map.ccp4"));
    ASSERT_TRUE(source.symops.empty())
        << "this fixture now carries its own symmetry records, so the block "
           "below would not be the one under test";

    ScratchPath out(".ccp4");
    write_map(out.Str(), *source.grid, semicolons);

    const std::string header = ReadBytesAt(out.Str(), 0, 1024);
    std::int32_t nsymbt = 0;
    std::memcpy(&nsymbt, header.data() + (24 - 1) * 4, 4);
    EXPECT_EQ(nsymbt, 160)
        << "the two operators were not written as two 80-byte records";

    const std::string written = ReadBytesAt(out.Str(), 1024, 160);
    EXPECT_EQ(written.substr(0, 80),
              std::string("x,y,z") + std::string(75, ' '));
    EXPECT_EQ(written.substr(80, 80),
              std::string("-x,y+1/2,-z") + std::string(69, ' '));
}

TEST(MapIoWriteTest, WritesToEachDocumentedExtension) {
    // The three spellings the header and the refusal message put in front of a
    // caller. They dispatch to the same CCP4 writer, so this is not a format
    // test: it is what shows the temporary carrying whichever extension it was
    // given rather than a hardcoded .ccp4. '.map' is here because it was named
    // in both places and written by no test.
    const MapFile source = read_map(AssetPath("1d26_2fofc.ccp4"));
    for (const char* extension : {".ccp4", ".mrc", ".map"}) {
        SCOPED_TRACE(extension);
        ScratchPath out(extension);
        write_map(out.Str(), *source.grid, source.symops);
        ExpectRoundTrip(*source.grid, out.Str(), source.symops);
    }
}

TEST(MapIoWriteTest, RoundTripsAGridWhoseNodeZeroIsAHalfIntegerOfSpacings) {
    // SetMid(0,0,0) with an even dim puts node 0 on a half integer, which
    // NCSTART cannot encode. The ORIGIN patch carries it and read_map prefers
    // ORIGIN, so the position survives. This is the regression test for the
    // integer-NCSTART refusal section 2.4 retired: if either half of that
    // mechanism regresses, this grid lands half a voxel out.
    OESystem::OESkewGrid grid;
    ASSERT_TRUE(grid.SetDim(20u, 20u, 20u));
    ASSERT_TRUE(grid.SetUnitCell(20.0f, 20.0f, 20.0f, 90.0f, 90.0f, 90.0f,
                                 40u, 40u, 40u));
    ASSERT_TRUE(grid.SetMid(0.0f, 0.0f, 0.0f));
    ASSERT_TRUE(grid.SetSpaceGroup(1u));
    std::vector<float> values(grid.GetSize());
    for (std::size_t i = 0; i < values.size(); ++i) {
        values[i] = static_cast<float>(i % 17u);
    }
    ASSERT_TRUE(grid.SetValues(values.data(),
                               static_cast<unsigned int>(values.size())));

    const GridParams before = get_grid_params(grid);
    ASSERT_NEAR(std::fmod(std::fabs(before.x_origin / before.x_spacing), 1.0),
                0.5, 1e-6)
        << "node 0 is not a half-integer number of spacings, so this test no "
           "longer exercises what it was written for";

    ScratchPath out(".ccp4");
    write_map(out.Str(), grid);
    ExpectRoundTrip(grid, out.Str(), "");
}

TEST(MapIoWriteTest, PreservesASpaceGroupTheInputCarries) {
    const MapFile source = read_map(AssetPath("3q9g_2fofc.ccp4"));
    ScratchPath out(".ccp4");
    write_map(out.Str(), *source.grid, source.symops);
    EXPECT_EQ(read_map(out.Str()).grid->GetSpaceGroup(), 98u);
}

TEST(MapIoWriteTest, DefaultsAnAbsentSpaceGroupToP1) {
    // Deliberately not identity: section 4 defaults to P1 because an unset
    // space group makes OEWriteGrid double the cell and regrid.
    const MapFile source = read_map(AssetPath("390_emd_30342_A_z4.mrc"));
    ASSERT_EQ(source.grid->GetSpaceGroup(), 0u);
    ScratchPath out(".mrc");
    write_map(out.Str(), *source.grid);
    const MapFile back = read_map(out.Str());
    EXPECT_TRUE(back.grid->HasSpaceGroup());
    EXPECT_EQ(back.grid->GetSpaceGroup(), 1u);
}

TEST(MapIoWriteTest, RefusesAGridWhoseCellEqualsItsSampledExtent) {
    // The geometric condition behind the wrap_and_pad_grid refusal, reached
    // without a molecule: no test under tests/cpp reads a molecule from the
    // assets, and this needs no new fixture plumbing to exercise the same
    // branch. MakeEmptyGrid declares cell = n * spacing with n divisions, so
    // the written NX equals the written NC and the re-read node count comes
    // back one higher per axis. The Python suite covers the same refusal
    // through wrap_and_pad_grid itself on a real asset.
    //
    // Measured on this exact construction, with the P1 default applied as
    // step 2 applies it: dim (21,21,21), cell 10.5, spacing 0.5 writes and
    // reads back as dim (22,22,22), size 9261 -> 10648. Both OEWriteGrid and
    // OEReadGrid return true; the header is self-consistent and describes a
    // different grid, which is what the verify exists to catch.
    OESystem::OESkewGrid grid = MakeEmptyGrid(5.0, 0.5);
    const GridParams gp = get_grid_params(grid);
    const UnitCellParams cell = get_unit_cell(grid);
    ASSERT_NEAR(cell.a, gp.x_dim * gp.x_spacing, 1e-6)
        << "this grid's cell is no longer its sampled extent, so it no longer "
           "reaches the branch this test exists for";

    ScratchPath out(".ccp4");
    try {
        write_map(out.Str(), grid);
        FAIL() << "expected GridError: the written header describes a 22^3 grid";
    } catch (const GridError& error) {
        EXPECT_NE(std::string(error.what()).find("dimensions differ"),
                  std::string::npos)
            << "refused by the wrong branch: " << error.what();
    }
    std::ifstream probe(out.Str(), std::ios::binary);
    EXPECT_FALSE(probe.good()) << "a refused write left a file behind";
}

TEST(MapIoWriteTest, AFailedWriteLeavesTheDestinationAlone) {
    const MapFile source = read_map(AssetPath("1d26_2fofc.ccp4"));
    ScratchPath out(".ccp4");
    write_map(out.Str(), *source.grid, source.symops);
    const MapFile before = read_map(out.Str());

    OESystem::OESkewGrid doomed = MakeEmptyGrid(5.0, 0.5);
    EXPECT_THROW(write_map(out.Str(), doomed), GridError);

    const MapFile after = read_map(out.Str());
    EXPECT_EQ(compare_written_map(*before.grid, *after.grid), "");
    EXPECT_EQ(after.symops, before.symops);
}

TEST(MapIoWriteTest, ARefusedWriteLeavesNoHiddenTemporaryBehind) {
    // The two refusal cases that reach the verify both assert on the
    // destination, which is not where the temporary is: write_map names it as a
    // hidden sibling in the same directory. So a TemporaryFile destructor that
    // stopped removing it would leak one file per refused write with every
    // other case in this file still green.
    //
    // Measured, rather than argued: neutering that destructor leaves MapIo at
    // 38/38 and two files behind in TempDir(), one per test process that threw
    // after the temporary was created. This case is the only thing that sees
    // it.
    OESystem::OESkewGrid doomed = MakeEmptyGrid(5.0, 0.5);
    ScratchPath out(".ccp4");

    // A counter that reports zero because it is pointed at the wrong directory,
    // or because it matches the wrong name shape, is indistinguishable from one
    // reporting a genuine absence. Show it can see a file of exactly the shape
    // and in exactly the place it is about to report the absence of.
    {
        const DecoySibling decoy(out.Str());
        ASSERT_EQ(CountTemporarySiblings(out.Str()), 1u)
            << "the sibling probe cannot see " << decoy.Str()
            << ", so its zero below would say nothing";
    }
    ASSERT_EQ(CountTemporarySiblings(out.Str()), 0u)
        << "the decoy outlived its scope, so the count below starts dirty";

    EXPECT_THROW(write_map(out.Str(), doomed), GridError);
    EXPECT_EQ(CountTemporarySiblings(out.Str()), 0u)
        << "a refused write left its temporary beside " << out.Str();
}

TEST(MapIoWriteTest, WritesPastASiblingOfTheTemporaryNameShape) {
    // write_map reserves its temporary's name with an exclusive create, so a
    // name already taken is skipped rather than truncated. This case does not
    // force that skip and does not claim to: the name carries 32 bits of
    // std::random_device, which is not injectable, so a decoy cannot be made to
    // collide with it. What it pins is the reachable half -- a file of the shape
    // the reservation draws from, sitting beside the destination, neither blocks
    // the write nor is written through.
    const MapFile source = read_map(DataPath("test_map.ccp4"));
    ScratchPath out(".ccp4");
    const DecoySibling decoy(out.Str());

    write_map(out.Str(), *source.grid, source.symops);
    ExpectRoundTrip(*source.grid, out.Str(), source.symops);

    EXPECT_EQ(ReadBytesAt(decoy.Str(), 0, 5), "decoy")
        << "the write went through " << decoy.Str();
    EXPECT_EQ(CountTemporarySiblings(out.Str()), 1u)
        << "the count beside " << out.Str()
        << " is not the planted decoy alone, so the write left a temporary";
}

TEST(MapIoWriteTest, LeavesTheDestinationAtTheModeAnOrdinaryCreateGives) {
    // The destination is published by renaming the temporary onto it, so the
    // mode the caller sees is the temporary's. write_map creates that temporary
    // with 0666 and lets the umask narrow it, which is what puts it where an
    // ordinary create lands. Creating it 0600 -- what mkstemps does -- would
    // publish every map readable only by its writer. That decision was recorded
    // in a comment on MakeTemporarySibling and asserted nowhere.
    //
    // The reference is an ordinary create in the same directory rather than a
    // hardcoded 0644, so the case says the same thing under whatever umask the
    // suite runs.
    const MapFile source = read_map(DataPath("test_map.ccp4"));

    ScratchPath reference(".reference");
    {
        std::ofstream ordinary(reference.Str(), std::ios::binary);
        ASSERT_TRUE(ordinary.good()) << "cannot create " << reference.Str();
        ordinary << "reference";
    }
    struct ::stat reference_info {};
    ASSERT_EQ(::stat(reference.Str().c_str(), &reference_info), 0)
        << "cannot stat " << reference.Str();

    ScratchPath out(".ccp4");
    write_map(out.Str(), *source.grid, source.symops);
    struct ::stat written_info {};
    ASSERT_EQ(::stat(out.Str().c_str(), &written_info), 0)
        << "cannot stat " << out.Str();

    EXPECT_EQ(written_info.st_mode & 07777u, reference_info.st_mode & 07777u)
        << "the written map's mode is not the one an ordinary create gives in "
           "this directory under this umask";
}

TEST(MapIoWriteTest, DoesNotAttributeAnotherProcessesTemporaryToThisWrite) {
    // Measured, not argued. With the temporary's name carrying no pid, a serial
    // run whose TemporaryFile cleanup had been neutered left four
    // '.maptitude-<hex>.ccp4' files behind, one per test process that threw
    // after creating one; a later test's leak assertion then reported 4 against
    // 0. Those four processes had already exited, so `ctest -j` was not what
    // exposed it and the window was not the milliseconds between OEWriteGrid
    // and the rename. A temporary that outlives its run is counted from then on
    // by every assertion whose destination carries the same extension.
    //
    // The pid in the name is what scopes the count to this process. The stale
    // file below is of exactly the counted shape apart from that pid.
    const MapFile source = read_map(DataPath("test_map.ccp4"));
    ScratchPath out(".ccp4");
    const StaleSibling stale(out.Str());

    // A count of zero at the end says nothing if the counter matches nothing at
    // all. Show it still sees this process's own shape in this directory.
    const std::size_t baseline = CountTemporarySiblings(out.Str());
    {
        const DecoySibling decoy(out.Str());
        ASSERT_EQ(CountTemporarySiblings(out.Str()), baseline + 1u)
            << "the sibling probe cannot see " << decoy.Str()
            << ", so its count below would say nothing";
    }

    write_map(out.Str(), *source.grid, source.symops);
    EXPECT_EQ(CountTemporarySiblings(out.Str()), 0u)
        << "the count attributed " << stale.Str()
        << " to this write, though its name carries another process's pid";
}

TEST(MapIoWriteTest, ReservesTheTemporaryBeforeHandingTheNameToOEWriteGrid) {
    // The reservation is an exclusive create, so a directory that cannot hold
    // the temporary fails at that create and the message names it. Without the
    // reservation the same call carried an unwritable name all the way to
    // OEWriteGrid, which reported only that it had failed -- so this separates
    // "the name was taken first" from "the name was merely chosen first".
    //
    // A parent directory that does not exist is the cheapest way to make the
    // create fail. OEWriteGrid returns false rather than aborting on such a
    // path, measured for a '.ccp4' name, so the pre-reservation behavior this
    // discriminates against is a catchable GridError and not a killed process.
    const MapFile source = read_map(DataPath("test_map.ccp4"));
    const std::filesystem::path absent =
        std::filesystem::path(::testing::TempDir()) /
        ("maptitude_absent_reserve_" + std::to_string(::getpid()));
    ASSERT_FALSE(std::filesystem::exists(absent))
        << "the parent directory exists, so the temporary's create would "
           "succeed: " << absent.string();
    const std::string out = (absent / "out.ccp4").string();

    try {
        write_map(out, *source.grid, source.symops);
        FAIL() << "expected GridError: the temporary cannot be created in a "
                  "directory that does not exist";
    } catch (const GridError& error) {
        const std::string message(error.what());
        EXPECT_NE(message.find(out), std::string::npos)
            << "the message does not name the destination: " << message;
        EXPECT_NE(message.find("temporary"), std::string::npos)
            << "the message does not name the temporary whose create failed: "
            << message;
        EXPECT_EQ(message.find("OEWriteGrid"), std::string::npos)
            << "the name reached OEWriteGrid, so it was not reserved ahead of "
               "it: " << message;
    }

    EXPECT_FALSE(std::filesystem::exists(absent))
        << "a refused write created the destination's parent directory";
}

TEST(MapIoWriteTest, BuildsTheTemporarysNameOutsideTheGlobalLocale) {
    // write_map builds its temporary's name in an ostringstream, where num_put
    // groups integers through the locale's numpunct. A stream left on the
    // global locale therefore follows whatever separator the process has
    // installed. Measured with that name built verbatim: en_US.UTF-8 gives
    // '.maptitude-5,656-de,adb,eef.ccp4' and de_DE.UTF-8 gives
    // '.maptitude-5.656-de.adb.eef.ccp4', against
    // '.maptitude-5656-deadbeef.ccp4' under C.
    //
    // Both halves of that matter. The de_DE name carries three dot components
    // into a shape MakeTemporarySibling measured without any, and neither name
    // is one CountTemporarySiblings above matches -- it builds its prefix with
    // std::to_string, which never groups. So under such a locale the counter
    // stops seeing the temporary a write made while still seeing the decoys
    // these cases plant, which carry a std::to_string pid too. Every one of the
    // eleven CountTemporarySiblings assertions in this file then passes without
    // having looked at a temporary the library made, in all three of the shapes
    // they take: the eight expecting 0 read zero having matched nothing, the
    // two decoy cases expecting 1 read the decoy alone, and the reservation
    // case's baseline + 1 still holds because both of its terms count decoys.
    // That last one is the trap -- it is the positive control those cases run
    // to show the counter can see something, and it keeps reporting success
    // while the counter is blind to the only file it exists to catch.
    //
    // A real grouping locale would not make that visible from here: the write
    // still succeeds and the vacuous count still reads zero. The facet
    // installed below groups with '/' instead, so the grouped name resolves
    // under a directory that does not exist and the create refuses it -- errno
    // 2, measured, against the errno 17 a taken name gives. That refusal is
    // what this case asserts the absence of.
    //
    // A pid of two digits or more takes a separator of its own, and the hex
    // beside it takes one per digit past its first. The assertion is on the
    // pid because the hex is drawn from random_device: a run that drew a
    // one-digit value would leave the name ungrouped without it.
    ASSERT_GT(::getpid(), 9)
        << "a single-digit pid takes no separator, so this case would have "
           "nothing certain for the facet to group";

    const MapFile source = read_map(DataPath("test_map.ccp4"));
    ScratchPath out(".ccp4");

    // The case discriminates over two links: grouping reaches the name, and a
    // grouped name cannot be created. The ASSERT above arms the first. This
    // arms the second, which otherwise holds only by inference -- if the
    // candidate stopped being composed as parent_path() / name, or if this
    // directory came to hold the chain a slash-bearing name resolves under, the
    // grouped create would succeed and this case would stay green with the
    // imbue deleted. Checked before the locale scope, because gtest formats its
    // own file and line through the global locale too.
    const std::filesystem::path slash_bearing =
        std::filesystem::path(out.Str()).parent_path() / ".maptitude-1/2-x.ccp4";
    std::ofstream probe(slash_bearing);
    ASSERT_FALSE(probe.is_open())
        << "a slash-bearing sibling can be created beside " << out.Str()
        << ", so a grouped name would not be refused there and this case would "
           "pass with the imbue removed";

    // The refusal is carried out of the scope rather than reported inside it,
    // because gtest formats its own file and line through the global locale
    // too: reported in place, the failure names 'test_map_io.cpp:1/2/8/3'.
    std::string refusal;
    {
        const ScopedGlobalLocale grouping(
            std::locale(std::locale::classic(), new SlashGrouping));
        try {
            write_map(out.Str(), *source.grid, source.symops);
        } catch (const GridError& error) {
            refusal = error.what();
        }
    }
    if (!refusal.empty()) {
        FAIL() << "the temporary's name was formatted through the global "
                  "locale, so its grouping separator reached the name: "
               << refusal;
    }

    // Read back outside the scope, so the locale under test is not also the
    // one this verification runs under.
    ExpectRoundTrip(*source.grid, out.Str(), source.symops);
    EXPECT_EQ(CountTemporarySiblings(out.Str()), 0u)
        << "the write left a temporary beside " << out.Str();
}

TEST(MapIoWriteTest, RoundTripsAPermutedAxisOrder) {
    // Section 9 leaves the write side of MAPC/MAPR/MAPS unmeasured: no measured
    // input makes OEWriteGrid emit a permuted header, and section 4 leaves
    // NC/NX/MAPC as OEWriteGrid wrote them. Reading a permuted file and writing
    // it back is the cheapest thing that exercises the ORIGIN patch against a
    // permuted input.
    //
    // If this fails, do not patch MAPC to make it pass -- that is the axis
    // recomputation section 3.1 forbids. Record the gap in the task report,
    // change the case to EXPECT_THROW(GridError) documenting the refusal, and
    // say so. A silently wrong permuted write is the outcome to avoid; a
    // refused one is acceptable.
    HeaderVariant variant(DataPath("test_map.ccp4"));
    variant.SetInt(17, 3).SetInt(18, 1).SetInt(19, 2);
    const MapFile source = read_map(variant.Write());

    ScratchPath out(".ccp4");
    write_map(out.Str(), *source.grid, source.symops);
    ExpectRoundTrip(*source.grid, out.Str(), source.symops);
}

TEST(MapIoWriteTest, RaisesOnAnExtensionOpenEyeDoesNotWrite) {
    // Section 2.4 measures the unguarded call to be a process abort, not a
    // false return. If the check in step 1 is missing or wrong, this test does
    // not fail -- it takes the whole binary down.
    const MapFile source = read_map(DataPath("test_map.ccp4"));
    for (const char* extension : {".dat", ""}) {
        SCOPED_TRACE(extension);
        ScratchPath out(extension);
        EXPECT_THROW(write_map(out.Str(), *source.grid), GridError);
        std::ifstream probe(out.Str(), std::ios::binary);
        EXPECT_FALSE(probe.good()) << "a refused write left a file behind";
    }
}

TEST(MapIoWriteTest, RefusesADestinationWithAnEmbeddedNul) {
    // The other side of the extension gate above. That gate reads
    // std::filesystem's view of the whole string, which reports the extension
    // of "victim.dat\0.ccp4" as ".ccp4" and admits it, while every filesystem
    // call downstream stops at the NUL. So the name the gate admitted and the
    // name the rename published were different names: measured before the
    // guard, this call returned successfully having replaced victim.dat with a
    // 38068-byte map. The defect's signature is that successful return, so the
    // bytes are what discriminates here and the refusal is recorded rather
    // than required, to keep the byte check running when the call does not
    // throw.
    const MapFile source = read_map(DataPath("test_map.ccp4"));
    ScratchPath victim(".dat");
    {
        std::ofstream out(victim.Str(), std::ios::binary);
        ASSERT_TRUE(out.good()) << "cannot create " << victim.Str();
        out << "ORIGINAL CONTENTS\n";
    }

    const std::string destination = victim.Str() + '\0' + ".ccp4";
    std::string refusal;
    try {
        write_map(destination, *source.grid, source.symops);
    } catch (const GridError& error) {
        refusal = error.what();
    }

    const std::string original = "ORIGINAL CONTENTS\n";
    std::ifstream probe(victim.Str(), std::ios::binary);
    ASSERT_TRUE(probe.good()) << "the write removed " << victim.Str();
    const std::string after((std::istreambuf_iterator<char>(probe)),
                            std::istreambuf_iterator<char>());
    // The size is asserted first only to keep the report readable: a map
    // written through the NUL differs in length, and comparing the contents
    // outright prints the whole 38068 bytes of it.
    ASSERT_EQ(after.size(), original.size())
        << victim.Str() << " was written through the NUL";
    EXPECT_EQ(after, original)
        << victim.Str() << " was written through the NUL";

    ASSERT_FALSE(refusal.empty()) << "the write was not refused";
    EXPECT_NE(refusal.find("embedded NUL"), std::string::npos)
        << "refused by the wrong branch: " << refusal;
    EXPECT_EQ(refusal.find('\0'), std::string::npos)
        << "the message carries the NUL it is refusing, so whatever prints it "
           "will stop there";
}

TEST(MapIoWriteTest, RaisesSymOpErrorBeforeTouchingTheFilesystem) {
    const MapFile source = read_map(DataPath("test_map.ccp4"));
    ScratchPath out(".ccp4");
    EXPECT_THROW(write_map(out.Str(), *source.grid, "not a symop"), SymOpError);
    std::ifstream probe(out.Str(), std::ios::binary);
    EXPECT_FALSE(probe.good());
}

TEST(MapIoWriteTest, RaisesCellErrorOnSamplingThatIsNotAxisAligned) {
    // The header documents CellError for this path. read_map can hand back a
    // grid with a non-90 cell angle, and neither write_map nor
    // compare_written_map can describe one: get_grid_params requires the
    // sampling axes to line up with the cartesian axes. Pin the type and the
    // clean refusal so the documented contract is not just prose.
    OESystem::OESkewGrid grid;
    ASSERT_TRUE(grid.SetDim(4u, 5u, 6u));
    ASSERT_TRUE(grid.SetUnitCell(4.0f, 5.0f, 6.0f, 90.0f, 90.0f, 120.0f,
                                 8u, 10u, 12u));
    ASSERT_TRUE(grid.SetMid(0.0f, 0.0f, 0.0f));
    ASSERT_TRUE(grid.SetSpaceGroup(1u));
    std::vector<float> values(grid.GetSize(), 1.0f);
    ASSERT_TRUE(grid.SetValues(values.data(),
                               static_cast<unsigned int>(values.size())));

    ScratchPath out(".ccp4");
    try {
        write_map(out.Str(), grid);
        FAIL() << "expected CellError: the sampling is not axis-aligned";
    } catch (const CellError& error) {
        EXPECT_NE(std::string(error.what()).find("axis-aligned"),
                  std::string::npos)
            << "the message does not name the axis alignment requirement: "
            << error.what();
    }

    std::ifstream probe(out.Str(), std::ios::binary);
    EXPECT_FALSE(probe.good()) << "a refused write left a file behind";
    EXPECT_EQ(CountTemporarySiblings(out.Str()), 0u)
        << "a refused write left its hidden temporary behind";
}

TEST(MapIoWriteTest, RefusesACompressedDestination) {
    // OEIsWriteableGrid accepts '.ccp4.gz' and OEWriteGrid really does emit a
    // gzip stream. The header patch then reads compressed bytes as a header and
    // blames NSYMBT, which describes neither the cause nor the remedy. The
    // guard has to name compression instead.
    const MapFile source = read_map(DataPath("test_map.ccp4"));
    ScratchPath out(".ccp4.gz");
    try {
        write_map(out.Str(), *source.grid, source.symops);
        FAIL() << "expected GridError: a compressed destination cannot carry "
                  "the records this path patches in";
    } catch (const GridError& error) {
        const std::string message(error.what());
        EXPECT_NE(message.find("compress"), std::string::npos)
            << "the message does not name compression: " << message;
        EXPECT_EQ(message.find("NSYMBT"), std::string::npos)
            << "the message still blames NSYMBT: " << message;
    }

    std::ifstream probe(out.Str(), std::ios::binary);
    EXPECT_FALSE(probe.good()) << "a refused write left a file behind";
}

TEST(MapIoWriteTest, WritesADestinationWhoseStemEndsInACompressionSuffix) {
    // The other side of RefusesACompressedDestination. '.ccp4.gz' names a
    // compressed file and is refused; 'x.gz.ccp4' names a plain CCP4 file and
    // has to be written.
    //
    // OEWriteGrid reads the whole filename, not just the final suffix.
    // Measured on four temporary-shaped names carrying the same grid,
    // '.x.gz-<hex>.ccp4' came out gzipped while '.x-<hex>.ccp4',
    // '.x.GZ-<hex>.ccp4' and '.x.y-<hex>.ccp4' came out plain. A temporary
    // that embedded the destination's stem therefore came out gzipped for both
    // of the stems below, and PatchHeaderRecords then read the compressed
    // bytes as a CCP4 header and refused the write against a garbage NSYMBT.
    // Both stems, because 'x.gzz.ccp4' was refused for the same reason as
    // 'x.gz.ccp4': the component that provokes this does not have to be
    // exactly 'gz'.
    const MapFile source = read_map(DataPath("test_map.ccp4"));
    for (const char* extension : {".gz.ccp4", ".gzz.ccp4"}) {
        SCOPED_TRACE(extension);
        ScratchPath out(extension);
        write_map(out.Str(), *source.grid, source.symops);
        ASSERT_TRUE(std::filesystem::exists(out.Str()))
            << "the write returned without producing " << out.Str();
        ExpectRoundTrip(*source.grid, out.Str(), source.symops);
        EXPECT_EQ(CountTemporarySiblings(out.Str()), 0u)
            << "the write left its hidden temporary beside " << out.Str();
    }
}

TEST(MapIoWriteTest, RefusesAGridFormatOpenEyeWritesButThisPathCannotPatch) {
    // '.grd' is a format OpenEye writes -- OEGetGridFileType("grd") is GRD and
    // OEIsWriteableGrid accepts it -- but not one whose header holds the CCP4
    // records this path splices back in. Under the old writeability gate it got
    // as far as OEWriteGrid and PatchHeaderRecords before the verify refused
    // it, and the refusal blamed the map for not reading back rather than the
    // extension for selecting a format this writer cannot patch. Same shape as
    // RefusesACompressedDestination: refuse up front, name the real reason.
    const MapFile source = read_map(DataPath("test_map.ccp4"));
    ScratchPath out(".grd");
    try {
        write_map(out.Str(), *source.grid, source.symops);
        FAIL() << "expected GridError: a GRD destination cannot carry the CCP4 "
                  "records this path patches in";
    } catch (const GridError& error) {
        const std::string message(error.what());
        EXPECT_NE(message.find("CCP4"), std::string::npos)
            << "the message does not name the format family this writer "
               "produces: " << message;
        EXPECT_EQ(message.find("does not read back"), std::string::npos)
            << "the message still blames the map for not reading back: "
            << message;
    }

    std::ifstream probe(out.Str(), std::ios::binary);
    EXPECT_FALSE(probe.good()) << "a refused write left a file behind";
    EXPECT_EQ(CountTemporarySiblings(out.Str()), 0u)
        << "a refused write left its hidden temporary behind";
}

TEST(MapIoWriteTest, RefusesAnUnpatchableExtensionBeforeItReachesOEWriteGrid) {
    // What makes the extension gate worth having is its position: it runs
    // ahead of MakeTemporarySibling and OEWriteGrid, so a '.grd' destination
    // never reaches the writer and no GRD file gets CCP4-offset bytes spliced
    // into it. RefusesAGridFormatOpenEyeWritesButThisPathCannotPatch cannot
    // pin that -- with the old OEIsWriteableGrid gate in place its
    // destination-absence and sibling-count assertions both still pass,
    // because the later refusal is atomic too. A parent directory that does
    // not exist separates the two by message: the gate refuses on the
    // extension, while anything that gets past it fails at OEWriteGrid, which
    // has nowhere to open a file.
    const MapFile source = read_map(DataPath("test_map.ccp4"));
    const std::filesystem::path absent =
        std::filesystem::path(::testing::TempDir()) /
        ("maptitude_absent_" + std::to_string(::getpid()));
    ASSERT_FALSE(std::filesystem::exists(absent))
        << "the parent directory exists, so this case no longer separates the "
           "gate from the write: " << absent.string();
    const std::string out = (absent / "out.grd").string();

    try {
        write_map(out, *source.grid, source.symops);
        FAIL() << "expected GridError: '.grd' is not an extension OpenEye maps "
                  "to CCP4";
    } catch (const GridError& error) {
        const std::string message(error.what());
        EXPECT_NE(message.find("CCP4"), std::string::npos)
            << "the message does not name the format family this writer "
               "produces: " << message;
        EXPECT_EQ(message.find("OEWriteGrid failed"), std::string::npos)
            << "the extension reached OEWriteGrid, so the gate did not run "
               "ahead of it: " << message;
    }

    EXPECT_FALSE(std::filesystem::exists(absent))
        << "a refused write created the destination's parent directory";
}

TEST(MapIoWriteTest, RefusesASymopRecordWiderThanTheOnDiskField) {
    // The on-disk record is a fixed 80 bytes and the block builder pads to it
    // with resize(), which truncates just as readily. Past 80 the record lost
    // its tail, and the failure surfaced from the verify's own read_map as a
    // component-count error -- blaming the caller's text for a shape ParseAll
    // had already accepted.
    std::string wide("x");
    for (int i = 0; i < 38; ++i) {
        wide += "+0";
    }
    wide += ",y,z";
    ASSERT_EQ(wide.size(), 81u)
        << "this record is no longer one byte over the 80-byte field, so it "
           "tests nothing";
    ASSERT_EQ(SymOp::ParseAll(wide).size(), 1u)
        << "the parser rejects this record, so write_map would refuse it "
           "before reaching the width guard";

    const MapFile source = read_map(DataPath("test_map.ccp4"));
    ScratchPath out(".ccp4");
    try {
        write_map(out.Str(), *source.grid, wide);
        FAIL() << "expected SymOpError: the record does not fit the 80-byte "
                  "field";
    } catch (const SymOpError& error) {
        const std::string message(error.what());
        EXPECT_NE(message.find("81"), std::string::npos)
            << "the message does not name the record's length: " << message;
        EXPECT_EQ(message.find("3 components"), std::string::npos)
            << "the message still blames the component count: " << message;
    }

    std::ifstream probe(out.Str(), std::ios::binary);
    EXPECT_FALSE(probe.good()) << "a refused write left a file behind";
    EXPECT_EQ(CountTemporarySiblings(out.Str()), 0u)
        << "the SymOpError unwind left the hidden temporary behind";
}

TEST(MapIoWriteTest, RoundTripsAGridCarryingANaNVoxel) {
    // NaN compares unequal to itself, so a bytewise-faithful round trip of a
    // masked voxel used to read as a difference and the write was refused.
    MapFile source = read_map(DataPath("test_map.ccp4"));
    float* values = source.grid->GetValues();
    values[100] = std::numeric_limits<float>::quiet_NaN();
    ASSERT_TRUE(std::isnan(source.grid->GetValues()[100]));

    ScratchPath out(".ccp4");
    write_map(out.Str(), *source.grid, source.symops);
    ExpectRoundTrip(*source.grid, out.Str(), source.symops);

    // Without this the verify could be passing because the NaN never reached
    // the file, which is the opposite of what the fix is for.
    const MapFile back = read_map(out.Str());
    EXPECT_TRUE(std::isnan(back.grid->GetValues()[100]))
        << "the NaN did not survive the round trip";
}

TEST(MapIoWriteTest, WritesAnEmMapWhoseTwoPlacementRecordsAgree) {
    // The written header states node 0 twice, and the verify re-reads with the
    // default ORIGIN tiebreak, which discards NxSTART on exactly the maps whose
    // NxSTART can be wrong. This is the only shipped fixture with an off-lattice
    // node 0, so it is the only one where the two records can come apart.
    const MapFile source = read_map(AssetPath("390_emd_30342_A_z4.mrc"));
    ScratchPath out(".mrc");
    write_map(out.Str(), *source.grid, source.symops);

    const MapFile by_origin = read_map(out.Str(), OriginSource::ORIGIN_RECORD);
    const MapFile by_nxstart = read_map(out.Str(), OriginSource::NXSTART);
    EXPECT_EQ(compare_placement_records(*by_origin.grid, *by_nxstart.grid), "");

    // Without this the agreement above would also hold on a file that carried
    // no NxSTART at all, or whose two reads were identical for some unrelated
    // reason -- neither of which is what the check is meant to establish.
    const GridParams origin_gp = get_grid_params(*by_origin.grid);
    const GridParams nxstart_gp = get_grid_params(*by_nxstart.grid);
    const double divergence =
        std::fabs(origin_gp.x_origin - nxstart_gp.x_origin);
    EXPECT_GT(divergence, 0.0)
        << "the two reads place node 0 identically, so this file cannot "
           "distinguish a correct NxSTART from an absent one";
    EXPECT_LT(divergence, 0.5 * origin_gp.x_spacing)
        << "the two placement records are further apart than the rounding an "
           "integer node index forces";
}

TEST(MapIoWriteTest, AcceptsANonDyadicSpacingWhosePlacementLandsInsideTheSlack) {
    // An even dim with SetMid(0,0,0) puts node 0 at a half-integer number of
    // spacings whatever the spacing is, so the two placement records sit on the
    // half-interval bound. With 47 divisions the spacing is not a dyadic
    // rational, NCSTART * spacing is not exactly representable in float32, and
    // the reconstructed NxSTART placement lands 9.4e-7 relative beyond the half.
    // MAP_PLACEMENT_SLACK is what keeps this a correct write: with the slack at
    // zero, write_map refuses this grid.
    OESystem::OESkewGrid grid;
    ASSERT_TRUE(grid.SetDim(20u, 20u, 20u));
    ASSERT_TRUE(grid.SetUnitCell(20.0f, 20.0f, 20.0f, 90.0f, 90.0f, 90.0f,
                                 47u, 47u, 47u));
    ASSERT_TRUE(grid.SetMid(0.0f, 0.0f, 0.0f));
    ASSERT_TRUE(grid.SetSpaceGroup(1u));
    std::vector<float> values(grid.GetSize());
    for (std::size_t i = 0; i < values.size(); ++i) {
        values[i] = static_cast<float>(i % 17u);
    }
    ASSERT_TRUE(grid.SetValues(values.data(),
                               static_cast<unsigned int>(values.size())));

    ScratchPath out(".ccp4");
    write_map(out.Str(), grid);

    // The premise, stated rather than assumed: this pair really does sit above
    // the bare half-interval. Without it the acceptance above would pass for an
    // implementation with no slack at all, which is the gap this test closes.
    const GridParams by_origin =
        get_grid_params(*read_map(out.Str(), OriginSource::ORIGIN_RECORD).grid);
    const GridParams by_nxstart =
        get_grid_params(*read_map(out.Str(), OriginSource::NXSTART).grid);
    const double divergence =
        std::fabs(by_origin.x_origin - by_nxstart.x_origin);
    const double half = 0.5 * by_origin.x_spacing;
    EXPECT_GT(divergence, half)
        << "the divergence no longer exceeds half the node interval, so this "
           "grid no longer distinguishes MAP_PLACEMENT_SLACK from zero";
    EXPECT_LT(divergence, half * (1.0 + MAP_PLACEMENT_SLACK));
}

namespace {

/// A small grid the seam tests mutate one quantity at a time.
OESystem::OESkewGrid SeamGrid() {
    OESystem::OESkewGrid grid;
    grid.SetDim(4u, 5u, 6u);
    grid.SetUnitCell(4.0f, 5.0f, 6.0f, 90.0f, 90.0f, 90.0f, 8u, 10u, 12u);
    grid.SetMid(0.0f, 0.0f, 0.0f);
    grid.SetSpaceGroup(1u);
    std::vector<float> values(grid.GetSize(), 1.0f);
    grid.SetValues(values.data(), static_cast<unsigned int>(values.size()));
    return grid;
}

}  // namespace

TEST(MapIoVerifySeamTest, AgreesWithItself) {
    const OESystem::OESkewGrid grid = SeamGrid();
    EXPECT_EQ(compare_written_map(grid, grid), "");
}

TEST(MapIoVerifySeamTest, NamesADimensionDifference) {
    const OESystem::OESkewGrid expected = SeamGrid();
    OESystem::OESkewGrid actual = SeamGrid();
    ASSERT_TRUE(actual.SetDim(4u, 5u, 7u));
    const std::string message = compare_written_map(expected, actual);
    EXPECT_NE(message.find("dimension"), std::string::npos) << message;
}

TEST(MapIoVerifySeamTest, NamesACellDifference) {
    const OESystem::OESkewGrid expected = SeamGrid();
    OESystem::OESkewGrid actual = SeamGrid();
    ASSERT_TRUE(actual.SetUnitCell(4.5f, 5.0f, 6.0f, 90.0f, 90.0f, 90.0f,
                                   8u, 10u, 12u));
    const std::string message = compare_written_map(expected, actual);
    EXPECT_NE(message.find("cell"), std::string::npos) << message;
}

TEST(MapIoVerifySeamTest, NamesASpacingDifference) {
    // Same dim and same declared cell edges, different division count -- so
    // the spacing moves while dim and cell do not. A verify that skipped
    // spacing would pass this pair.
    const OESystem::OESkewGrid expected = SeamGrid();
    OESystem::OESkewGrid actual = SeamGrid();
    ASSERT_TRUE(actual.SetUnitCell(4.0f, 5.0f, 6.0f, 90.0f, 90.0f, 90.0f,
                                   16u, 20u, 24u));

    // The pair rests on SetUnitCell's division count moving the derived
    // spacing. If it stops doing that, the assertion below still passes for the
    // wrong reason, so state the premise rather than assume it.
    const GridParams expected_gp = get_grid_params(expected);
    const GridParams actual_gp = get_grid_params(actual);
    ASSERT_NE(expected_gp.x_spacing, actual_gp.x_spacing)
        << "the two grids have the same spacing, so this pair no longer "
           "discriminates a verify that skips spacing";
    ASSERT_EQ(expected_gp.x_dim, actual_gp.x_dim);

    const std::string message = compare_written_map(expected, actual);
    EXPECT_NE(message.find("spacing"), std::string::npos) << message;
}

TEST(MapIoVerifySeamTest, NamesANodeZeroDifference) {
    const OESystem::OESkewGrid expected = SeamGrid();
    OESystem::OESkewGrid actual = SeamGrid();
    ASSERT_TRUE(actual.SetMid(1.0f, 0.0f, 0.0f));
    const std::string message = compare_written_map(expected, actual);
    EXPECT_NE(message.find("node 0"), std::string::npos) << message;
}

TEST(MapIoVerifySeamTest, NamesAVoxelDifference) {
    const OESystem::OESkewGrid expected = SeamGrid();
    OESystem::OESkewGrid actual = SeamGrid();
    std::vector<float> values(actual.GetSize(), 1.0f);
    values[7] = 2.0f;
    ASSERT_TRUE(actual.SetValues(values.data(),
                                 static_cast<unsigned int>(values.size())));
    const std::string message = compare_written_map(expected, actual);
    EXPECT_NE(message.find("voxel"), std::string::npos) << message;
}

TEST(MapIoPlacementSeamTest, AgreesWhenBothRecordsPlaceNodeZeroAlike) {
    const OESystem::OESkewGrid grid = SeamGrid();
    EXPECT_EQ(compare_placement_records(grid, grid), "");
}

TEST(MapIoPlacementSeamTest, AcceptsADisagreementOfExactlyHalfANodeInterval) {
    // The bound is inclusive. Measured at spacing 0.5: a node 0 of 2.75 makes
    // OpenEye round away from zero to NCSTART 6, putting the NxSTART placement
    // exactly half a node interval from the ORIGIN one -- and that write is
    // correct. A >= comparison would refuse it.
    const OESystem::OESkewGrid by_origin = SeamGrid();
    OESystem::OESkewGrid by_nxstart = SeamGrid();
    ASSERT_TRUE(by_nxstart.SetMid(0.25f, 0.0f, 0.0f));

    const GridParams origin_gp = get_grid_params(by_origin);
    const GridParams nxstart_gp = get_grid_params(by_nxstart);
    ASSERT_NEAR(std::fabs(origin_gp.x_origin - nxstart_gp.x_origin),
                0.5 * origin_gp.x_spacing, 1e-6)
        << "the shift is no longer exactly half a node interval, so this pair "
           "no longer sits on the bound it is meant to pin";

    EXPECT_EQ(compare_placement_records(by_origin, by_nxstart), "");
}

TEST(MapIoPlacementSeamTest, NamesAnAxisWhoseRecordsPlaceNodeZeroApart) {
    const OESystem::OESkewGrid by_origin = SeamGrid();
    OESystem::OESkewGrid by_nxstart = SeamGrid();
    ASSERT_TRUE(by_nxstart.SetMid(1.0f, 0.0f, 0.0f));

    const GridParams origin_gp = get_grid_params(by_origin);
    const GridParams nxstart_gp = get_grid_params(by_nxstart);
    ASSERT_NEAR(std::fabs(origin_gp.x_origin - nxstart_gp.x_origin) /
                    origin_gp.x_spacing,
                2.0, 1e-6)
        << "the disagreement is no longer about two node intervals";

    const std::string message =
        compare_placement_records(by_origin, by_nxstart);
    EXPECT_NE(message.find("node 0 x"), std::string::npos) << message;
}

TEST(MapIoVerifySeamTest, TreatsTwoNaNVoxelsAsAgreeing) {
    OESystem::OESkewGrid expected = SeamGrid();
    OESystem::OESkewGrid actual = SeamGrid();
    std::vector<float> values(expected.GetSize(), 1.0f);
    values[7] = std::numeric_limits<float>::quiet_NaN();
    ASSERT_TRUE(expected.SetValues(values.data(),
                                   static_cast<unsigned int>(values.size())));
    ASSERT_TRUE(actual.SetValues(values.data(),
                                 static_cast<unsigned int>(values.size())));
    EXPECT_EQ(compare_written_map(expected, actual), "");
}

TEST(MapIoVerifySeamTest, StillNamesANaNAgainstAFiniteVoxel) {
    // The NaN-tolerant comparison must stay one-sided: a masked voxel that came
    // back as a number, or the reverse, is a real difference.
    const OESystem::OESkewGrid expected = SeamGrid();
    OESystem::OESkewGrid actual = SeamGrid();
    std::vector<float> values(actual.GetSize(), 1.0f);
    values[7] = std::numeric_limits<float>::quiet_NaN();
    ASSERT_TRUE(actual.SetValues(values.data(),
                                 static_cast<unsigned int>(values.size())));
    const std::string message = compare_written_map(expected, actual);
    EXPECT_NE(message.find("voxel"), std::string::npos) << message;
}
