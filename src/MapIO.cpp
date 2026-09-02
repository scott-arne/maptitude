#include "maptitude/MapIO.h"

#include "maptitude/Error.h"
#include "maptitude/Grid.h"
#include "maptitude/SymOp.h"

#include <oechem.h>
#include <oegrid.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <random>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace Maptitude {

namespace {

/// A CCP4/MRC header is exactly 1024 bytes, followed by NSYMBT bytes of
/// symmetry records and then the payload.
constexpr std::size_t CCP4_HEADER_BYTES = 1024;

/// Symmetry records are fixed-width and space-padded, with no separator.
constexpr std::size_t SYMOP_RECORD_BYTES = 80;

/// Maximum number of symmetry records allowed. The largest crystallographic
/// space group has 192 general positions; this cap is twenty times that while
/// keeping the allocation bounded at 320 KB.
constexpr std::size_t MAX_SYMOP_RECORDS = 4096;

// Word indices are 1-based, as the CCP4 specification numbers them.
constexpr std::size_t WORD_NCSTART = 5;
constexpr std::size_t WORD_NSYMBT = 24;
constexpr std::size_t WORD_ORIGIN = 50;
constexpr std::size_t WORD_MACHST = 54;

/// The four records OEReadGrid discards or never exposes.
struct MapHeader {
    std::int32_t nxstart[3] = {0, 0, 0};
    std::int32_t nsymbt = 0;
    float origin[3] = {0.0f, 0.0f, 0.0f};
    bool byte_swapped = false;  ///< File byte order differs from this host's.
};

bool HostIsLittleEndian() {
    const std::uint32_t probe = 1u;
    unsigned char first = 0u;
    std::memcpy(&first, &probe, 1);
    return first == 1u;
}

std::uint32_t Swap32(const std::uint32_t value) {
    return ((value & 0x000000FFu) << 24) | ((value & 0x0000FF00u) << 8) |
           ((value & 0x00FF0000u) >> 8) | ((value & 0xFF000000u) >> 24);
}

std::uint32_t RawWord(const std::string& header, const std::size_t word) {
    std::uint32_t value = 0u;
    std::memcpy(&value, header.data() + (word - 1u) * 4u, 4);
    return value;
}

std::int32_t WordAsInt(const std::string& header, const std::size_t word,
                       const bool swapped) {
    const std::uint32_t raw = swapped ? Swap32(RawWord(header, word))
                                      : RawWord(header, word);
    std::int32_t value = 0;
    std::memcpy(&value, &raw, 4);
    return value;
}

float WordAsFloat(const std::string& header, const std::size_t word,
                  const bool swapped) {
    const std::uint32_t raw = swapped ? Swap32(RawWord(header, word))
                                      : RawWord(header, word);
    float value = 0.0f;
    std::memcpy(&value, &raw, 4);
    return value;
}

/// Read the fixed 1024-byte header, or throw naming what was wrong with it.
std::string ReadRawHeader(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw GridError("Cannot open map file '" + path + "'");
    }
    std::string header(CCP4_HEADER_BYTES, '\0');
    in.read(&header[0], static_cast<std::streamsize>(CCP4_HEADER_BYTES));
    if (in.gcount() != static_cast<std::streamsize>(CCP4_HEADER_BYTES)) {
        throw GridError("Map file '" + path + "' is shorter than a " +
                        std::to_string(CCP4_HEADER_BYTES) +
                        "-byte CCP4 header");
    }
    return header;
}

/// MACHST's first byte is 0x44 for a little-endian writer and 0x11 for a
/// big-endian one. Files written by tools that leave it zero are read in host
/// order, which is what the toolkit does with them too.
bool FileIsBigEndian(const std::string& header) {
    const unsigned char stamp =
        static_cast<unsigned char>(header[(WORD_MACHST - 1u) * 4u]);
    return stamp == 0x11u;
}

MapHeader ParseHeader(const std::string& raw) {
    MapHeader header;
    header.byte_swapped = FileIsBigEndian(raw) == HostIsLittleEndian();
    for (std::size_t axis = 0; axis < 3; ++axis) {
        header.nxstart[axis] =
            WordAsInt(raw, WORD_NCSTART + axis, header.byte_swapped);
        header.origin[axis] =
            WordAsFloat(raw, WORD_ORIGIN + axis, header.byte_swapped);
    }
    header.nsymbt = WordAsInt(raw, WORD_NSYMBT, header.byte_swapped);
    return header;
}

bool AnyNonzero(const std::int32_t (&values)[3]) {
    return values[0] != 0 || values[1] != 0 || values[2] != 0;
}

bool AnyNonzero(const float (&values)[3]) {
    return values[0] != 0.0f || values[1] != 0.0f || values[2] != 0.0f;
}

std::string Strip(const std::string& text) {
    const std::size_t first = text.find_first_not_of(" \t\r\n\f\v");
    if (first == std::string::npos) {
        return std::string();
    }
    const std::size_t last = text.find_last_not_of(" \t\r\n\f\v");
    return text.substr(first, last - first + 1u);
}

/// Turn the raw fixed-width block into the newline-delimited text SymOp::ParseAll
/// accepts. The parser rejects the raw block outright, so the join is required
/// rather than cosmetic.
std::string ReadSymopBlock(const std::string& path, const MapHeader& header) {
    if (header.nsymbt == 0) {
        return std::string();
    }
    if (header.nsymbt < 0 ||
        header.nsymbt % static_cast<std::int32_t>(SYMOP_RECORD_BYTES) != 0) {
        throw GridError("Map file '" + path + "' declares NSYMBT " +
                        std::to_string(header.nsymbt) +
                        ", which is not a non-negative multiple of " +
                        std::to_string(SYMOP_RECORD_BYTES));
    }

    const std::size_t record_count =
        static_cast<std::size_t>(header.nsymbt) / SYMOP_RECORD_BYTES;
    if (record_count > MAX_SYMOP_RECORDS) {
        const std::size_t cap_bytes = MAX_SYMOP_RECORDS * SYMOP_RECORD_BYTES;
        throw GridError("Map file '" + path + "' declares " +
                        std::to_string(record_count) +
                        " symmetry records (" +
                        std::to_string(header.nsymbt) + " bytes), exceeds cap of " +
                        std::to_string(MAX_SYMOP_RECORDS) + " records (" +
                        std::to_string(cap_bytes) + " bytes)");
    }

    const std::size_t block_bytes = static_cast<std::size_t>(header.nsymbt);
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw GridError("Cannot reopen map file '" + path +
                        "' to read its symmetry block");
    }

    // Check that the file is large enough before allocating block_bytes.
    in.seekg(0, std::ios::end);
    const std::streampos file_size_pos = in.tellg();
    if (!in || file_size_pos < 0) {
        throw GridError("Cannot determine size of map file '" + path + "'");
    }
    const std::size_t file_size = static_cast<std::size_t>(file_size_pos);
    const std::size_t required_size = CCP4_HEADER_BYTES + block_bytes;
    if (file_size < required_size) {
        throw GridError("Map file '" + path + "' declares NSYMBT " +
                        std::to_string(header.nsymbt) + " (" +
                        std::to_string(block_bytes) +
                        " bytes) but file size is " +
                        std::to_string(file_size) + " bytes, need at least " +
                        std::to_string(required_size) + " bytes");
    }

    in.seekg(static_cast<std::streamoff>(CCP4_HEADER_BYTES));
    std::string block(block_bytes, '\0');
    in.read(&block[0], static_cast<std::streamsize>(block_bytes));
    if (in.gcount() != static_cast<std::streamsize>(block_bytes)) {
        throw GridError("Map file '" + path + "' declares NSYMBT " +
                        std::to_string(header.nsymbt) +
                        " but ends before the symmetry block does");
    }

    std::string joined;
    for (std::size_t offset = 0; offset < block_bytes;
         offset += SYMOP_RECORD_BYTES) {
        const std::string record =
            Strip(block.substr(offset, SYMOP_RECORD_BYTES));
        if (record.empty()) {
            continue;
        }
        if (!joined.empty()) {
            joined += '\n';
        }
        joined += record;
    }
    return joined;
}

/// True when two quantities agree to within @p tol scaled by their magnitude.
///
/// The same shape as the predicate in Grid.cpp, which lives in an anonymous
/// namespace there and so cannot be linked from this translation unit.
bool NearlyEqualRelative(const double p, const double q, const double tol) {
    const double scale = std::max(1.0, std::max(std::fabs(p), std::fabs(q)));
    return std::fabs(p - q) <= tol * scale;
}

std::string Describe(const char* label, const double expected,
                     const double actual) {
    std::ostringstream out;
    out << label << " differs: wrote " << actual << ", expected " << expected;
    return out.str();
}

/// Deletes its file unless released, so every failure path between the write
/// and the rename leaves nothing behind.
class TemporaryFile {
public:
    explicit TemporaryFile(std::filesystem::path path)
        : path_(std::move(path)) {}
    ~TemporaryFile() {
        if (!released_) {
            std::error_code ignored;
            std::filesystem::remove(path_, ignored);
        }
    }
    TemporaryFile(const TemporaryFile&) = delete;
    TemporaryFile& operator=(const TemporaryFile&) = delete;

    const std::filesystem::path& Path() const { return path_; }
    void Release() { released_ = true; }

private:
    std::filesystem::path path_;
    bool released_ = false;
};

/// Whether @p path names a file OpenEye would compress on the way out.
///
/// OEIsWriteableGrid accepts ".ccp4.gz" and OEWriteGrid honours it, emitting a
/// real gzip stream. Everything this file does afterwards -- the NSYMBT and
/// ORIGIN patch and the raw-header verify -- reads the result as a plain CCP4
/// header, so it ends up blaming whatever the compressed bytes decode to
/// instead of naming the format it cannot patch.
///
/// Only ".gz" is listed, because it is the only compression suffix measured to
/// get past OEIsWriteableGrid; any other one is already refused there, with a
/// message about the extension rather than about compression.
bool IsCompressedPath(const std::string& path) {
    std::string extension = std::filesystem::path(path).extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](const unsigned char c) {
                       return static_cast<char>(std::tolower(c));
                   });
    return extension == ".gz";
}

/// A hidden sibling of @p dest that keeps its extension and does not yet exist.
///
/// Same directory, so the rename stays within one filesystem where POSIX
/// rename is atomic. Same extension, because OpenEye's format dispatch reads
/// the name it is handed: a ".tmp" or ".ccp4.tmp" temporary is the case that
/// aborts the process, and the check on the caller's path would not cover it.
///
/// The existence check is a filter, not a reservation. It closes the case that
/// actually happens -- a temporary left behind by a crashed or killed run,
/// which would otherwise be silently truncated and could collide again on the
/// next call if std::random_device degrades to a deterministic sequence. It
/// does not close a TOCTOU race against a concurrent writer; nothing short of
/// an atomic create would, and mkstemps -- the obvious way to get one -- costs
/// a permissions regression. Measured on this machine under umask 0022:
/// mkstemps creates its file at mode 0600, rename carries that 0600 onto the
/// destination, and an ordinary create lands at 0644. So every map written
/// through mkstemps would come out readable only by the user who wrote it.
/// Two concurrent write_map calls to the same destination are already
/// unsafe at the rename regardless of how the temporary is named, so the
/// reservation would buy nothing the caller can rely on and would cost a
/// visible permissions regression.
std::filesystem::path MakeTemporarySibling(const std::filesystem::path& dest) {
    constexpr int MAX_ATTEMPTS = 8;
    std::random_device entropy;
    for (int attempt = 0; attempt < MAX_ATTEMPTS; ++attempt) {
        std::ostringstream name;
        name << "." << dest.stem().string() << "-" << std::hex << entropy()
             << dest.extension().string();
        std::filesystem::path candidate = dest.parent_path() / name.str();
        std::error_code ignored;
        if (!std::filesystem::exists(candidate, ignored)) {
            return candidate;
        }
    }
    throw GridError("Cannot write '" + dest.string() +
                    "': no free temporary name beside it after " +
                    std::to_string(MAX_ATTEMPTS) + " attempts");
}

/// Normalize symop text to the one-triplet-per-line form read_map returns.
std::string CanonicalSymops(const std::string& symops) {
    std::string joined;
    std::size_t start = 0;
    while (start <= symops.size() && !symops.empty()) {
        const std::size_t end = symops.find('\n', start);
        const std::string record = Strip(
            end == std::string::npos ? symops.substr(start)
                                     : symops.substr(start, end - start));
        if (!record.empty()) {
            if (!joined.empty()) {
                joined += '\n';
            }
            joined += record;
        }
        if (end == std::string::npos) {
            break;
        }
        start = end + 1;
    }
    return joined;
}

/// The fixed-width, space-padded, separator-free on-disk form.
std::string SymopBlockBytes(const std::string& canonical) {
    std::string block;
    if (canonical.empty()) {
        return block;
    }
    std::size_t start = 0;
    while (start <= canonical.size()) {
        const std::size_t end = canonical.find('\n', start);
        const std::string record =
            end == std::string::npos ? canonical.substr(start)
                                     : canonical.substr(start, end - start);
        // resize() truncates as readily as it pads, and the tail it drops is
        // usually a whole component. The truncated block then fails the verify's
        // own read_map with a component-count error against text ParseAll had
        // already accepted, so the caller is told their symop is malformed when
        // the real fault is that it does not fit the format's fixed field.
        // Records that came from read_map cannot reach this: it strips them out
        // of 80-byte fields to begin with. Hand-built symops can.
        if (record.size() > SYMOP_RECORD_BYTES) {
            throw SymOpError(
                "Symmetry record is " + std::to_string(record.size()) +
                " characters, over the " + std::to_string(SYMOP_RECORD_BYTES) +
                "-character CCP4 record field: '" + record + "'");
        }
        std::string padded = record;
        padded.resize(SYMOP_RECORD_BYTES, ' ');
        block += padded;
        if (end == std::string::npos) {
            break;
        }
        start = end + 1;
    }
    return block;
}

/// Read the @p count bytes of symmetry records that follow the 1024-byte
/// header.
///
/// ReadRawHeader stops at the header, which is all the read path needs; the
/// verify step needs the block itself, so it gets its own reader rather than
/// widening that one and making every read_map call carry the extra bytes.
std::string ReadSymopBlock(const std::filesystem::path& file,
                           const std::size_t count) {
    std::string block(count, '\0');
    if (count == 0u) {
        return block;
    }
    std::ifstream in(file, std::ios::binary);
    if (!in) {
        throw GridError("Cannot reopen '" + file.string() +
                        "' to verify its symmetry block");
    }
    in.seekg(static_cast<std::streamoff>(CCP4_HEADER_BYTES));
    in.read(&block[0], static_cast<std::streamsize>(count));
    if (in.gcount() != static_cast<std::streamsize>(count)) {
        throw GridError("'" + file.string() + "' declares " +
                        std::to_string(count) +
                        " bytes of symmetry records but holds fewer");
    }
    return block;
}

/// Write one 4-byte header word, counting words from 1 as the format does.
///
/// The mirror of WordAsInt's read side, and the only thing in this file that
/// writes inside the 1024-byte header -- PatchHeaderRecords also splices the
/// symop block, but that starts past it. It takes the word already in file byte
/// order; the caller swaps.
void SetRawWord(std::string& bytes, const std::size_t word,
                const std::uint32_t raw) {
    std::memcpy(&bytes[(word - 1u) * 4u], &raw, 4);
}

/// Splice NSYMBT, the symop block and ORIGIN into a file OEWriteGrid produced.
///
/// NX/NY/NZ is deliberately not patched: section 2.4 measures the dim - 1
/// patch an earlier draft applied here as a no-op on whole-cell grids and a
/// silent corruption on a sub-box, where it preserves dim, cell and voxel
/// count while moving the spacing and the grid's position.
void PatchHeaderRecords(const std::filesystem::path& file,
                        const std::string& canonical_symops,
                        const OESystem::OESkewGrid& grid) {
    std::string bytes;
    {
        std::ifstream in(file, std::ios::binary);
        if (!in) {
            throw GridError("Cannot reopen '" + file.string() +
                            "' to restore its header records");
        }
        bytes.assign(std::istreambuf_iterator<char>(in),
                     std::istreambuf_iterator<char>());
    }
    if (bytes.size() < CCP4_HEADER_BYTES) {
        throw GridError("OEWriteGrid produced '" + file.string() +
                        "' with a header shorter than " +
                        std::to_string(CCP4_HEADER_BYTES) + " bytes");
    }

    // Match whatever byte order OEWriteGrid emitted rather than assuming one.
    const bool swapped = FileIsBigEndian(bytes) == HostIsLittleEndian();
    auto put_int = [&bytes, swapped](const std::size_t word,
                                     const std::int32_t value) {
        std::uint32_t raw = 0u;
        std::memcpy(&raw, &value, 4);
        SetRawWord(bytes, word, swapped ? Swap32(raw) : raw);
    };
    auto put_float = [&bytes, swapped](const std::size_t word,
                                       const float value) {
        std::uint32_t raw = 0u;
        std::memcpy(&raw, &value, 4);
        SetRawWord(bytes, word, swapped ? Swap32(raw) : raw);
    };

    const std::int32_t existing_nsymbt = WordAsInt(bytes, WORD_NSYMBT, swapped);
    if (existing_nsymbt < 0 ||
        static_cast<std::size_t>(existing_nsymbt) >
            bytes.size() - CCP4_HEADER_BYTES) {
        throw GridError("OEWriteGrid produced '" + file.string() +
                        "' with an unusable NSYMBT of " +
                        std::to_string(existing_nsymbt));
    }

    const std::string block = SymopBlockBytes(canonical_symops);
    bytes.replace(CCP4_HEADER_BYTES, static_cast<std::size_t>(existing_nsymbt),
                  block);
    put_int(WORD_NSYMBT, static_cast<std::int32_t>(block.size()));

    float node_x = 0.0f, node_y = 0.0f, node_z = 0.0f;
    if (!grid.ElementToSpatialCoord(0u, node_x, node_y, node_z)) {
        throw GridError("Cannot read the position of element 0 to write the "
                        "ORIGIN record of '" + file.string() + "'");
    }
    if (node_x != 0.0f || node_y != 0.0f || node_z != 0.0f) {
        put_float(WORD_ORIGIN + 0u, node_x);
        put_float(WORD_ORIGIN + 1u, node_y);
        put_float(WORD_ORIGIN + 2u, node_z);
    }

    std::ofstream out(file, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw GridError("Cannot rewrite '" + file.string() +
                        "' with its restored header records");
    }
    out.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    out.close();
    if (!out) {
        throw GridError("Failed while rewriting '" + file.string() + "'");
    }
}

}  // namespace

MapFile::MapFile() = default;
MapFile::MapFile(MapFile&& other) noexcept = default;
MapFile& MapFile::operator=(MapFile&& other) noexcept = default;
MapFile::~MapFile() = default;

MapFile read_map(const std::string& path, const OriginSource tiebreak) {
    MapFile result;
    result.grid.reset(new OESystem::OESkewGrid());
    if (!OESystem::OEReadGrid(path, *result.grid)) {
        throw GridError("OEReadGrid could not read map file '" + path + "'");
    }

    const std::string raw = ReadRawHeader(path);
    const MapHeader header = ParseHeader(raw);

    // Spacing, node count and axis order are left exactly as the reader
    // produced them: section 2.1 measures that it already applies NxSTART
    // through the MAPC/MAPR/MAPS permutation. The only record it discards is
    // ORIGIN, so that is the only correction made here.
    const bool origin_set = AnyNonzero(header.origin);
    const bool nxstart_set = AnyNonzero(header.nxstart);
    const bool apply_origin =
        origin_set && (!nxstart_set || tiebreak == OriginSource::ORIGIN_RECORD);
    if (apply_origin) {
        float node_x = 0.0f, node_y = 0.0f, node_z = 0.0f;
        if (!result.grid->ElementToSpatialCoord(0u, node_x, node_y, node_z)) {
            throw GridError("Cannot read the position of element 0 of '" +
                            path + "' to apply its ORIGIN record");
        }
        float mid_x = 0.0f, mid_y = 0.0f, mid_z = 0.0f;
        if (!result.grid->GetMid(mid_x, mid_y, mid_z)) {
            throw GridError("Cannot read the midpoint of '" + path +
                            "' to apply its ORIGIN record");
        }
        if (!result.grid->SetMid(mid_x + (header.origin[0] - node_x),
                                 mid_y + (header.origin[1] - node_y),
                                 mid_z + (header.origin[2] - node_z))) {
            throw GridError("Cannot move '" + path +
                            "' onto the origin its header declares");
        }
    }

    result.symops = ReadSymopBlock(path, header);
    // Validation by the parser the library already has, rather than a second
    // private notion of a valid symop. An empty block parses to zero operators.
    SymOp::ParseAll(result.symops);
    return result;
}

std::string compare_written_map(const OESystem::OESkewGrid& expected,
                                const OESystem::OESkewGrid& actual) {
    // Node counts first: a dim mismatch makes every later comparison
    // meaningless, and it is the one the sub-box case trips.
    if (expected.GetXDim() != actual.GetXDim() ||
        expected.GetYDim() != actual.GetYDim() ||
        expected.GetZDim() != actual.GetZDim()) {
        std::ostringstream out;
        out << "dimensions differ: wrote " << actual.GetXDim() << "x"
            << actual.GetYDim() << "x" << actual.GetZDim() << ", expected "
            << expected.GetXDim() << "x" << expected.GetYDim() << "x"
            << expected.GetZDim();
        return out.str();
    }
    if (expected.GetSize() != actual.GetSize()) {
        return Describe("voxel count", expected.GetSize(), actual.GetSize());
    }

    if (expected.HasUnitCell() != actual.HasUnitCell()) {
        return std::string("unit cell presence differs: wrote ") +
               (actual.HasUnitCell() ? "one" : "none") + ", expected " +
               (expected.HasUnitCell() ? "one" : "none");
    }
    if (expected.HasUnitCell()) {
        const UnitCellParams want = get_unit_cell(expected);
        const UnitCellParams got = get_unit_cell(actual);
        const struct {
            const char* label;
            double want;
            double got;
        } cell_terms[] = {
            {"cell edge a", want.a, got.a},
            {"cell edge b", want.b, got.b},
            {"cell edge c", want.c, got.c},
            {"cell angle alpha", want.alpha, got.alpha},
            {"cell angle beta", want.beta, got.beta},
            {"cell angle gamma", want.gamma, got.gamma},
        };
        for (const auto& term : cell_terms) {
            if (!NearlyEqualRelative(term.want, term.got, MAP_VERIFY_TOL)) {
                return Describe(term.label, term.want, term.got);
            }
        }
    }

    const GridParams want_gp = get_grid_params(expected);
    const GridParams got_gp = get_grid_params(actual);
    const struct {
        const char* label;
        double want;
        double got;
    } geometry_terms[] = {
        {"spacing along x", want_gp.x_spacing, got_gp.x_spacing},
        {"spacing along y", want_gp.y_spacing, got_gp.y_spacing},
        {"spacing along z", want_gp.z_spacing, got_gp.z_spacing},
        {"node 0 x", want_gp.x_origin, got_gp.x_origin},
        {"node 0 y", want_gp.y_origin, got_gp.y_origin},
        {"node 0 z", want_gp.z_origin, got_gp.z_origin},
    };
    for (const auto& term : geometry_terms) {
        if (!NearlyEqualRelative(term.want, term.got, MAP_VERIFY_TOL)) {
            return Describe(term.label, term.want, term.got);
        }
    }

    const float* want_values = expected.GetValues();
    const float* got_values = actual.GetValues();
    for (unsigned int i = 0; i < expected.GetSize(); ++i) {
        // NaN compares unequal to itself, so the exact comparison below reads a
        // faithfully round-tripped masked voxel as a difference. Two NaNs are
        // the same value for this purpose. The tolerance is one-sided: a NaN
        // against a number, either way round, is still a difference.
        if (std::isnan(want_values[i]) && std::isnan(got_values[i])) {
            continue;
        }
        if (want_values[i] != got_values[i]) {
            std::ostringstream out;
            out << "voxel " << i << " differs: wrote " << got_values[i]
                << ", expected " << want_values[i];
            return out.str();
        }
    }
    return std::string();
}

std::string compare_placement_records(const OESystem::OESkewGrid& by_origin,
                                      const OESystem::OESkewGrid& by_nxstart) {
    const GridParams origin_gp = get_grid_params(by_origin);
    const GridParams nxstart_gp = get_grid_params(by_nxstart);
    const struct {
        const char* label;
        double spacing;
        double from_origin;
        double from_nxstart;
    } axes[] = {
        {"node 0 x", origin_gp.x_spacing, origin_gp.x_origin,
         nxstart_gp.x_origin},
        {"node 0 y", origin_gp.y_spacing, origin_gp.y_origin,
         nxstart_gp.y_origin},
        {"node 0 z", origin_gp.z_spacing, origin_gp.z_origin,
         nxstart_gp.z_origin},
    };
    for (const auto& axis : axes) {
        const double interval = std::fabs(axis.spacing);
        const double divergence = std::fabs(axis.from_origin -
                                            axis.from_nxstart);
        if (divergence > 0.5 * interval * (1.0 + MAP_PLACEMENT_SLACK)) {
            std::ostringstream out;
            out << axis.label << " differs between the header's two placement "
                << "records by " << divergence << " A, more than half the "
                << interval << " A node interval: ORIGIN puts it at "
                << axis.from_origin << ", NxSTART at " << axis.from_nxstart;
            return out.str();
        }
    }
    return std::string();
}

void write_map(const std::string& path, const OESystem::OESkewGrid& grid,
               const std::string& symops) {
    // Step 1: validate both arguments before touching the filesystem, so a bad
    // call fails without leaving a partial file.
    SymOp::ParseAll(symops);
    if (IsCompressedPath(path)) {
        throw GridError("Cannot write '" + path +
                        "': this names a compressed file, and the CCP4 header "
                        "records this writer restores after OEWriteGrid cannot "
                        "be spliced into a compressed stream");
    }
    if (!OESystem::OEIsWriteableGrid(path)) {
        // Section 2.4 measures that OEWriteGrid on an extension OpenEye does
        // not recognize prints a fatal error and exits rather than returning
        // false. This check is the only thing between the caller and a killed
        // process; once OEWriteGrid has the path there is nothing to catch.
        throw GridError("Cannot write '" + path +
                        "': OpenEye writes no grid format with this "
                        "file extension");
    }

    // Step 2: copy, and default the space group. Section 2.3 shows an unset
    // space group is what makes OEWriteGrid double the cell and regrid.
    OESystem::OESkewGrid out(grid);
    if (!out.HasSpaceGroup() && !out.SetSpaceGroup(1u)) {
        throw GridError("Cannot default the space group to P1 before "
                        "writing '" + path + "'");
    }

    const std::filesystem::path dest(path);
    TemporaryFile temporary(MakeTemporarySibling(dest));

    // Step 3.
    if (!OESystem::OEWriteGrid(temporary.Path().string(), out)) {
        throw GridError("OEWriteGrid failed while writing '" + path + "'");
    }

    // Step 4.
    const std::string canonical = CanonicalSymops(symops);
    PatchHeaderRecords(temporary.Path(), canonical, out);

    // Step 5: verify against a read_map re-read. The round trip is closed
    // under read_map, not under OEReadGrid: the bare reader never consults
    // ORIGIN, so a map written with a nonzero one lands wrong by construction.
    const MapFile back = read_map(temporary.Path().string());
    const std::string difference = compare_written_map(out, *back.grid);
    if (!difference.empty()) {
        throw GridError("Refusing to write '" + path +
                        "': the map read back from disk " + difference);
    }
    if (back.symops != canonical) {
        throw GridError("Refusing to write '" + path +
                        "': the symmetry block read back from disk does not "
                        "match the one written");
    }

    // Both remaining checks read the raw bytes, because read_map normalizes
    // away exactly what they are checking. Read the header once, and
    // unconditionally: gating it on a nonzero node 0 would gate the NSYMBT
    // check with it. The symop check below opens the file a second time for the
    // block itself, which starts past what ReadRawHeader returns.
    const std::string raw = ReadRawHeader(temporary.Path().string());
    const MapHeader written = ParseHeader(raw);

    // back.symops above compares the parsed text. This compares the bytes:
    // NSYMBT and the fixed-width record shape are what a consumer that is not
    // read_map will read, and the text comparison passes on a block whose
    // padding or record count is wrong but which strips to the same triplets.
    const std::string expected_block = SymopBlockBytes(canonical);
    if (static_cast<std::size_t>(written.nsymbt) != expected_block.size()) {
        throw GridError("Refusing to write '" + path + "': NSYMBT is " +
                        std::to_string(written.nsymbt) + " but the symmetry "
                        "block is " + std::to_string(expected_block.size()) +
                        " bytes");
    }
    if (ReadSymopBlock(temporary.Path(), expected_block.size()) !=
        expected_block) {
        throw GridError("Refusing to write '" + path +
                        "': the symmetry records on disk are not the "
                        "80-byte space-padded block that was written");
    }

    // The ORIGIN record is checked separately from the placement above,
    // because a zero ORIGIN whose NxSTART happens to place node 0 correctly
    // would pass the placement comparison while losing the record.
    float node_x = 0.0f, node_y = 0.0f, node_z = 0.0f;
    if (!out.ElementToSpatialCoord(0u, node_x, node_y, node_z)) {
        // Checked rather than ignored: on failure the three stay 0.0f, which is
        // precisely the value that skips the comparison below. PatchHeaderRecords
        // already made this call succeed on this grid, so reaching here means
        // something changed underneath us.
        throw GridError("Refusing to write '" + path +
                        "': cannot read the position of element 0 to verify "
                        "the ORIGIN record");
    }
    if (node_x != 0.0f || node_y != 0.0f || node_z != 0.0f) {
        if (!NearlyEqualRelative(node_x, written.origin[0], MAP_VERIFY_TOL) ||
            !NearlyEqualRelative(node_y, written.origin[1], MAP_VERIFY_TOL) ||
            !NearlyEqualRelative(node_z, written.origin[2], MAP_VERIFY_TOL)) {
            throw GridError("Refusing to write '" + path +
                            "': the ORIGIN record read back from disk does "
                            "not match the grid's node 0");
        }
    }

    // The verify above re-read with the default ORIGIN tiebreak, which discards
    // NxSTART whenever ORIGIN is nonzero -- so nothing so far has looked at
    // NxSTART on precisely the maps where it can be wrong. read_map under
    // OriginSource::NXSTART is public API, and it reads that record.
    //
    // This check does not fire on anything measured: OpenEye derives NCSTART by
    // rounding ORIGIN/spacing to the nearest node, which bounds the two
    // placements half a node interval apart, and the shipped fixtures stay well
    // inside that. It is here to hold the bound as a guarantee rather than as
    // an observation about one toolkit version. The failures it would catch are
    // a toolkit that truncates where it rounds, an axis permutation pairing one
    // axis's NCSTART with another's spacing, and -- the case worth the cost --
    // a toolkit that stops emitting NCSTART, which would misplace a written EM
    // map by its whole origin offset for every NXSTART reader with nothing else
    // noticing.
    const MapFile by_nxstart =
        read_map(temporary.Path().string(), OriginSource::NXSTART);
    const std::string placement =
        compare_placement_records(*back.grid, *by_nxstart.grid);
    if (!placement.empty()) {
        throw GridError("Refusing to write '" + path + "': " + placement);
    }

    // Step 6: the only step that modifies the destination, and it runs only
    // after step 5 has passed.
    std::error_code rename_error;
    std::filesystem::rename(temporary.Path(), dest, rename_error);
    if (rename_error) {
        throw GridError("Verified map could not be renamed onto '" + path +
                        "': " + rename_error.message());
    }
    temporary.Release();
}

}  // namespace Maptitude
