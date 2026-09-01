#include "maptitude/MapIO.h"

#include "maptitude/Error.h"
#include "maptitude/Grid.h"
#include "maptitude/SymOp.h"

#include <oechem.h>
#include <oegrid.h>

#include <cstdint>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace Maptitude {

namespace {

/// A CCP4/MRC header is exactly 1024 bytes, followed by NSYMBT bytes of
/// symmetry records and then the payload.
constexpr std::size_t CCP4_HEADER_BYTES = 1024;

/// Symmetry records are fixed-width and space-padded, with no separator.
constexpr std::size_t SYMOP_RECORD_BYTES = 80;

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

}  // namespace Maptitude
