#include "maptitude/MapIO.h"

#include "maptitude/Error.h"
#include "maptitude/Grid.h"
#include "maptitude/SymOp.h"

#include <oechem.h>
#include <oegrid.h>

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <locale>
#include <random>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

// The exclusive create in MakeTemporarySibling takes no platform guard:
// std::fopen's C11 "x" mode is standard C. getpid does take one -- MSVC has no
// <unistd.h> and declares _getpid() in <process.h> -- and ProcessId below is
// where it is settled. Both spellings are reached, because CMakeLists.txt
// builds this file into libmaptitude unconditionally and the release workflow
// builds a Windows wheel with MSVC.
#ifdef _WIN32
#include <process.h>
#else
#include <unistd.h>
#endif

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
constexpr std::size_t WORD_LSKFLG = 25;
constexpr std::size_t WORD_SKWTRN = 35;
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
///
/// Record boundaries are decided here and nowhere else: a record still carrying
/// a ';' or a newline once its padding is stripped is refused, for the reason
/// given at that check.
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
        // A separator that survives Strip sits inside one fixed-width record,
        // and is refused here rather than split. The SymOp::ParseAll call
        // read_map makes on the joined text cannot see it: that parser splits
        // on ';' as well as newline, so a separator inside a record reads to it
        // as one between records. Measured on tests/assets/mapq/1d26_2fofc.ccp4
        // with record 0 respliced as "x,y,z;-x,-y,z" and NSYMBT left at 640,
        // read_map returned that text verbatim as one of eight operators and
        // ParseAll passed it; the pristine file returns one operator per record
        // over those same eight records, so what the spliced read returned is
        // the splice's doing.
        //
        // Splitting the record here instead of refusing it would not close
        // that: CanonicalSymops splits on ';' too, so that spliced file round
        // trips as nine records whichever place does the splitting -- NSYMBT
        // 640 in, 720 out. Only refusing leaves that file's record count alone.
        // Refusing is also the read half of the position CanonicalSymops states
        // for the write half: a consumer reading the fixed-width records
        // literally gets one malformed operator out of a record holding two
        // triplets, not two well-formed ones.
        const std::size_t separator = record.find_first_of("\n;");
        if (separator != std::string::npos) {
            throw SymOpError(
                "Map file '" + path + "' holds " +
                (record[separator] == ';' ? "a ';'" : "a newline") +
                " inside symmetry record " +
                std::to_string(offset / SYMOP_RECORD_BYTES) + " of " +
                std::to_string(record_count) + " (counting from 0): a CCP4 "
                "symmetry record is a fixed " +
                std::to_string(SYMOP_RECORD_BYTES) +
                "-byte field carrying one operator, so a separator inside one "
                "separates nothing on disk");
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

/// Refuse a path carrying an embedded NUL, before anything opens it.
///
/// std::filesystem reads the whole std::string, while the calls that open,
/// write and rename the file take it as a NUL-terminated name and stop at the
/// first NUL. So "victim.dat\0.ccp4" reaches write_map's extension gate as a
/// ".ccp4" name and is then written through as "victim.dat": the name the gate
/// admitted and the name the rename published are different names, so passing
/// the gate says nothing about the file that is replaced. Measured before this
/// guard, a write to that name returned successfully having replaced an
/// existing victim.dat with a 38068-byte map. read_map has no gate to escalate
/// past, but the same split reaches it: measured before this guard, reading a
/// real map's path with "\0.ccp4" appended returned that map, under a name the
/// caller never gave.
///
/// The message names only the prefix: one interpolating the raw path would
/// carry the NUL too, and whatever printed it would stop there.
void RejectEmbeddedNul(const std::string& path, const std::string& verb) {
    const std::size_t nul = path.find('\0');
    if (nul == std::string::npos) {
        return;
    }
    throw GridError("Cannot " + verb + " '" + path.substr(0, nul) +
                    "...': the path contains an embedded NUL, so the name "
                    "checked here and the name the operating system would open "
                    "are not the same name");
}

/// Whether @p path names a file OpenEye would compress on the way out.
///
/// OEIsWriteableGrid accepts ".ccp4.gz" and OEWriteGrid honours it, emitting a
/// real gzip stream. Everything this file does afterwards -- the NSYMBT and
/// ORIGIN patch and the raw-header verify -- reads the result as a plain CCP4
/// header, so it ends up blaming whatever the compressed bytes decode to
/// instead of naming the format it cannot patch.
///
/// Only ".gz" is listed. Six other compression suffixes were measured against
/// OEIsWriteableGrid -- ".bz2", ".Z", ".zst", ".gzip" and ".xz" behind a
/// ".ccp4" stem, and a bare ".bz2" -- and every one is refused there, with a
/// message about the extension rather than about compression. Of the suffixes
/// measured, ".gz" is the only one that gets past it, in any letter case and
/// behind any stem: ".ccp4.gz", ".map.gz", ".mrc.gz" and ".ccp4.GZ" are all
/// accepted there and all caught here.
bool IsCompressedPath(const std::string& path) {
    std::string extension = std::filesystem::path(path).extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](const unsigned char c) {
                       return static_cast<char>(std::tolower(c));
                   });
    return extension == ".gz";
}

/// Re-word a nested read's message to name the destination, not the temporary.
///
/// The verify reads the file this writer just produced, which lives under a
/// hidden sibling name the caller never asked for and cannot look at. Every
/// GridError the read path raises interpolates its path argument unmodified, so
/// replacing that exact substring keeps the diagnosis and drops the leak.
std::string MessageAgainstDestination(const std::string& message,
                                      const std::string& temporary,
                                      const std::string& destination) {
    std::string out = message;
    for (std::size_t at = out.find(temporary); at != std::string::npos;
         at = out.find(temporary, at + destination.size())) {
        out.replace(at, temporary.size(), destination);
    }
    return out;
}

/// This process's id, under the spelling its platform gives the call.
///
/// The temporary's name carries the pid, and this file is compiled for the
/// MSVC wheel as well as the POSIX one: MSVC declares _getpid() in
/// <process.h>, POSIX declares getpid() in <unistd.h>. Both are the caller's
/// process, and neither has a failure return.
long ProcessId() {
#ifdef _WIN32
    return static_cast<long>(::_getpid());
#else
    return static_cast<long>(::getpid());
#endif
}

/// A hidden, empty sibling of @p dest that keeps its extension, created here
/// and owned by the TemporaryFile returned.
///
/// Same directory, so the rename stays within one filesystem where POSIX
/// rename is atomic. Same extension, because OpenEye's format dispatch reads
/// the name it is handed: a ".tmp" or ".ccp4.tmp" temporary is the case that
/// aborts the process, and the check on the caller's path would not cover it.
///
/// The basename is not built from the destination's stem, because that
/// dispatch reads the whole name and not only the final suffix. Measured with
/// bare OEWriteGrid over four temporary-shaped names carrying one grid:
/// ".x.gz-<hex>.ccp4" came out gzipped, while ".x-<hex>.ccp4",
/// ".x.GZ-<hex>.ccp4" and ".x.y-<hex>.ccp4" came out plain. So a destination
/// named "x.gz.ccp4" -- a plain CCP4 file by its own extension, and one the
/// compression guard on the caller's path passes -- got a gzipped temporary
/// out of the old scheme, and PatchHeaderRecords then read those compressed
/// bytes as a CCP4 header and refused the write against a garbage NSYMBT.
/// ".maptitude-<hex>.ccp4" was measured plain for that same destination, and
/// the file it produced read back through read_map. The pid this name now also
/// carries is decimal digits, so it adds no dot component to the shape those
/// four measured; the ".gz.ccp4" and ".gzz.ccp4" write cases exercise the name
/// that results.
///
/// The classic locale on that stream is what holds the previous sentence true.
/// num_put groups integers through the locale's numpunct, so a stream left on
/// the global locale follows whatever separator the process installed:
/// measured with this name built verbatim, en_US.UTF-8 gives
/// ".maptitude-5,656-de,adb,eef.ccp4" and de_DE.UTF-8 gives
/// ".maptitude-5.656-de.adb.eef.ccp4" -- three dot components in a shape
/// measured without any. The test tree's leak counter builds its prefix with
/// std::to_string, which never groups, so it would stop matching the temporary
/// this function made while still matching the decoys those cases plant. Every
/// assertion it backs would then pass without having seen the file it is there
/// to catch -- including the positive control those cases run to show the
/// counter can see something, which counts only decoys and so keeps reporting
/// success.
///
/// The name is taken, not merely chosen. The C11 "x" mode claims a name in the
/// same step that tests it, and the pid gives each process its own space of
/// names to draw from. An exists() call ahead of the create cannot do that --
/// it is a filter, and a second writer can slip between the test and the use.
/// A name already held comes back as EEXIST, measured here and what POSIX
/// specifies for that mode, and the loop draws another; the GridError at the
/// bottom is what a run of MAX_ATTEMPTS collisions raises. Any other errno is
/// not a collision, so it is reported as itself instead of being retried into
/// that message.
///
/// A symlink sitting at the candidate name does not redirect the create. POSIX
/// specifies that O_CREAT|O_EXCL fails on an existing path even where that path
/// is a symbolic link, whatever the link resolves to, so the attempt comes back
/// EEXIST and the loop draws a different name. Measured under both shapes: a
/// dangling link and a link resolving onto an existing file are each refused
/// with errno 17, and the file the second pointed at is not truncated.
///
/// A temporary left behind by a crashed or killed run is covered by the same
/// EEXIST -- skipped rather than truncated -- which is what the old existence
/// check was there for. Pids are reused, so the name space a leftover sits in
/// does not retire with the process that made it; the exclusive create, not the
/// pid, is what makes that harmless.
///
/// Under POSIX, std::fopen creates at 0666 before the umask, so the process
/// umask narrows the temporary as it narrows an ordinary create. ISO C fixes no
/// permissions for fopen, and Windows has no 0666 to create at: the permission
/// bits it documents for a created file are _S_IREAD and _S_IWRITE. That is why
/// this reaches for the "x" mode rather than mkstemps, the obvious way to get
/// an atomic create: mkstemps moves the permissions rather than fixing them.
/// Neither scheme preserves the destination's mode, because the rename replaces
/// its inode either way. Measured on this machine under umask 0022: this create
/// and an ordinary one both land at 0644, mkstemps at 0600, and the rename
/// carries whichever one onto the destination. So mkstemps would publish every
/// map readable only by its writer, while the scheme used here silently widens
/// a destination the caller had narrowed -- 0600 back to 0644 on rewrite.
/// Trading a visible regression on every write for an invisible one on rewrites
/// of a restricted file is the choice made, not a permissions-preserving option
/// that mkstemps lacks; the caller-facing consequences are stated on write_map.
///
/// OEWriteGrid is therefore handed a name that already exists. Measured against
/// this reservation on tests/data/test_map.ccp4 under umask 0022: it returned
/// true, filled the empty file, left the mode at the 0644 the create had
/// produced, and the result read back through read_map as the grid written.
/// The destination that rename published carried the same 0644 an unreserved
/// OEWriteGrid create produced at a destination beside it.
///
/// Two things the reservation does not close.
///
/// The first is two write_map calls naming one destination. Each now draws a
/// temporary the other's exclusive create is refused at, so the exposure
/// recorded here before is gone: two calls selecting one name, writing one file
/// between them, and one of them publishing the other's bytes at this call's
/// destination. The two renames still land in some order, and the destination
/// keeps whichever went last. The temporary's name has no bearing on that; the
/// contention there is over the destination.
///
/// The second is substitution at the temporary's own name. The create claims
/// the name; the fclose gives up the file. Every step after it -- OEWriteGrid,
/// the header patch, each verification re-read, the rename -- addresses the
/// temporary by path, because OEWriteGrid takes a path and not a descriptor, so
/// nothing holds the file across them. A process able to write the
/// destination's directory can therefore unlink the entry this function
/// verified, put its own file at that name, and have the rename publish those
/// bytes while write_map returns success. Refusing a second exclusive create is
/// what the reservation buys; it does not keep the name bound to the file it
/// created. Such a process can already create, replace and remove the
/// destination itself, so what the window adds is the success return on bytes
/// the verification never saw -- which is a statement about what write_map
/// promises its caller, and is made there.
TemporaryFile MakeTemporarySibling(const std::filesystem::path& dest) {
    constexpr int MAX_ATTEMPTS = 8;
    std::random_device entropy;
    for (int attempt = 0; attempt < MAX_ATTEMPTS; ++attempt) {
        std::ostringstream name;
        name.imbue(std::locale::classic());
        name << ".maptitude-" << ProcessId() << '-' << std::hex << entropy()
             << dest.extension().string();
        std::filesystem::path candidate = dest.parent_path() / name.str();
        // Narrowed into a named local rather than a temporary, so that no
        // destructor -- and so no deallocation, of exactly the class the next
        // comment is about -- runs between the create and the errno read.
        const std::string candidate_name = candidate.string();
        std::FILE* const reserved = std::fopen(candidate_name.c_str(), "wx");
        // errno is read once, into a local: the throw below builds its message
        // with operator+, whose operands are evaluated in an unspecified order
        // and whose allocations can overwrite errno before a second read of it
        // reaches the message.
        const int failure = errno;
        if (reserved != nullptr) {
            std::fclose(reserved);
            return TemporaryFile(std::move(candidate));
        }
        if (failure != EEXIST) {
            throw GridError(
                "Cannot write '" + dest.string() +
                "': cannot create a temporary beside it: " +
                std::error_code(failure, std::generic_category()).message());
        }
    }
    throw GridError("Cannot write '" + dest.string() +
                    "': no free temporary name beside it after " +
                    std::to_string(MAX_ATTEMPTS) + " attempts");
}

/// Normalize symop text to the one-triplet-per-line form read_map returns.
///
/// The boundaries have to be the ones SymOp::ParseAll uses, which are newline
/// *and* semicolon -- write_map validates the caller's text through that parser,
/// and SymopBlockBytes emits one 80-byte CCP4 record per line this function
/// hands back. Splitting on newlines alone accepts "x,y,z;-x,y,-z" as two
/// operators and then writes it as a single record with an embedded semicolon.
/// A consumer that reads the fixed-width records literally gets one malformed
/// operator out of that record, not two well-formed ones.
///
/// Normalizing rather than refusing, because the semicolon spelling is one the
/// library takes in, not one it emits. Both consumers checked accept it:
/// SymOp::ParseAll splits on it, and fc_density documents its symops argument
/// as accepting it, then consumes the string through ParseAll and returns a
/// grid. No semicolon join exists under src/, include/, python/ or swig/;
/// read_map joins with '\n'. Refusing the spelling here alone would therefore
/// make the writer stricter than both of those consumers over a spelling
/// nothing on those four paths produces. The cost is that the spelling does
/// not survive the round trip: a caller who writes "a;b" reads back "a\nb".
/// The operator set does.
std::string CanonicalSymops(const std::string& symops) {
    std::string joined;
    std::size_t start = 0;
    while (start <= symops.size() && !symops.empty()) {
        const std::size_t end = symops.find_first_of("\n;", start);
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
///
/// @p destination is the caller's path, named in failures instead of @p file:
/// @p file is the randomly named hidden sibling, which the caller never chose
/// and whose destructor has removed it by the time the message is read.
std::string ReadSymopBlock(const std::filesystem::path& file,
                           const std::size_t count,
                           const std::string& destination) {
    std::string block(count, '\0');
    if (count == 0u) {
        return block;
    }
    std::ifstream in(file, std::ios::binary);
    if (!in) {
        throw GridError("Cannot reopen the map written for '" + destination +
                        "' to verify its symmetry block");
    }
    in.seekg(static_cast<std::streamoff>(CCP4_HEADER_BYTES));
    in.read(&block[0], static_cast<std::streamsize>(count));
    if (in.gcount() != static_cast<std::streamsize>(count)) {
        throw GridError("The map written for '" + destination + "' declares " +
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

/// Splice NSYMBT, the symop block and ORIGIN into a file OEWriteGrid produced,
/// and zero the SKWTRN words it leaves uninitialized.
///
/// NX/NY/NZ is deliberately not patched: section 2.4 measures the dim - 1
/// patch an earlier draft applied here as a no-op on whole-cell grids and a
/// silent corruption on a sub-box, where it preserves dim, cell and voxel
/// count while moving the spacing and the grid's position.
///
/// @p destination is the caller's path, named in failures instead of @p file,
/// for the reason given on ReadSymopBlock.
void PatchHeaderRecords(const std::filesystem::path& file,
                        const std::string& canonical_symops,
                        const OESystem::OESkewGrid& grid,
                        const std::string& destination) {
    std::string bytes;
    {
        std::ifstream in(file, std::ios::binary);
        if (!in) {
            throw GridError("Cannot reopen the map written for '" +
                            destination + "' to restore its header records");
        }
        bytes.assign(std::istreambuf_iterator<char>(in),
                     std::istreambuf_iterator<char>());
    }
    if (bytes.size() < CCP4_HEADER_BYTES) {
        throw GridError("OEWriteGrid produced a map for '" + destination +
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
        throw GridError("OEWriteGrid produced a map for '" + destination +
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
                        "ORIGIN record of the map for '" + destination + "'");
    }
    if (node_x != 0.0f || node_y != 0.0f || node_z != 0.0f) {
        put_float(WORD_ORIGIN + 0u, node_x);
        put_float(WORD_ORIGIN + 1u, node_y);
        put_float(WORD_ORIGIN + 2u, node_z);
    }

    // OEWriteGrid leaves SKWTRN, words 35-37, holding whatever was in the
    // memory behind them. Measured across three processes writing one asset,
    // word 35 carried the low half of a heap pointer under ASLR, so two writes
    // of the same grid produced different files and four bytes of the writing
    // process's address space landed in a file the caller may share. It is not
    // something this module puts there: a grid handed straight from OEReadGrid
    // comes back with those words filled with spaces, while a bare
    // copy-construct of that same grid, with no maptitude code involved, shows
    // the drift -- and write_map copies the caller's grid to default its space
    // group. Zero is the semantically correct SKWTRN and is what the four
    // source assets under tests/assets/mapq carry, so a map read from one of
    // them and written back keeps those words unchanged. tests/data/test_map.ccp4
    // carries the same junk pattern at those words and is normalized to zero on
    // the way through.
    //
    // Guarded on LSKFLG rather than unconditional, so a future path that does
    // write a real skew translation does not have it silently erased. The guard
    // does not cost this path its fix: LSKFLG was 0 in every file measured, and
    // write_map refuses with CellError the non-axis-aligned sampling a nonzero
    // one would describe.
    if (WordAsInt(bytes, WORD_LSKFLG, swapped) == 0) {
        put_float(WORD_SKWTRN + 0u, 0.0f);
        put_float(WORD_SKWTRN + 1u, 0.0f);
        put_float(WORD_SKWTRN + 2u, 0.0f);
    }

    std::ofstream out(file, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw GridError("Cannot rewrite the map written for '" + destination +
                        "' with its restored header records");
    }
    out.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    out.close();
    if (!out) {
        throw GridError("Failed while rewriting the map written for '" +
                        destination + "'");
    }
}

}  // namespace

MapFile::MapFile() = default;
MapFile::MapFile(MapFile&& other) noexcept = default;
MapFile& MapFile::operator=(MapFile&& other) noexcept = default;
MapFile::~MapFile() = default;

MapFile read_map(const std::string& path, const OriginSource tiebreak) {
    RejectEmbeddedNul(path, "read");
    if (IsCompressedPath(path)) {
        throw GridError("Cannot read '" + path +
                        "': this names a compressed file, and OEReadGrid "
                        "decompresses it while the header and symmetry block "
                        "are read from the bytes on disk, so the two halves of "
                        "this function would be reading two different streams");
    }
    // Extension handling as the write gate does it, for the reasons measured
    // there: no leading dot, and no lowercasing.
    std::string extension = std::filesystem::path(path).extension().string();
    if (!extension.empty()) {
        extension.erase(0, 1u);
    }
    if (OESystem::OEGetGridFileType(extension.c_str()) !=
        OESystem::OEGridFileType::CCP4) {
        // OEReadGrid alone would take this file: it reads every format OpenEye
        // reads. What cannot is everything after it -- ReadRawHeader and
        // ParseHeader interpret the file's own first 1024 bytes as a CCP4
        // header, and ReadSymopBlock reads NSYMBT bytes past that as 80-byte
        // symmetry records. Nothing downstream catches the result: the read
        // path has no counterpart to compare_placement_records, which lives
        // only in the write verify.
        //
        // The NSYMBT sanity check is what refuses such a file today, and it
        // refuses on an accident of the bytes rather than on the format. Each
        // of the three formats OpenEye writes was written from the grid
        // read_map returns for tests/assets/mapq/1d26_2fofc.ccp4 and read back:
        // '.grd', '.agd' and '.phi' all land on a word 24 that is not a
        // multiple of 80, so all three are refused there. Zeroing that one word
        // in the '.phi' -- and nothing else -- was enough for read_map to
        // accept it and hand back the 65x65x65 grid OEReadGrid made of the
        // Grasp stream, 19.5 A from the node 0 of the map it was written from,
        // with no error. What the header half does when its bytes are not zero
        // is visible in the '.agd' written from that same grid, whose words
        // 50-52 are the text "502e-02\n-5.0" and read as an ORIGIN of
        // (5.3e22, 8.6e-33, 6.3e-10).
        //
        // Which grid was written matters for the byte figures above, so it is
        // stated rather than left to the asset's name: a bare OEReadGrid of
        // that same asset writes a different '.agd' -- "7e-02\n2.7225" at those
        // words -- and a '.grd' whose word 24 is 1049322889 rather than 1.
        // Neither is a multiple of 80 either, so the sentence about the three
        // formats holds for both grids; the quoted bytes do not.
        //
        // UNDEFINED is not CCP4, so an extension OpenEye does not recognize is
        // covered by the same comparison.
        throw GridError("Cannot read '" + path +
                        "': this reader parses the file's own first 1024 bytes "
                        "as a CCP4 header, so the extension has to be one "
                        "OpenEye maps to CCP4; read '.ccp4', '.mrc' or '.map'");
    }

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
    RejectEmbeddedNul(path, "write");
    SymOp::ParseAll(symops);
    if (IsCompressedPath(path)) {
        throw GridError("Cannot write '" + path +
                        "': this names a compressed file, and the CCP4 header "
                        "records this writer restores after OEWriteGrid cannot "
                        "be spliced into a compressed stream");
    }
    // OEGetGridFileType wants the extension with no leading dot: over an
    // eleven-spelling sweep of leading-dot forms every one read as UNDEFINED,
    // ".ccp4" and ".mrc" among them, while "ccp4." read as CCP4. It is
    // case-insensitive over the mixed-case spellings swept -- "ccpQ", "MAPX"
    // and "MaPq" all read as CCP4 -- so lowercasing here would be dead code.
    std::string extension = std::filesystem::path(path).extension().string();
    if (!extension.empty()) {
        extension.erase(0, 1u);
    }
    if (OESystem::OEGetGridFileType(extension.c_str()) !=
        OESystem::OEGridFileType::CCP4) {
        // Section 2.4 measures that OEWriteGrid on an extension OpenEye does
        // not recognize prints a fatal error and exits rather than returning
        // false. This check is still the only thing between the caller and a
        // killed process; once OEWriteGrid has the path there is nothing to
        // catch. UNDEFINED is not CCP4, so it covers that case as the old
        // OEIsWriteableGrid test did.
        //
        // It is narrower than that test, deliberately. A format OpenEye writes
        // but this path cannot patch -- ".grd" is the measured one -- passed
        // the old gate, so OEWriteGrid produced a GRD file and
        // PatchHeaderRecords spliced CCP4-offset bytes into it before the
        // verify refused. The refusal was atomic, but it blamed the map for not
        // reading back when the fault was the extension.
        throw GridError("Cannot write '" + path +
                        "': this writer restores CCP4 header records after "
                        "OEWriteGrid, so the extension has to be one OpenEye "
                        "maps to CCP4; write '.ccp4', '.mrc' or '.map'");
    }

    // Step 2: copy, and default the space group. Section 2.3 shows an unset
    // space group is what makes OEWriteGrid double the cell and regrid.
    OESystem::OESkewGrid out(grid);
    if (!out.HasSpaceGroup() && !out.SetSpaceGroup(1u)) {
        throw GridError("Cannot default the space group to P1 before "
                        "writing '" + path + "'");
    }

    const std::filesystem::path dest(path);
    // MakeTemporarySibling creates the file it names, so it hands back the
    // owner rather than a path a caller has to remember to wrap.
    TemporaryFile temporary = MakeTemporarySibling(dest);

    // Step 3.
    if (!OESystem::OEWriteGrid(temporary.Path().string(), out)) {
        throw GridError("OEWriteGrid failed while writing '" + path + "'");
    }

    // Step 4.
    const std::string canonical = CanonicalSymops(symops);
    PatchHeaderRecords(temporary.Path(), canonical, out, path);

    // Step 5: verify against a read_map re-read. The round trip is closed
    // under read_map, not under OEReadGrid: the bare reader never consults
    // ORIGIN, so a map written with a nonzero one lands wrong by construction.
    //
    // read_map names the file it was handed, which here is the temporary, so
    // its GridErrors are re-raised against the destination with that path
    // substituted out of the reason. The outer clause is what tells the reader
    // the inner message is about the map just written rather than about a
    // destination that may not exist yet.
    MapFile back;
    try {
        back = read_map(temporary.Path().string());
    } catch (const GridError& error) {
        throw GridError("Refusing to write '" + path +
                        "': the map written for it does not read back as a "
                        "map: " +
                        MessageAgainstDestination(error.what(),
                                                  temporary.Path().string(),
                                                  path));
    }
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
    std::string raw;
    try {
        raw = ReadRawHeader(temporary.Path().string());
    } catch (const GridError& error) {
        // Same treatment as the read_map above: the two failures ReadRawHeader
        // reports both quote the temporary's path.
        throw GridError("Refusing to write '" + path +
                        "': the header of the map written for it cannot be "
                        "re-read: " +
                        MessageAgainstDestination(error.what(),
                                                  temporary.Path().string(),
                                                  path));
    }
    const MapHeader written = ParseHeader(raw);

    // back.symops above compares the parsed text. This compares the bytes:
    // NSYMBT and the fixed-width record shape are what a consumer that is not
    // read_map will read, and the text comparison passes on a block whose
    // padding or record count is wrong but which strips to the same triplets.
    //
    // Both sides below are derived from `canonical`, so this guard is only ever
    // as good as CanonicalSymops: it checks that the writer emitted the record
    // boundaries that function chose, never that those boundaries are the right
    // ones. That is why CanonicalSymops has to split on exactly what
    // SymOp::ParseAll splits on, and why the test for it asserts on raw bytes
    // written from hand-built symop text rather than through this check.
    const std::string expected_block = SymopBlockBytes(canonical);
    if (static_cast<std::size_t>(written.nsymbt) != expected_block.size()) {
        throw GridError("Refusing to write '" + path + "': NSYMBT is " +
                        std::to_string(written.nsymbt) + " but the symmetry "
                        "block is " + std::to_string(expected_block.size()) +
                        " bytes");
    }
    if (ReadSymopBlock(temporary.Path(), expected_block.size(), path) !=
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
    MapFile by_nxstart;
    try {
        by_nxstart = read_map(temporary.Path().string(),
                              OriginSource::NXSTART);
    } catch (const GridError& error) {
        throw GridError("Refusing to write '" + path +
                        "': the map written for it cannot be re-read under "
                        "the NxSTART placement: " +
                        MessageAgainstDestination(error.what(),
                                                  temporary.Path().string(),
                                                  path));
    }
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
