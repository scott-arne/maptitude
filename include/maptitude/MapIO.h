/**
 * @file MapIO.h
 * @brief CCP4/MRC map reading and writing that preserves two header records
 *        the OpenEye grid I/O discards: the MRC2000 ORIGIN placement and the
 *        symmetry block.
 *
 * OEReadGrid returns the payload and the cell but drops the MRC2000 ORIGIN
 * record and the symmetry block; OEWriteGrid emits neither. Both are needed to
 * place an EM map correctly and to keep a crystallographic map's symmetry with
 * it, so this module reads and writes them alongside the toolkit calls.
 */
#ifndef MAPTITUDE_MAPIO_H
#define MAPTITUDE_MAPIO_H

#include <memory>
#include <string>

namespace OESystem {
class OESkewGrid;
}

namespace Maptitude {

/// Which header record wins when both encode a nonzero origin.
enum class OriginSource {
    ORIGIN_RECORD,  ///< Header words 50-52, the MRC2000 ORIGIN.
    NXSTART         ///< Words 5-7, as the reader already applied them.
};

/// The contents of a map file: the density carrier and its symmetry text.
///
/// The carrier holds cell edges, cell angles and the space group natively,
/// so this bundle adds only what OEReadGrid discards.
struct MapFile {
    MapFile();
    MapFile(MapFile&& other) noexcept;
    MapFile& operator=(MapFile&& other) noexcept;
    ~MapFile();

    std::unique_ptr<OESystem::OESkewGrid> grid;
    std::string symops;  ///< Empty when the file carries no symmetry block.
};

/**
 * @brief Read a CCP4 or MRC map, preserving the ORIGIN record and symmetry
 *        block OEReadGrid drops.
 *
 * @param path Map file to read.
 * @param tiebreak Which record wins when ORIGIN and NxSTART both encode a
 *        nonzero origin, whether or not the two agree -- nothing here compares
 *        them. Ignored when at most one is nonzero.
 * @return The grid, positioned at the file's origin, and its symmetry text.
 * @throws GridError if the file cannot be read, if its header cannot be parsed,
 *         or if @p path contains an embedded NUL, since the calls that open the
 *         file stop at the NUL, so the name this function is given and the name
 *         it would open are not the same name.
 * @throws SymOpError if the symmetry block is present and does not parse, or if
 *         one of its 80-byte records still holds a ';' or a newline once its
 *         padding is stripped. A CCP4 record carries one operator, so a
 *         separator inside one separates nothing on disk; a record holding two
 *         triplets that way is refused rather than read as two operators, which
 *         would give a write_map round trip one more record than the file has.
 */
MapFile read_map(const std::string& path,
                 OriginSource tiebreak = OriginSource::ORIGIN_RECORD);

/// Relative tolerance for the geometry comparison write_map verifies with.
///
/// The ORIGIN record is float32, so a node-0 position that is not exactly
/// representable comes back quantized. Section 2.4 measures the worst residual
/// over five placements at 1.34e-7 relative, which this clears by a factor of
/// about seven. An exact comparison here would refuse grids this path writes
/// correctly.
constexpr double MAP_VERIFY_TOL = 1e-6;

/**
 * @brief Compare a re-read map against the grid it was written from.
 *
 * The predicate behind write_map's verification step, exposed so a test can
 * drive it directly with two deliberately different grids rather than having to
 * manufacture a file that fails on each quantity.
 *
 * Node counts and voxel values compare exactly; cell parameters, per-axis
 * spacing and the position of node 0 compare relatively at MAP_VERIFY_TOL. Two
 * NaN voxels count as agreeing, since a masked voxel that round-tripped
 * faithfully is not a difference; a NaN against a number is.
 *
 * Two kinds of input are reported by throwing rather than by the return value,
 * because they leave the function unable to describe either grid rather than
 * describing a difference between them.
 *
 * @param expected The grid handed to write_map.
 * @param actual The grid read back from the written file.
 * @return An empty string when the two agree, otherwise a sentence naming the
 *         first quantity that differed and both of its values.
 * @throws CellError if either grid's sampling is not axis-aligned, which a
 *         cell angle away from 90 degrees produces. read_map can return such a
 *         grid; this comparison cannot express one.
 * @throws GridError if either grid has an axis with fewer than two nodes, so
 *         its spacing is undefined.
 */
std::string compare_written_map(const OESystem::OESkewGrid& expected,
                                const OESystem::OESkewGrid& actual);

/// Fractional slack above the half-node-interval placement bound.
///
/// The bound below which the two placement records may legitimately disagree is
/// exactly half a node interval, and it is inclusive: OpenEye rounds
/// ORIGIN/spacing to the nearest node with ties away from zero, so a node 0 at
/// 2.75 with a spacing of 0.5 lands exactly on the half and is a correct write.
///
/// Comparing against the half exactly refuses a subset of those correct writes.
/// The NxSTART placement is reconstructed as NCSTART times the spacing and held
/// in float32; where that product is not exactly representable the reconstructed
/// position lands up to half a float32 ULP off, which pushes the divergence to
/// either side of the half. Measured over fifteen division counts on a dim-20,
/// cell-20 grid whose node 0 is a half-integer of spacings: six land above the
/// half, the worst at 9.4e-7 relative, and all six are refused with this slack
/// at zero. The ORIGIN record is not the source -- at those placements it
/// round-trips node 0 bit-exactly. This slack clears the worst measured
/// excursion by about a factor of a hundred, and stays four orders below the
/// smallest genuine disagreement, which is a whole node interval.
constexpr double MAP_PLACEMENT_SLACK = 1e-4;

/**
 * @brief Compare the two placements a written map's header records encode.
 *
 * A CCP4/MRC header states node 0 twice: exactly, in the MRC2000 ORIGIN record,
 * and as an integer node count, in NxSTART. read_map returns one or the other
 * depending on its OriginSource argument, so the two must agree to within the
 * rounding an integer index forces -- half a node interval -- or the same file
 * places its density differently for two callers who both used the public API.
 *
 * The predicate behind write_map's placement check, exposed so a test can drive
 * it with two deliberately disagreeing grids rather than having to manufacture
 * a header whose NxSTART is wrong.
 *
 * Both arguments must describe the same sampling; the node interval compared
 * against is taken from @p by_origin.
 *
 * @param by_origin The file read back under OriginSource::ORIGIN_RECORD.
 * @param by_nxstart The same file read back under OriginSource::NXSTART.
 * @return An empty string when every axis agrees to within half a node interval
 *         (inclusive, plus MAP_PLACEMENT_SLACK), otherwise a sentence naming
 *         the first axis that did not and both of its positions.
 * @throws CellError if either grid's sampling is not axis-aligned.
 * @throws GridError if either grid has an axis with fewer than two nodes.
 */
std::string compare_placement_records(const OESystem::OESkewGrid& by_origin,
                                      const OESystem::OESkewGrid& by_nxstart);

/**
 * @brief Write a grid as CCP4/MRC, restoring the ORIGIN record and symmetry
 *        block OEWriteGrid drops.
 *
 * Writes and verifies a temporary file, then renames it onto @p path. A grid
 * this function cannot write faithfully raises with @p path untouched, rather
 * than leaving a silently wrong map on disk or destroying a good one.
 *
 * A successful write replaces the destination's inode rather than rewriting it
 * in place, so three properties of an existing destination do not survive. Its
 * mode resets to whatever an ordinary create gives, which under POSIX is 0666
 * narrowed by the process umask, widening a permission the caller had
 * narrowed. Hard links to it keep the old contents under their own names. A
 * destination that was a symlink becomes a regular file, with its former
 * target left untouched.
 *
 * The verification covers the bytes this function wrote, not the bytes that
 * arrive at @p path. It writes and checks a temporary beside the destination
 * and then renames that onto @p path, and every step addresses the temporary by
 * path: the toolkit's writer takes a path rather than a descriptor, and the
 * rename that publishes it takes one too. So a process able to write the
 * destination's directory can replace the temporary between the last check and
 * the rename, and this function will publish its bytes and return successfully.
 * Such a process can already create, replace and remove @p path itself; what
 * this adds is that a successful return stops implying the published bytes are
 * the ones that were verified.
 *
 * A destination directory only the caller can write closes that. Confining the
 * temporary to a directory only the caller can enter would narrow it to the
 * same set, and is not what this function does. Nothing closes it outright: the
 * publication step is a rename, and rename names its source by path.
 *
 * Node 0 is written twice, into the MRC2000 ORIGIN record exactly and into
 * NxSTART as an integer node count, and the write is refused unless the two
 * agree on every axis to within half a node interval. So a map this function
 * wrote reproduces node 0 to float32 under the default
 * OriginSource::ORIGIN_RECORD, and only to half a node interval per axis under
 * OriginSource::NXSTART. That gap is not a defect to be closed: an integer
 * lattice index cannot encode an off-lattice origin. The refusal is what makes
 * the bound a guarantee instead of an observation about one toolkit version.
 *
 * That bound is per axis, and there is no second one over the three together:
 * each axis is compared against its own node interval. So the straight-line
 * distance between the two placements can reach the root-sum-square of the
 * three half-intervals -- sqrt(3) times half a node interval on a grid sampled
 * equally on all three axes -- and a caller budgeting one distance rather than
 * three per-axis bounds needs that larger figure.
 *
 * @param path Destination. The format follows the extension, though not
 *        because OEWriteGrid reads the extension: it dispatches on the whole
 *        filename, and a measured name carrying a "gz" dot-component ahead of
 *        ".ccp4" came out as a gzip stream. The format follows the extension
 *        because the gate here reads the extension, and the temporary handed
 *        to OEWriteGrid keeps @p path's extension and takes nothing else of
 *        its name, so that whole-filename dispatch lands on the format the
 *        gate already admitted. This writer admits whatever OEGetGridFileType
 *        maps to CCP4. In a thirty-spelling sweep that class was every
 *        spelling whose first three characters are "ccp", "map" or "mrc",
 *        compared case-insensitively, and nothing else: ".ccp4junk" and
 *        ".mrcs" are admitted and written as CCP4 exactly as ".ccp4" is.
 *        ".ccp4", ".mrc" and ".map" are the spellings this writer is tested
 *        on, not the accepted set. An extension OpenEye maps elsewhere is
 *        rejected before anything is written, because the header records this
 *        function restores after OEWriteGrid are at CCP4 offsets; so is a
 *        compressed destination such as ".ccp4.gz", whose stream has no place
 *        to splice them. The extension is whatever std::filesystem reports, so
 *        a filename that is nothing but an extension -- ".ccp4" -- has none
 *        and is refused.
 * @param grid Grid to write.
 * @param symops Symmetry text in the form read_map returns: one triplet per
 *        line, newline-separated. Empty writes no block. The semicolon
 *        separator SymOp::ParseAll also accepts is normalized to a newline, so
 *        the operator set survives the round trip but that spelling does not:
 *        writing "a;b" reads back as "a\nb".
 * @throws GridError if @p path's extension is not one OpenEye maps to the CCP4
 *         format or names a compressed file, if @p path contains an embedded
 *         NUL, since std::filesystem reads the whole string while the calls
 *         that write and publish the file stop at the NUL, so the name checked
 *         here and the name written are not the same name, if the temporary
 *         this function writes beside @p path cannot be created, which an
 *         unwritable or absent destination directory produces, if eight
 *         attempts at a temporary name beside @p path all collide with an
 *         existing file, if the write fails, if the re-read map differs from
 *         the grid in dimensions, cell, per-axis spacing, node 0, or any voxel,
 *         if either the grid or the re-read map has an axis with fewer than two
 *         nodes, so its spacing is undefined, if on any axis its ORIGIN and
 *         NxSTART records place node 0 more than half that axis's node interval
 *         apart, or if the temporary file cannot be renamed onto @p path. That
 *         list is not closed: this function raises GridError from further
 *         internal checks, among them the header re-reads and the NSYMBT and
 *         symop-block byte comparisons. Catch the class rather than switching
 *         on the list.
 * @throws SymOpError if symops is non-empty and does not parse, or if any of
 *         its records is longer than the format's 80-character field.
 * @throws CellError if the grid's sampling is not axis-aligned, which a cell
 *         angle away from 90 degrees produces. read_map can return such a grid;
 *         this writer cannot express one.
 */
void write_map(const std::string& path,
               const OESystem::OESkewGrid& grid,
               const std::string& symops = "");

}  // namespace Maptitude

#endif  // MAPTITUDE_MAPIO_H
