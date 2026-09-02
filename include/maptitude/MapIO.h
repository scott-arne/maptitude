/**
 * @file MapIO.h
 * @brief CCP4/MRC map reading and writing that preserves the header records
 *        the OpenEye grid I/O discards.
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
 * @brief Read a CCP4 or MRC map, preserving the records OEReadGrid drops.
 *
 * @param path Map file to read.
 * @param tiebreak Which record wins when ORIGIN and NxSTART both encode a
 *        nonzero, differing origin. Ignored when at most one is nonzero.
 * @return The grid, positioned at the file's origin, and its symmetry text.
 * @throws GridError if the file cannot be read or its header cannot be parsed.
 * @throws SymOpError if the symmetry block is present and does not parse.
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
/// Comparing against the half exactly would put that case on the boundary,
/// where float32 noise in the ORIGIN record decides it. This widens the bound
/// by a hundredth of a percent, which is far below any genuine disagreement --
/// the smallest of those is a whole node interval.
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
 * @brief Write a grid as CCP4/MRC, restoring the records OEWriteGrid drops.
 *
 * Writes and verifies a temporary file, then renames it onto @p path. A grid
 * this function cannot write faithfully raises with @p path untouched, rather
 * than leaving a silently wrong map on disk or destroying a good one.
 *
 * Node 0 is written twice, into the MRC2000 ORIGIN record exactly and into
 * NxSTART as an integer node count, and the write is refused unless the two
 * agree on every axis to within half a node interval. So a map this function
 * wrote reproduces node 0 to float32 under the default
 * OriginSource::ORIGIN_RECORD, and only to half a node interval under
 * OriginSource::NXSTART. That gap is not a defect to be closed: an integer
 * lattice index cannot encode an off-lattice origin. The refusal is what makes
 * the bound a guarantee instead of an observation about one toolkit version.
 *
 * @param path Destination. The format follows the extension, as OpenEye
 *        dispatches on it; an extension OpenEye does not recognize is
 *        rejected before anything is written, as is a compressed destination
 *        such as ".ccp4.gz", whose stream has no place to splice the header
 *        records this function restores.
 * @param grid Grid to write.
 * @param symops Symmetry text in the form read_map returns: one triplet per
 *        line, newline-separated. Empty writes no block.
 * @throws GridError if @p path's extension is not a grid format OpenEye
 *         writes or names a compressed file, if the write fails, if the re-read
 *         map differs from the grid in dimensions, cell, per-axis spacing,
 *         node 0, or any voxel, if its ORIGIN and NxSTART records place node 0
 *         more than half a node interval apart, or if the temporary file
 *         cannot be renamed onto @p path.
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
