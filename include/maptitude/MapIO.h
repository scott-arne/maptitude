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
 * spacing and the position of node 0 compare relatively at MAP_VERIFY_TOL.
 *
 * @param expected The grid handed to write_map.
 * @param actual The grid read back from the written file.
 * @return An empty string when the two agree, otherwise a sentence naming the
 *         first quantity that differed and both of its values.
 */
std::string compare_written_map(const OESystem::OESkewGrid& expected,
                                const OESystem::OESkewGrid& actual);

/**
 * @brief Write a grid as CCP4/MRC, restoring the records OEWriteGrid drops.
 *
 * Writes and verifies a temporary file, then renames it onto @p path. A grid
 * this function cannot write faithfully raises with @p path untouched, rather
 * than leaving a silently wrong map on disk or destroying a good one.
 *
 * @param path Destination. The format follows the extension, as OpenEye
 *        dispatches on it; an extension OpenEye does not recognize is
 *        rejected before anything is written.
 * @param grid Grid to write.
 * @param symops Symmetry text in the form read_map returns: one triplet per
 *        line, newline-separated. Empty writes no block.
 * @throws GridError if @p path's extension is not a grid format OpenEye
 *         writes, if the write fails, if the re-read map differs from the
 *         grid in dimensions, cell, per-axis spacing, node 0, or any voxel, or
 *         if the temporary file cannot be renamed onto @p path.
 * @throws SymOpError if symops is non-empty and does not parse.
 */
void write_map(const std::string& path,
               const OESystem::OESkewGrid& grid,
               const std::string& symops = "");

}  // namespace Maptitude

#endif  // MAPTITUDE_MAPIO_H
