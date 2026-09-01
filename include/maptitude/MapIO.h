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

}  // namespace Maptitude

#endif  // MAPTITUDE_MAPIO_H
