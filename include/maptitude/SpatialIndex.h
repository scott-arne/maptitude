/**
 * @file SpatialIndex.h
 * @brief K-d tree spatial index for efficient radius queries.
 *
 * SpatialIndex wraps nanoflann to provide O(log n) radius queries
 * used by Q-score point isolation and distance-based operations.
 */

#ifndef MAPTITUDE_SPATIALINDEX_H
#define MAPTITUDE_SPATIALINDEX_H

#include <cstddef>
#include <memory>
#include <vector>

namespace OEChem {
class OEMolBase;
class OEAtomBase;
}

namespace Maptitude {

/**
 * @brief K-d tree based spatial index for efficient atom distance queries.
 *
 * Builds a 3D k-d tree from atom coordinates for fast radius searches.
 * Used internally by Q-score for isolating shell sample points from
 * neighboring atoms.
 *
 * @code
 * SpatialIndex idx(mol);
 * auto nearby = idx.FindWithinRadius(10.0, 20.0, 30.0, 2.5);
 * @endcode
 */
class SpatialIndex {
public:
    /**
     * @brief Construct spatial index from molecule coordinates.
     *
     * Builds a k-d tree from all atom positions. O(n log n) construction.
     *
     * @param mol The molecule to index.
     */
    explicit SpatialIndex(OEChem::OEMolBase& mol);

    /// @brief Destructor.
    ~SpatialIndex();

    // Non-copyable
    SpatialIndex(const SpatialIndex&) = delete;
    SpatialIndex& operator=(const SpatialIndex&) = delete;

    /**
     * @brief Find all atoms within radius of a point.
     *
     * @param x X coordinate of query point.
     * @param y Y coordinate of query point.
     * @param z Z coordinate of query point.
     * @param radius Maximum distance in Angstroms.
     * @return Vector of atom indices within the radius.
     */
    [[nodiscard]] std::vector<unsigned int> FindWithinRadius(
        double x, double y, double z, double radius) const;

    /**
     * @brief Find all atoms within radius of another atom.
     *
     * @param atom Reference atom for the query.
     * @param radius Maximum distance in Angstroms.
     * @return Vector of atom indices within the radius.
     */
    [[nodiscard]] std::vector<unsigned int> FindWithinRadius(
        const OEChem::OEAtomBase& atom, double radius) const;

    /**
     * @brief Get the number of indexed atoms.
     * @return Number of atoms in the index.
     */
    [[nodiscard]] size_t Size() const;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

}  // namespace Maptitude

#endif  // MAPTITUDE_SPATIALINDEX_H
