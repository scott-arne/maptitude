#include "maptitude/SpatialIndex.h"

#include <oechem.h>
#include <nanoflann.hpp>

#include <vector>

namespace Maptitude {

// Point cloud adaptor for nanoflann
struct PointCloud {
    std::vector<double> coords;  // Flat: x0,y0,z0, x1,y1,z1, ...
    size_t num_points = 0;

    inline size_t kdtree_get_point_count() const { return num_points; }

    inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
        return coords[idx * 3 + dim];
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false; }
};

using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, PointCloud>,
    PointCloud, 3>;

struct SpatialIndex::Impl {
    PointCloud cloud;
    std::unique_ptr<KDTree> tree;
    std::vector<unsigned int> atom_indices;  // Maps tree index -> atom index
};

SpatialIndex::SpatialIndex(OEChem::OEMolBase& mol)
    : pimpl_(std::make_unique<Impl>()) {
    // Count atoms and allocate
    unsigned int n = mol.NumAtoms();
    pimpl_->cloud.coords.reserve(n * 3);
    pimpl_->atom_indices.reserve(n);

    float coords[3];
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        mol.GetCoords(&(*atom), coords);
        pimpl_->cloud.coords.push_back(coords[0]);
        pimpl_->cloud.coords.push_back(coords[1]);
        pimpl_->cloud.coords.push_back(coords[2]);
        pimpl_->atom_indices.push_back(atom->GetIdx());
    }
    pimpl_->cloud.num_points = pimpl_->atom_indices.size();

    // Build k-d tree
    pimpl_->tree = std::make_unique<KDTree>(
        3, pimpl_->cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    pimpl_->tree->buildIndex();
}

SpatialIndex::~SpatialIndex() = default;

std::vector<unsigned int> SpatialIndex::FindWithinRadius(
    double x, double y, double z, double radius) const {
    double query_pt[3] = {x, y, z};
    double search_radius_sq = radius * radius;

    nanoflann::SearchParameters params;
    params.sorted = false;

    std::vector<nanoflann::ResultItem<uint32_t, double>> matches;
    pimpl_->tree->radiusSearch(query_pt, search_radius_sq, matches, params);

    std::vector<unsigned int> result;
    result.reserve(matches.size());
    for (const auto& match : matches) {
        result.push_back(pimpl_->atom_indices[match.first]);
    }

    return result;
}

std::vector<unsigned int> SpatialIndex::FindWithinRadius(
    const OEChem::OEAtomBase& atom, double radius) const {
    // Look up stored coordinates by atom index
    unsigned int target_idx = atom.GetIdx();
    for (size_t i = 0; i < pimpl_->atom_indices.size(); ++i) {
        if (pimpl_->atom_indices[i] == target_idx) {
            return FindWithinRadius(
                pimpl_->cloud.coords[i * 3],
                pimpl_->cloud.coords[i * 3 + 1],
                pimpl_->cloud.coords[i * 3 + 2],
                radius);
        }
    }
    return {};
}

size_t SpatialIndex::Size() const {
    return pimpl_->cloud.num_points;
}

}  // namespace Maptitude
