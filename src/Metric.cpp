#include "maptitude/Metric.h"
#include "maptitude/Error.h"
#include "maptitude/Grid.h"
#include "maptitude/SpatialIndex.h"
#include "maptitude/detail/ScoringHelpers.h"

#include <oechem.h>
#include <oegrid.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <vector>

namespace Maptitude {

// Import detail helpers into this TU
using detail::MapStats;
using detail::PearsonCorrelation;
using detail::EDIAmSigmoid;
using detail::FibonacciSpherePoints;

// ---- OE-aware wrappers around detail helpers ----

static MapStats ComputeMapStats(const OESystem::OEScalarGrid& grid) {
    unsigned int size = grid.GetSize();
    if (size == 0) return {};
    std::vector<double> values(size);
    for (unsigned int i = 0; i < size; ++i) {
        values[i] = grid[i];
    }
    return detail::ComputeMapStats(values);
}

static void GetMapNormalization(const OESystem::OEScalarGrid& grid,
                                double& A, double& B) {
    unsigned int size = grid.GetSize();
    if (size == 0) { A = 1.0; B = 0.0; return; }
    std::vector<double> values(size);
    for (unsigned int i = 0; i < size; ++i) {
        values[i] = grid[i];
    }
    detail::GetMapNormalization(values.data(), values.size(), A, B);
}

// ---- Helper: collect atoms grouped by residue ----

static std::map<Residue, std::vector<const OEChem::OEAtomBase*>>
CollectAtomsByResidue(
    OEChem::OEMolBase& mol,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask) {
    std::map<Residue, std::vector<const OEChem::OEAtomBase*>> result;

    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() == 1) continue;
        if (mask && !(*mask)(*atom)) continue;

        Residue res = Residue::FromAtom(*atom);
        result[res].push_back(&(*atom));
    }

    return result;
}

// ---- Helper: get atom coordinates ----

static void GetAtomCoords(const OEChem::OEMolBase& mol,
                          const OEChem::OEAtomBase& atom,
                          double& x, double& y, double& z) {
    float coords[3];
    mol.GetCoords(&atom, coords);
    x = coords[0];
    y = coords[1];
    z = coords[2];
}

// ---- Helper: adaptive scoring radius (OE-aware wrapper) ----

static double ScoringRadius(const OEChem::OEAtomBase& atom, double resolution) {
    OEChem::OEResidue res = OEChem::OEAtomGetResidue(&atom);
    return detail::ScoringRadius(res.GetBFactor(), resolution);
}

// ---- Helper: assign Bondi VDW radii if missing ----

static void PrepareStructure(OEChem::OEMolBase& mol) {
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(OEChem::OEIsHeavy()); atom; ++atom) {
        if (atom->GetRadius() > 0.0) return;  // radii already assigned
    }
    OEChem::OEAssignBondiVdWRadii(mol);
}

// ==== Density scoring functions ====

DensityScoreResult RSCC(
    OEChem::OEMolBase& mol,
    const OESystem::OEScalarGrid& grid,
    double resolution,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask,
    const OESystem::OEScalarGrid* calc_grid,
    const RSCCOptions& options) {
    if (resolution <= 0.0) {
        throw GridError("Resolution must be positive");
    }
    PrepareStructure(mol);

    auto residue_atoms = CollectAtomsByResidue(mol, mask);
    if (residue_atoms.empty()) {
        throw StructureError("No scorable heavy atoms after applying mask");
    }

    // If no calc_grid provided, we would compute one via DensityCalculator.
    // For now, require calc_grid to be provided.
    if (calc_grid == nullptr) {
        throw GridError("calc_grid is required (auto-generation not yet supported)");
    }

    DensityScoreResult result;
    std::vector<double> all_obs, all_calc;

    for (const auto& [res, atoms] : residue_atoms) {
        std::vector<double> res_obs, res_calc;

        for (const auto* atom : atoms) {
            double x, y, z;
            GetAtomCoords(mol, *atom, x, y, z);

            if (!grid.IsInGrid(static_cast<float>(x),
                               static_cast<float>(y),
                               static_cast<float>(z))) {
                result.by_atom[atom->GetIdx()] =
                    std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            double radius;
            switch (options.GetAtomRadiusMethod()) {
                case AtomRadius::Fixed:
                    radius = options.GetFixedAtomRadius();
                    break;
                case AtomRadius::Scaled:
                    radius = atom->GetRadius() * options.GetAtomRadiusScaling();
                    if (radius < 0.1) radius = 1.5;
                    break;
                case AtomRadius::Binned:
                default:
                    radius = detail::BinnedAtomRadius(resolution);
                    break;
            }

            auto pts = GetAtomGridPoints(grid, x, y, z, radius);
            if (pts.empty()) {
                result.by_atom[atom->GetIdx()] =
                    std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            std::vector<double> obs_vals, calc_vals;
            obs_vals.reserve(pts.size());
            calc_vals.reserve(pts.size());

            for (unsigned int idx : pts) {
                obs_vals.push_back(grid[idx]);
                calc_vals.push_back((*calc_grid)[idx]);
            }

            // Per-atom RSCC
            result.by_atom[atom->GetIdx()] =
                PearsonCorrelation(obs_vals, calc_vals);

            res_obs.insert(res_obs.end(), obs_vals.begin(), obs_vals.end());
            res_calc.insert(res_calc.end(),
                calc_vals.begin(), calc_vals.end());
        }

        // Per-residue RSCC over union of grid points
        if (!res_obs.empty()) {
            result.by_residue[res] = PearsonCorrelation(res_obs, res_calc);
            all_obs.insert(all_obs.end(), res_obs.begin(), res_obs.end());
            all_calc.insert(all_calc.end(),
                res_calc.begin(), res_calc.end());
        } else {
            result.by_residue[res] =
                std::numeric_limits<double>::quiet_NaN();
        }
    }

    // Overall RSCC
    if (!all_obs.empty()) {
        result.overall = PearsonCorrelation(all_obs, all_calc);
    } else {
        result.overall = std::numeric_limits<double>::quiet_NaN();
    }

    return result;
}

DensityScoreResult RSR(
    OEChem::OEMolBase& mol,
    const OESystem::OEScalarGrid& grid,
    double resolution,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask,
    const OESystem::OEScalarGrid* calc_grid) {
    if (resolution <= 0.0) {
        throw GridError("Resolution must be positive");
    }
    PrepareStructure(mol);

    auto residue_atoms = CollectAtomsByResidue(mol, mask);
    if (residue_atoms.empty()) {
        throw StructureError("No scorable heavy atoms after applying mask");
    }

    if (calc_grid == nullptr) {
        throw GridError("calc_grid is required (auto-generation not yet supported)");
    }

    DensityScoreResult result;
    std::vector<double> all_obs, all_calc;

    for (const auto& [res, atoms] : residue_atoms) {
        std::vector<double> res_obs, res_calc;

        for (const auto* atom : atoms) {
            double x, y, z;
            GetAtomCoords(mol, *atom, x, y, z);

            if (!grid.IsInGrid(static_cast<float>(x),
                               static_cast<float>(y),
                               static_cast<float>(z))) {
                result.by_atom[atom->GetIdx()] =
                    std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            double radius = ScoringRadius(*atom, resolution);
            auto pts = GetAtomGridPoints(grid, x, y, z, radius);
            if (pts.empty()) {
                result.by_atom[atom->GetIdx()] =
                    std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            std::vector<double> obs_vals, calc_vals;
            obs_vals.reserve(pts.size());
            calc_vals.reserve(pts.size());

            for (unsigned int idx : pts) {
                obs_vals.push_back(grid[idx]);
                calc_vals.push_back((*calc_grid)[idx]);
            }

            // Per-atom RSR: sum|obs-calc| / sum|obs+calc|
            double num = 0.0, den = 0.0;
            for (size_t i = 0; i < obs_vals.size(); ++i) {
                num += std::abs(obs_vals[i] - calc_vals[i]);
                den += std::abs(obs_vals[i] + calc_vals[i]);
            }
            result.by_atom[atom->GetIdx()] =
                (den > 0.0) ? (num / den)
                            : std::numeric_limits<double>::quiet_NaN();

            res_obs.insert(res_obs.end(), obs_vals.begin(), obs_vals.end());
            res_calc.insert(res_calc.end(),
                calc_vals.begin(), calc_vals.end());
        }

        // Per-residue RSR over union of grid points
        if (!res_obs.empty()) {
            double num = 0.0, den = 0.0;
            for (size_t i = 0; i < res_obs.size(); ++i) {
                num += std::abs(res_obs[i] - res_calc[i]);
                den += std::abs(res_obs[i] + res_calc[i]);
            }
            result.by_residue[res] =
                (den > 0.0) ? (num / den)
                            : std::numeric_limits<double>::quiet_NaN();
            all_obs.insert(all_obs.end(), res_obs.begin(), res_obs.end());
            all_calc.insert(all_calc.end(),
                res_calc.begin(), res_calc.end());
        } else {
            result.by_residue[res] =
                std::numeric_limits<double>::quiet_NaN();
        }
    }

    // Overall RSR
    if (!all_obs.empty()) {
        double num = 0.0, den = 0.0;
        for (size_t i = 0; i < all_obs.size(); ++i) {
            num += std::abs(all_obs[i] - all_calc[i]);
            den += std::abs(all_obs[i] + all_calc[i]);
        }
        result.overall =
            (den > 0.0) ? (num / den)
                        : std::numeric_limits<double>::quiet_NaN();
    } else {
        result.overall = std::numeric_limits<double>::quiet_NaN();
    }

    return result;
}

DensityScoreResult QScore(
    OEChem::OEMolBase& mol,
    const OESystem::OEScalarGrid& grid,
    double resolution,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask,
    const QScoreOptions& options) {
    if (resolution <= 0.0) {
        throw GridError("Resolution must be positive");
    }
    PrepareStructure(mol);

    auto residue_atoms = CollectAtomsByResidue(mol, mask);
    if (residue_atoms.empty()) {
        throw StructureError("No scorable heavy atoms after applying mask");
    }

    // Map normalization
    double A, B;
    if (options.GetNormalizeMap()) {
        GetMapNormalization(grid, A, B);
    } else {
        A = 1.0;
        B = 0.0;
    }
    double sigma = options.GetSigma();

    // Build spatial index for point isolation
    std::unique_ptr<SpatialIndex> spatial_idx;
    if (options.GetIsolatePoints()) {
        spatial_idx = std::make_unique<SpatialIndex>(mol);
    }

    // Pre-compute unit sphere offsets for fixed mode
    std::unordered_map<int, std::vector<std::array<double, 3>>> unit_spheres;
    if (options.GetRadialSampling() == RadialSampling::Fixed) {
        double R = options.GetRadialStep();
        while (R < options.GetMaxRadius() + 0.01) {
            int rkey = static_cast<int>(std::round(R * 1e6));
            unit_spheres[rkey] = FibonacciSpherePoints(
                0.0, 0.0, 0.0, R, options.GetNumPoints());
            R += options.GetRadialStep();
        }
    }

    // Pre-compute reference Gaussian values per shell
    std::unordered_map<int, double> ref_by_shell;
    if (options.GetRadialSampling() == RadialSampling::Fixed) {
        double R = options.GetRadialStep();
        while (R < options.GetMaxRadius() + 0.01) {
            int rkey = static_cast<int>(std::round(R * 1e6));
            if (options.GetNormalizeMap()) {
                ref_by_shell[rkey] =
                    A * std::exp(-0.5 * (R / sigma) * (R / sigma)) + B;
            } else {
                ref_by_shell[rkey] =
                    std::exp(-(R * R) / (2.0 * sigma * sigma));
            }
            R += options.GetRadialStep();
        }
    }

    constexpr int MIN_SHELLS = 7;
    double grid_spacing = grid.GetSpacing();

    DensityScoreResult result;
    std::vector<double> all_q;

    for (const auto& [res, atoms] : residue_atoms) {
        std::vector<double> res_q;

        for (const auto* atom : atoms) {
            double x, y, z;
            GetAtomCoords(mol, *atom, x, y, z);

            if (!grid.IsInGrid(static_cast<float>(x),
                               static_cast<float>(y),
                               static_cast<float>(z))) {
                result.by_atom[atom->GetIdx()] =
                    std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            // Determine radial parameters
            double step, max_r;
            if (options.GetRadialSampling() == RadialSampling::Adaptive) {
                step = std::min(static_cast<double>(grid_spacing),
                                resolution / MIN_SHELLS);
                max_r = atom->GetRadius() * 2.0;
            } else {
                step = options.GetRadialStep();
                max_r = options.GetMaxRadius();
            }

            // Collect sample points and reference values
            std::vector<double> sample_x, sample_y, sample_z;
            std::vector<double> ref_vals;

            // Center point: replicate num_points times for equal weighting
            double ref_at_center = options.GetNormalizeMap() ? (A + B) : 1.0;
            for (unsigned int p = 0; p < options.GetNumPoints(); ++p) {
                sample_x.push_back(x);
                sample_y.push_back(y);
                sample_z.push_back(z);
                ref_vals.push_back(ref_at_center);
            }

            // Radial shells
            double R = step;
            while (R < max_r + 0.01) {
                int rkey = static_cast<int>(std::round(R * 1e6));

                // Generate sphere points
                std::vector<std::array<double, 3>> shell_pts;
                auto it = unit_spheres.find(rkey);
                if (it != unit_spheres.end()) {
                    // Use pre-computed offsets translated to atom center
                    shell_pts.reserve(it->second.size());
                    for (const auto& off : it->second) {
                        shell_pts.push_back(
                            {x + off[0], y + off[1], z + off[2]});
                    }
                } else {
                    shell_pts = FibonacciSpherePoints(
                        x, y, z, R, options.GetNumPoints());
                }

                // Point isolation: remove points closer to neighbor atoms
                if (options.GetIsolatePoints() && spatial_idx) {
                    double threshold = 0.9 * R;
                    unsigned int target_idx = atom->GetIdx();

                    std::vector<std::array<double, 3>> filtered;
                    for (const auto& pt : shell_pts) {
                        auto neighbors = spatial_idx->FindWithinRadius(
                            pt[0], pt[1], pt[2], threshold);
                        bool has_other = false;
                        for (unsigned int ni : neighbors) {
                            if (ni != target_idx) {
                                has_other = true;
                                break;
                            }
                        }
                        if (!has_other) {
                            filtered.push_back(pt);
                        }
                    }

                    // Retry with more points if too few survived
                    if (filtered.size() < options.GetNumPoints()) {
                        for (int attempt = 1; attempt < 50; ++attempt) {
                            int n_gen = options.GetNumPoints() + attempt * 2;
                            auto retry_pts = FibonacciSpherePoints(
                                x, y, z, R, n_gen);
                            filtered.clear();
                            for (const auto& pt : retry_pts) {
                                auto nbrs = spatial_idx->FindWithinRadius(
                                    pt[0], pt[1], pt[2], threshold);
                                bool other = false;
                                for (unsigned int ni : nbrs) {
                                    if (ni != target_idx) {
                                        other = true;
                                        break;
                                    }
                                }
                                if (!other) {
                                    filtered.push_back(pt);
                                }
                            }
                            if (filtered.size() >= options.GetNumPoints()) break;
                        }
                    }
                    shell_pts = std::move(filtered);
                }

                if (!shell_pts.empty()) {
                    double ref_at_R;
                    auto rit = ref_by_shell.find(rkey);
                    if (rit != ref_by_shell.end()) {
                        ref_at_R = rit->second;
                    } else if (options.GetNormalizeMap()) {
                        ref_at_R = A * std::exp(
                            -0.5 * (R / sigma) * (R / sigma)) + B;
                    } else {
                        ref_at_R = std::exp(
                            -(R * R) / (2.0 * sigma * sigma));
                    }

                    for (const auto& pt : shell_pts) {
                        sample_x.push_back(pt[0]);
                        sample_y.push_back(pt[1]);
                        sample_z.push_back(pt[2]);
                        ref_vals.push_back(ref_at_R);
                    }
                }

                R += step;
            }

            // Interpolate density at all sample points
            std::vector<double> map_vals;
            std::vector<double> map_refs;
            for (size_t i = 0; i < sample_x.size(); ++i) {
                double val = InterpolateDensity(
                    grid, sample_x[i], sample_y[i], sample_z[i],
                    std::numeric_limits<double>::quiet_NaN());
                if (!std::isnan(val)) {
                    map_vals.push_back(val);
                    map_refs.push_back(ref_vals[i]);
                }
            }

            // Pearson correlation between observed and reference profiles
            double q;
            if (map_vals.size() >= 3) {
                q = PearsonCorrelation(map_vals, map_refs);
            } else {
                q = std::numeric_limits<double>::quiet_NaN();
            }

            result.by_atom[atom->GetIdx()] = q;
            if (!std::isnan(q)) {
                res_q.push_back(q);
            }
        }

        // Per-residue: mean of atom Q-scores
        if (!res_q.empty()) {
            result.by_residue[res] = std::accumulate(
                res_q.begin(), res_q.end(), 0.0) / res_q.size();
            all_q.insert(all_q.end(), res_q.begin(), res_q.end());
        } else {
            result.by_residue[res] =
                std::numeric_limits<double>::quiet_NaN();
        }
    }

    // Overall: mean of all atom Q-scores
    if (!all_q.empty()) {
        result.overall = std::accumulate(
            all_q.begin(), all_q.end(), 0.0) / all_q.size();
    } else {
        result.overall = std::numeric_limits<double>::quiet_NaN();
    }

    return result;
}

DensityScoreResult EDIAm(
    OEChem::OEMolBase& mol,
    const OESystem::OEScalarGrid& grid,
    double resolution,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask) {
    if (resolution <= 0.0) {
        throw GridError("Resolution must be positive");
    }
    PrepareStructure(mol);

    auto residue_atoms = CollectAtomsByResidue(mol, mask);
    if (residue_atoms.empty()) {
        throw StructureError("No scorable heavy atoms after applying mask");
    }

    double rho_expected = 1.5 / (resolution * resolution);

    DensityScoreResult result;
    std::vector<double> all_scores;

    for (const auto& [res, atoms] : residue_atoms) {
        std::vector<double> res_scores;

        for (const auto* atom : atoms) {
            double x, y, z;
            GetAtomCoords(mol, *atom, x, y, z);

            if (!grid.IsInGrid(static_cast<float>(x),
                               static_cast<float>(y),
                               static_cast<float>(z))) {
                result.by_atom[atom->GetIdx()] =
                    std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            // Sample at atom center
            std::vector<double> sigmoid_vals;
            double rho = InterpolateDensity(grid, x, y, z,
                                            std::numeric_limits<double>::quiet_NaN());
            if (!std::isnan(rho) && rho_expected > 0.0) {
                sigmoid_vals.push_back(EDIAmSigmoid(rho / rho_expected));
            }

            // Sample at bond midpoints to heavy neighbors
            for (OESystem::OEIter<OEChem::OEBondBase> bond = atom->GetBonds();
                 bond; ++bond) {
                const OEChem::OEAtomBase* nbr = bond->GetNbr(atom);
                if (nbr->GetAtomicNum() == 1) continue;  // skip H

                double nx, ny, nz;
                GetAtomCoords(mol, *nbr, nx, ny, nz);
                double mx = (x + nx) / 2.0;
                double my = (y + ny) / 2.0;
                double mz = (z + nz) / 2.0;

                double mid_rho = InterpolateDensity(grid, mx, my, mz,
                    std::numeric_limits<double>::quiet_NaN());
                if (!std::isnan(mid_rho) && rho_expected > 0.0) {
                    sigmoid_vals.push_back(
                        EDIAmSigmoid(mid_rho / rho_expected));
                }
            }

            double score;
            if (!sigmoid_vals.empty()) {
                score = std::accumulate(sigmoid_vals.begin(),
                    sigmoid_vals.end(), 0.0) / sigmoid_vals.size();
            } else {
                score = std::numeric_limits<double>::quiet_NaN();
            }

            result.by_atom[atom->GetIdx()] = score;
            if (!std::isnan(score)) {
                res_scores.push_back(score);
            }
        }

        if (!res_scores.empty()) {
            result.by_residue[res] = std::accumulate(res_scores.begin(),
                res_scores.end(), 0.0) / res_scores.size();
            all_scores.insert(all_scores.end(),
                res_scores.begin(), res_scores.end());
        } else {
            result.by_residue[res] =
                std::numeric_limits<double>::quiet_NaN();
        }
    }

    if (!all_scores.empty()) {
        result.overall = std::accumulate(all_scores.begin(),
            all_scores.end(), 0.0) / all_scores.size();
    } else {
        result.overall = std::numeric_limits<double>::quiet_NaN();
    }

    return result;
}

DensityScoreResult Coverage(
    OEChem::OEMolBase& mol,
    const OESystem::OEScalarGrid& grid,
    double sigma,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask) {
    PrepareStructure(mol);
    auto residue_atoms = CollectAtomsByResidue(mol, mask);
    if (residue_atoms.empty()) {
        throw StructureError("No scorable heavy atoms after applying mask");
    }

    MapStats stats = ComputeMapStats(grid);
    double threshold = stats.mean + sigma * stats.stddev;

    DensityScoreResult result;
    std::vector<double> all_scores;

    for (const auto& [res, atoms] : residue_atoms) {
        std::vector<double> res_scores;

        for (const auto* atom : atoms) {
            double x, y, z;
            GetAtomCoords(mol, *atom, x, y, z);

            if (!grid.IsInGrid(static_cast<float>(x),
                               static_cast<float>(y),
                               static_cast<float>(z))) {
                result.by_atom[atom->GetIdx()] =
                    std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            double rho = InterpolateDensity(grid, x, y, z);
            if (std::isnan(rho)) {
                result.by_atom[atom->GetIdx()] =
                    std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            double score = (rho >= threshold) ? 1.0 : 0.0;
            result.by_atom[atom->GetIdx()] = score;
            res_scores.push_back(score);
        }

        if (!res_scores.empty()) {
            double res_mean = std::accumulate(res_scores.begin(),
                res_scores.end(), 0.0) / res_scores.size();
            result.by_residue[res] = res_mean;
            all_scores.insert(all_scores.end(),
                res_scores.begin(), res_scores.end());
        } else {
            result.by_residue[res] =
                std::numeric_limits<double>::quiet_NaN();
        }
    }

    if (!all_scores.empty()) {
        result.overall = std::accumulate(all_scores.begin(),
            all_scores.end(), 0.0) / all_scores.size();
    } else {
        result.overall = std::numeric_limits<double>::quiet_NaN();
    }

    return result;
}

}  // namespace Maptitude
