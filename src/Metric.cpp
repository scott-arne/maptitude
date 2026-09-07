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
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace Maptitude {

// Import detail helpers into this TU
using detail::MapStats;
using detail::pearson_correlation;
using detail::ediam_sigmoid;
using detail::fibonacci_sphere_points;

// ---- Radial sweep validation ----

/// Explain why a radial step cannot drive the shell loop, or return an empty string when it
/// can.
///
/// Split out from the sweep check because the step is loop-invariant in both modes: FIXED
/// reads it from the options and ADAPTIVE derives it from the resolution argument and the
/// grid spacing. Neither depends on the atom, so an unusable step is always an error about
/// the call and is reported before any per-atom work begins.
static std::string DescribeUnusableStep(const char* label, double step) {
    if (!std::isfinite(step) || step < MIN_RADIAL_STEP) {
        std::ostringstream message;
        message << label << " radial sweep needs a step of at least " << MIN_RADIAL_STEP
                << " A (got " << step << "); below that the shell accumulator does not advance";
        return message.str();
    }
    return {};
}

/// Explain why a maximum radius on its own cannot drive a radial sweep, or return an empty
/// string when it can.
///
/// Split from the shell-count check because these two bounds read the radius alone. In
/// ADAPTIVE mode the radius is twice the atom's own, so a failure here can differ between
/// atoms of the same call whatever the step is, and belongs to the atom.
///
/// The lower and upper bounds are reported separately because they fail for unrelated
/// reasons: the lower bound is reached by an atom carrying no radius -- a fact about the
/// molecule, which the int-overflow ceiling does not describe.
static std::string DescribeUnusableRadius(const char* label, bool adaptive, double max_r) {
    std::ostringstream message;
    if (!std::isfinite(max_r) || max_r <= 0.0) {
        message << label << " radial sweep needs a positive maximum radius (got " << max_r << ")";
        if (adaptive) {
            message << "; the adaptive radius is twice the atom's own, so this atom has none assigned";
        }
        return message.str();
    }
    if (max_r > MAX_RADIUS_LIMIT) {
        message << label << " radial sweep needs a maximum radius of at most " << MAX_RADIUS_LIMIT
                << " A (got " << max_r << "); beyond that the shell key overflows int";
        return message.str();
    }
    return {};
}

/// Explain why a step and a maximum radius cannot produce a usable shell count together, or
/// return an empty string when they can.
///
/// Split from the radius check because every bound here reads both quantities, so in
/// ADAPTIVE mode a failure is attributable to neither on its own: it says the step this call
/// derived does not fit this atom's radius. `qscore` resolves that by quantifying over the
/// molecule -- see `RequireSomeAtomCanSweep`.
///
/// Expects a step and a radius that have already passed their own checks; the shell count
/// below is only meaningful once both are finite and positive.
static std::string DescribeUnusableShellCount(const char* label, double step, double max_r,
                                              unsigned int num_points) {
    std::ostringstream message;
    if (step >= max_r + 0.01) {
        message << label << " radial sweep produces no shells: step " << step
                << " A is not smaller than the maximum radius " << max_r << " A";
        return message.str();
    }
    const double shells = (max_r + 0.01) / step;
    if (shells > static_cast<double>(MAX_SHELLS)) {
        message << label << " radial sweep would run " << static_cast<long long>(shells)
                << " shells, over the " << MAX_SHELLS
                << " limit; raise the step or lower the maximum radius";
        return message.str();
    }
    // Shells and points are each bounded on their own, but it is their product that
    // allocates: FIXED precomputes one num_points-element sphere per shell, and both
    // modes push one sample per point per shell into four parallel vectors. Bounding
    // only the factors admits 510000 shells of 10000 points -- over 100 GB before any
    // scoring happens.
    if (shells * static_cast<double>(num_points) > static_cast<double>(MAX_TOTAL_SAMPLES)) {
        message << label << " radial sweep would take "
                << static_cast<long long>(shells * static_cast<double>(num_points))
                << " samples per atom (" << static_cast<long long>(shells) << " shells x "
                << num_points << " points), over the " << MAX_TOTAL_SAMPLES
                << " limit; raise the step, lower the maximum radius, or use fewer points";
        return message.str();
    }
    return {};
}

/// Explain why a radial sweep cannot terminate or cannot produce a score, or return an
/// empty string when it can.
///
/// Both sampling modes converge on the same shell loop, but only FIXED mode's step and
/// radius arrive through QScoreOptions' validated setters. ADAPTIVE derives its own from the
/// grid spacing, the resolution, and the atom, and consults no option, so the same invariants
/// have to be re-established here.
///
/// This is the whole check, in the order the three parts fail. FIXED evaluates it once, as
/// one loop-invariant statement about the call; ADAPTIVE evaluates the parts separately,
/// because only the radius part is a statement about a single atom.
static std::string DescribeUnusableSweep(RadialSampling mode, double step, double max_r,
                                         unsigned int num_points) {
    const bool adaptive = (mode == RadialSampling::ADAPTIVE);
    const char* label = adaptive ? "Adaptive" : "Fixed";
    const std::string step_problem = DescribeUnusableStep(label, step);
    if (!step_problem.empty()) {
        return step_problem;
    }
    const std::string radius_problem = DescribeUnusableRadius(label, adaptive, max_r);
    if (!radius_problem.empty()) {
        return radius_problem;
    }
    return DescribeUnusableShellCount(label, step, max_r, num_points);
}

/// Reject a radial sweep that cannot terminate or cannot produce a score.
static void RequireUsableSweep(RadialSampling mode, double step, double max_r,
                               unsigned int num_points) {
    const std::string reason = DescribeUnusableSweep(mode, step, max_r, num_points);
    if (!reason.empty()) {
        throw GridError(reason);
    }
}

// ---- OE-aware wrappers around detail helpers ----

static MapStats compute_map_stats(const OESystem::OESkewGrid& grid) {
    const unsigned int size = grid.GetSize();
    if (size == 0) return {};
    const float* raw = grid.GetValues();
    std::vector<double> values(size);
    for (unsigned int i = 0; i < size; ++i) {
        values[i] = raw[i];
    }
    return detail::compute_map_stats(values);
}

static void get_map_normalization(const OESystem::OESkewGrid& grid,
                                  double& A, double& B) {
    const unsigned int size = grid.GetSize();
    if (size == 0) { A = 1.0; B = 0.0; return; }
    const float* raw = grid.GetValues();
    std::vector<double> values(size);
    for (unsigned int i = 0; i < size; ++i) {
        values[i] = raw[i];
    }
    detail::get_map_normalization(values.data(), values.size(), A, B);
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

        const Residue res = Residue::FromAtom(*atom);
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

static double scoring_radius(const OEChem::OEAtomBase& atom, const double resolution) {
    const OEChem::OEResidue res = OEChem::OEAtomGetResidue(&atom);
    return detail::scoring_radius(res.GetBFactor(), resolution);
}

// ---- Helper: assign Bondi VDW radii if missing ----

void detail::assign_missing_radii(OEChem::OEMolBase& mol) {
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(OEChem::OEIsHeavy()); atom; ++atom) {
        if (atom->GetRadius() > 0.0) return;  // radii already assigned
    }
    OEChem::OEAssignBondiVdWRadii(mol);
}

/// Reject `AtomRadius::ADAPTIVE` for RSCC.
///
/// RSCC has no adaptive radius model and never had one. Its radius switch listed FIXED,
/// SCALED, and `BINNED: default:`, so ADAPTIVE fell through to the binned radius and
/// returned the binned score under another name -- bit-identical to it, with no exception,
/// no warning, and no field on the result recording which model ran. Giving RSCC a real
/// adaptive radius is an accuracy change; until then an exception is preferable to a
/// plausible wrong answer. `std::invalid_argument` surfaces in Python as `RuntimeError`,
/// the documented channel for option-value errors.
[[noreturn]] static void RejectAdaptiveRsccRadius() {
    throw std::invalid_argument(
        "rscc has no adaptive atom-radius model; use AtomRadius::FIXED, AtomRadius::SCALED, "
        "or AtomRadius::BINNED. AtomRadius::ADAPTIVE is supported by rsr only");
}

/// Reject an `AtomRadius` value the enum does not declare, at the point of use.
///
/// Unreachable through the public API: both option classes validate in
/// `SetAtomRadiusMethod` and hold the value privately, so no caller can present an
/// undeclared one here. It exists to make the two selectors below total functions --
/// every path returns a radius or throws -- instead of falling out of a switch with the
/// result still uninitialized. That shape is what made `static_cast<AtomRadius>(42)` an
/// indeterminate read, and a `default:` label cannot replace it without costing the
/// -Wswitch warning these switches are written without one to get.
[[noreturn]] static void RejectUndeclaredAtomRadius(const char* metric, AtomRadius method) {
    std::ostringstream message;
    message << metric << " received an AtomRadius value the enum does not declare ("
            << static_cast<int>(method) << "); this is a bug in the caller's construction of the "
            << "options object, which validates in SetAtomRadiusMethod";
    throw std::invalid_argument(message.str());
}

/// Select the RSCC scoring radius for one atom.
static double rscc_atom_radius(const RsccOptions& options, const OEChem::OEAtomBase& atom,
                               const double resolution) {
    switch (options.GetAtomRadiusMethod()) {
        case AtomRadius::FIXED:
            return options.GetFixedAtomRadius();
        case AtomRadius::SCALED: {
            const double scaled = atom.GetRadius() * options.GetAtomRadiusScaling();
            return (scaled < 0.1) ? 1.5 : scaled;
        }
        case AtomRadius::BINNED:
            return detail::binned_atom_radius(resolution);
        case AtomRadius::ADAPTIVE:
            // Unreachable: rejected before the loop. Listed anyway so the switch is
            // exhaustive without a `default:` label, which is what makes a future
            // enumerator a compiler warning here instead of another silent
            // fall-through into whichever model happens to be last.
            RejectAdaptiveRsccRadius();
    }
    RejectUndeclaredAtomRadius("rscc", options.GetAtomRadiusMethod());
}

/// Select the RSR scoring radius for one atom.
static double rsr_atom_radius(const RsrOptions& options, const OEChem::OEAtomBase& atom,
                              const double resolution) {
    switch (options.GetAtomRadiusMethod()) {
        case AtomRadius::FIXED:
            return options.GetFixedAtomRadius();
        case AtomRadius::SCALED: {
            const double scaled = atom.GetRadius() * options.GetAtomRadiusScaling();
            return (scaled < 0.1) ? 1.5 : scaled;
        }
        case AtomRadius::BINNED:
            return detail::binned_atom_radius(resolution);
        case AtomRadius::ADAPTIVE:
            return scoring_radius(atom, resolution);
    }
    RejectUndeclaredAtomRadius("rsr", options.GetAtomRadiusMethod());
}

// ==== Density scoring functions ====

DensityScoreResult rscc(
    OEChem::OEMolBase& mol,
    const OESystem::OESkewGrid& grid,
    const double resolution,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask,
    const OESystem::OESkewGrid* calc_grid,
    const RsccOptions& options) {
    require_usable_resolution(resolution);
    if (options.GetAtomRadiusMethod() == AtomRadius::ADAPTIVE) {
        RejectAdaptiveRsccRadius();
    }
    detail::assign_missing_radii(mol);

    auto residue_atoms = CollectAtomsByResidue(mol, mask);
    if (residue_atoms.empty()) {
        throw StructureError("No scorable heavy atoms after applying mask");
    }

    // If no calc_grid provided, we would compute one via DensityCalculator.
    // For now, require calc_grid to be provided.
    if (calc_grid == nullptr) {
        throw GridError("calc_grid is required (auto-generation not yet supported)");
    }

    const GridParams gp = get_grid_params(grid);
    const float* obs_values = grid.GetValues();
    const float* calc_values = calc_grid->GetValues();

    DensityScoreResult result;
    std::vector<double> all_obs, all_calc;

    for (const auto& [res, atoms] : residue_atoms) {
        std::vector<double> res_obs, res_calc;

        for (const auto* atom : atoms) {
            double x, y, z;
            GetAtomCoords(mol, *atom, x, y, z);

            if (!grid_contains(gp, x, y, z)) {
                result.by_atom[atom->GetIdx()] =
                    std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            const double radius = rscc_atom_radius(options, *atom, resolution);

            auto pts = get_atom_grid_points(grid, x, y, z, radius);
            if (pts.empty()) {
                result.by_atom[atom->GetIdx()] =
                    std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            std::vector<double> obs_vals, calc_vals;
            obs_vals.reserve(pts.size());
            calc_vals.reserve(pts.size());

            for (unsigned int idx : pts) {
                obs_vals.push_back(obs_values[idx]);
                calc_vals.push_back(calc_values[idx]);
            }

            // Per-atom RSCC
            result.by_atom[atom->GetIdx()] =
                pearson_correlation(obs_vals, calc_vals);

            res_obs.insert(res_obs.end(), obs_vals.begin(), obs_vals.end());
            res_calc.insert(res_calc.end(),
                calc_vals.begin(), calc_vals.end());
        }

        // Per-residue RSCC over union of grid points
        if (!res_obs.empty()) {
            result.by_residue[res] = pearson_correlation(res_obs, res_calc);
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
        result.overall = pearson_correlation(all_obs, all_calc);
    } else {
        result.overall = std::numeric_limits<double>::quiet_NaN();
    }

    return result;
}

DensityScoreResult rsr(
    OEChem::OEMolBase& mol,
    const OESystem::OESkewGrid& grid,
    double resolution,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask,
    const OESystem::OESkewGrid* calc_grid,
    const RsrOptions& options) {
    require_usable_resolution(resolution);
    detail::assign_missing_radii(mol);

    auto residue_atoms = CollectAtomsByResidue(mol, mask);
    if (residue_atoms.empty()) {
        throw StructureError("No scorable heavy atoms after applying mask");
    }

    if (calc_grid == nullptr) {
        throw GridError("calc_grid is required (auto-generation not yet supported)");
    }

    const GridParams gp = get_grid_params(grid);
    const float* obs_values = grid.GetValues();
    const float* calc_values = calc_grid->GetValues();

    DensityScoreResult result;
    std::vector<double> all_obs, all_calc;

    for (const auto& [res, atoms] : residue_atoms) {
        std::vector<double> res_obs, res_calc;

        for (const auto* atom : atoms) {
            double x, y, z;
            GetAtomCoords(mol, *atom, x, y, z);

            if (!grid_contains(gp, x, y, z)) {
                result.by_atom[atom->GetIdx()] =
                    std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            const double radius = rsr_atom_radius(options, *atom, resolution);

            auto pts = get_atom_grid_points(grid, x, y, z, radius);
            if (pts.empty()) {
                result.by_atom[atom->GetIdx()] =
                    std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            std::vector<double> obs_vals, calc_vals;
            obs_vals.reserve(pts.size());
            calc_vals.reserve(pts.size());

            for (unsigned int idx : pts) {
                obs_vals.push_back(obs_values[idx]);
                calc_vals.push_back(calc_values[idx]);
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

/// Reject an adaptive step that no atom in this molecule can sweep with.
///
/// The shell-count bounds read the step and the maximum radius jointly, so in ADAPTIVE mode
/// neither is on its own responsible for a failure. Which of the two it belongs to is
/// settled by quantifying over the molecule: if some atom can be scored, an atom that cannot
/// differs from it only in its own radius and is scored NaN in the loop below, as an atom
/// with no radius at all is; if no atom can, the statement no longer mentions any particular
/// atom and is a fact about the step, which comes from the resolution argument and the grid
/// spacing. A caller error must throw wherever it is evaluated, so that case throws here.
///
/// This matters most for the branch that fails when the radius is too small for the step.
/// At a grid spacing of 4 A no atom in the periodic table has a radius large enough, so
/// leaving it per-atom turned a plainly unusable resolution-and-spacing pair into a result
/// object full of NaN. It is genuinely atom-dependent in form -- a larger atom would pass --
/// which is why the quantifier decides it rather than a classification fixed in advance.
///
/// Atoms that cannot be scored for reasons of their own are skipped rather than counted as
/// failures: an out-of-grid atom is NaN before the sweep is consulted, and an atom with no
/// radius fails the radius check whatever the step is. Neither supports a conclusion about
/// the step, so a molecule of nothing but those still returns per-atom NaN.
static void RequireSomeAtomCanSweep(
    const OEChem::OEMolBase& mol, const OESystem::OESkewGrid& grid,
    const std::map<Residue, std::vector<const OEChem::OEAtomBase*>>& residue_atoms, double step,
    unsigned int num_points) {
    const GridParams gp = get_grid_params(grid);
    std::string first_failure;
    for (const auto& [res, atoms] : residue_atoms) {
        for (const auto* atom : atoms) {
            double x, y, z;
            GetAtomCoords(mol, *atom, x, y, z);
            if (!grid_contains(gp, x, y, z)) {
                continue;
            }
            const double max_r = atom->GetRadius() * 2.0;
            if (!DescribeUnusableRadius("Adaptive", true, max_r).empty()) {
                continue;
            }
            const std::string shell_problem =
                DescribeUnusableShellCount("Adaptive", step, max_r, num_points);
            if (shell_problem.empty()) {
                return;
            }
            if (first_failure.empty()) {
                first_failure = shell_problem;
            }
        }
    }
    if (!first_failure.empty()) {
        throw GridError(first_failure +
                        "; no atom in this molecule can be swept at this step, which comes from "
                        "the resolution and the grid spacing rather than from any atom");
    }
}

DensityScoreResult qscore(
    OEChem::OEMolBase& mol,
    const OESystem::OESkewGrid& grid,
    double resolution,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask,
    const QScoreOptions& options) {
    require_usable_resolution(resolution);
    detail::assign_missing_radii(mol);

    auto residue_atoms = CollectAtomsByResidue(mol, mask);
    if (residue_atoms.empty()) {
        throw StructureError("No scorable heavy atoms after applying mask");
    }

    // Map normalization
    double A, B;
    if (options.GetNormalizeMap()) {
        get_map_normalization(grid, A, B);
    } else {
        A = 1.0;
        B = 0.0;
    }
    const double sigma = options.GetSigma();

    // Build spatial index for point isolation
    std::unique_ptr<SpatialIndex> spatial_idx;
    if (options.GetIsolatePoints()) {
        spatial_idx = std::make_unique<SpatialIndex>(mol);
    }

    constexpr int MIN_SHELLS = 7;
    const GridParams gp = get_grid_params(grid);
    // The finest sampled direction sets the shell interval: a step taken from a
    // coarser axis would skip nodes along the fine one.
    const double grid_spacing = std::min({gp.x_spacing, gp.y_spacing, gp.z_spacing});

    // Reject every failure the call is responsible for here, before any per-atom work.
    // FIXED's whole sweep qualifies, and the two precompute loops below would otherwise
    // hang on an unusable one. ADAPTIVE's step qualifies on its own -- it comes from the
    // resolution argument and the grid, not from any atom -- and its shell count qualifies
    // when no atom in the molecule can satisfy it. Only the maximum radius is left to the
    // loop, where it is a statement about one atom.
    if (options.GetRadialSampling() == RadialSampling::FIXED) {
        RequireUsableSweep(RadialSampling::FIXED, options.GetRadialStep(), options.GetMaxRadius(),
                           options.GetNumPoints());
    } else {
        const double step = std::min(grid_spacing, resolution / MIN_SHELLS);
        const std::string reason = DescribeUnusableStep("Adaptive", step);
        if (!reason.empty()) {
            throw GridError(reason);
        }
        RequireSomeAtomCanSweep(mol, grid, residue_atoms, step, options.GetNumPoints());
    }

    // Pre-compute unit sphere offsets for fixed mode
    std::unordered_map<int, std::vector<std::array<double, 3>>> unit_spheres;
    if (options.GetRadialSampling() == RadialSampling::FIXED) {
        double R = options.GetRadialStep();
        while (R < options.GetMaxRadius() + 0.01) {
            int rkey = static_cast<int>(std::round(R * 1e6));
            unit_spheres[rkey] = fibonacci_sphere_points(
                0.0, 0.0, 0.0, R, static_cast<int>(options.GetNumPoints()));
            R += options.GetRadialStep();
        }
    }

    // Pre-compute reference Gaussian values per shell
    std::unordered_map<int, double> ref_by_shell;
    if (options.GetRadialSampling() == RadialSampling::FIXED) {
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

    const float* values = grid.GetValues();

    DensityScoreResult result;
    std::vector<double> all_q;

    for (const auto& [res, atoms] : residue_atoms) {
        std::vector<double> res_q;

        for (const auto* atom : atoms) {
            double x, y, z;
            GetAtomCoords(mol, *atom, x, y, z);

            if (!grid_contains(gp, x, y, z)) {
                result.by_atom[atom->GetIdx()] =
                    std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            // Determine radial parameters. This test keys on FIXED, as the pre-loop
            // validation above does, so a value takes matching arms in both. When they
            // disagreed, an enum value the type can hold but does not declare took the
            // `else` of each: validated as adaptive, then executed as fixed, so the fixed
            // sweep's parameters reached the loop unchecked. The setter rejects those
            // values now; keeping the two tests in the same form means the pairing does
            // not depend on it.
            double step, max_r;
            if (options.GetRadialSampling() == RadialSampling::FIXED) {
                // No per-atom check: these parameters are loop-invariant and were
                // validated before the loop.
                step = options.GetRadialStep();
                max_r = options.GetMaxRadius();
            } else {
                step = std::min(grid_spacing, resolution / MIN_SHELLS);
                max_r = atom->GetRadius() * 2.0;
                // The adaptive radius comes from the atom, so an unusable sweep is a fact
                // about this atom and not about the request. Score it NaN and carry on, the
                // way an out-of-grid atom is handled above; throwing here would discard the
                // scores of every other atom in the molecule.
                if (!DescribeUnusableSweep(RadialSampling::ADAPTIVE, step, max_r,
                                           options.GetNumPoints()).empty()) {
                    result.by_atom[atom->GetIdx()] =
                        std::numeric_limits<double>::quiet_NaN();
                    continue;
                }
            }

            // Collect sample points and reference values
            std::vector<double> sample_x, sample_y, sample_z;
            std::vector<double> ref_vals;

            // Center point: replicate num_points times for equal weighting
            const double ref_at_center = options.GetNormalizeMap() ? (A + B) : 1.0;
            for (unsigned int p = 0; p < options.GetNumPoints(); ++p) {
                sample_x.push_back(x);
                sample_y.push_back(y);
                sample_z.push_back(z);
                ref_vals.push_back(ref_at_center);
            }

            // Radial shells
            double R = step;
            while (R < max_r + 0.01) {
                const int rkey = static_cast<int>(std::round(R * 1e6));

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
                    shell_pts = fibonacci_sphere_points(
                        x, y, z, R, static_cast<int>(options.GetNumPoints()));
                }

                // Point isolation: remove points closer to neighbor atoms
                if (options.GetIsolatePoints() && spatial_idx) {
                    const double threshold = 0.9 * R;
                    const unsigned int target_idx = atom->GetIdx();

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
                    if (filtered.size() < static_cast<size_t>(options.GetNumPoints())) {
                        for (int attempt = 1; attempt < 50; ++attempt) {
                            const int n_gen = static_cast<int>(options.GetNumPoints()) + attempt * 2;
                            auto retry_pts = fibonacci_sphere_points(
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
                            if (filtered.size() >= static_cast<size_t>(options.GetNumPoints())) break;
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
                const double val = interpolate_density_at(
                    gp, values, sample_x[i], sample_y[i], sample_z[i],
                    std::numeric_limits<double>::quiet_NaN());
                if (!std::isnan(val)) {
                    map_vals.push_back(val);
                    map_refs.push_back(ref_vals[i]);
                }
            }

            // Pearson correlation between observed and reference profiles
            double q;
            if (map_vals.size() >= 3) {
                q = pearson_correlation(map_vals, map_refs);
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
                res_q.begin(), res_q.end(), 0.0) / static_cast<double>(res_q.size());
            all_q.insert(all_q.end(), res_q.begin(), res_q.end());
        } else {
            result.by_residue[res] =
                std::numeric_limits<double>::quiet_NaN();
        }
    }

    // Overall: mean of all atom Q-scores
    if (!all_q.empty()) {
        result.overall = std::accumulate(
            all_q.begin(), all_q.end(), 0.0) / static_cast<double>(all_q.size());
    } else {
        result.overall = std::numeric_limits<double>::quiet_NaN();
    }

    return result;
}

DensityScoreResult ediam(
    OEChem::OEMolBase& mol,
    const OESystem::OESkewGrid& grid,
    const double resolution,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask) {
    require_usable_resolution(resolution);
    detail::assign_missing_radii(mol);

    auto residue_atoms = CollectAtomsByResidue(mol, mask);
    if (residue_atoms.empty()) {
        throw StructureError("No scorable heavy atoms after applying mask");
    }

    const double rho_expected = 1.5 / (resolution * resolution);

    const GridParams gp = get_grid_params(grid);
    const float* values = grid.GetValues();

    DensityScoreResult result;
    std::vector<double> all_scores;

    for (const auto& [res, atoms] : residue_atoms) {
        std::vector<double> res_scores;

        for (const auto* atom : atoms) {
            double x, y, z;
            GetAtomCoords(mol, *atom, x, y, z);

            if (!grid_contains(gp, x, y, z)) {
                result.by_atom[atom->GetIdx()] =
                    std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            // Sample at atom center
            std::vector<double> sigmoid_vals;
            const double rho = interpolate_density_at(gp, values, x, y, z,
                                                  std::numeric_limits<double>::quiet_NaN());
            if (!std::isnan(rho) && rho_expected > 0.0) {
                sigmoid_vals.push_back(ediam_sigmoid(rho / rho_expected));
            }

            // Sample at bond midpoints to heavy neighbors
            for (OESystem::OEIter<OEChem::OEBondBase> bond = atom->GetBonds();
                 bond; ++bond) {
                const OEChem::OEAtomBase* nbr = bond->GetNbr(atom);
                if (nbr->GetAtomicNum() == 1) continue;  // skip H

                double nx, ny, nz;
                GetAtomCoords(mol, *nbr, nx, ny, nz);
                const double mx = (x + nx) / 2.0;
                const double my = (y + ny) / 2.0;
                const double mz = (z + nz) / 2.0;

                const double mid_rho = interpolate_density_at(gp, values, mx, my, mz,
                    std::numeric_limits<double>::quiet_NaN());
                if (!std::isnan(mid_rho) && rho_expected > 0.0) {
                    sigmoid_vals.push_back(
                        ediam_sigmoid(mid_rho / rho_expected));
                }
            }

            double score;
            if (!sigmoid_vals.empty()) {
                score = std::accumulate(sigmoid_vals.begin(),
                    sigmoid_vals.end(), 0.0) / static_cast<double>(sigmoid_vals.size());
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
                res_scores.end(), 0.0) / static_cast<double>(res_scores.size());
            all_scores.insert(all_scores.end(),
                res_scores.begin(), res_scores.end());
        } else {
            result.by_residue[res] =
                std::numeric_limits<double>::quiet_NaN();
        }
    }

    if (!all_scores.empty()) {
        result.overall = std::accumulate(all_scores.begin(),
            all_scores.end(), 0.0) / static_cast<double>(all_scores.size());
    } else {
        result.overall = std::numeric_limits<double>::quiet_NaN();
    }

    return result;
}

DensityScoreResult coverage(
    OEChem::OEMolBase& mol,
    const OESystem::OESkewGrid& grid,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask,
    const CoverageOptions& options) {
    detail::assign_missing_radii(mol);
    auto residue_atoms = CollectAtomsByResidue(mol, mask);
    if (residue_atoms.empty()) {
        throw StructureError("No scorable heavy atoms after applying mask");
    }

    const double sigma = options.GetSigma();
    const MapStats stats = compute_map_stats(grid);
    const double threshold = stats.mean + sigma * stats.stddev;

    const GridParams gp = get_grid_params(grid);
    const float* values = grid.GetValues();

    DensityScoreResult result;
    std::vector<double> all_scores;

    for (const auto& [res, atoms] : residue_atoms) {
        std::vector<double> res_scores;

        for (const auto* atom : atoms) {
            double x, y, z;
            GetAtomCoords(mol, *atom, x, y, z);

            if (!grid_contains(gp, x, y, z)) {
                result.by_atom[atom->GetIdx()] =
                    std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            const double rho = interpolate_density_at(gp, values, x, y, z);
            if (std::isnan(rho)) {
                result.by_atom[atom->GetIdx()] =
                    std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            const double score = (rho >= threshold) ? 1.0 : 0.0;
            result.by_atom[atom->GetIdx()] = score;
            res_scores.push_back(score);
        }

        if (!res_scores.empty()) {
            const double res_mean = std::accumulate(res_scores.begin(),
                res_scores.end(), 0.0) / static_cast<double>(res_scores.size());
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
            all_scores.end(), 0.0) / static_cast<double>(all_scores.size());
    } else {
        result.overall = std::numeric_limits<double>::quiet_NaN();
    }

    return result;
}

}  // namespace Maptitude
