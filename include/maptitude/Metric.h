/**
 * @file Metric.h
 * @brief Real-space density scoring metrics.
 *
 * Provides functions for computing RSCC, RSR, Q-score, EDIAm,
 * and Coverage metrics that assess how well a molecular model fits
 * an experimental electron density map.
 */

#ifndef MAPTITUDE_METRIC_H
#define MAPTITUDE_METRIC_H

#include "maptitude/DensityScoreResult.h"
#include "maptitude/QScoreOptions.h"
#include "maptitude/RSCCOptions.h"

namespace OEChem {
class OEMolBase;
class OEAtomBase;
}

namespace OESystem {
class OEScalarGrid;
}

namespace OESystem {
template <class T> class OEUnaryPredicate;
}

namespace Maptitude {

/**
 * @brief Real-Space Correlation Coefficient.
 *
 * Pearson correlation between observed and calculated density at
 * grid points within a scaled vdW radius of each atom.
 *
 * @param mol Input molecule.
 * @param grid Observed electron density map.
 * @param resolution Resolution in Angstroms.
 * @param mask Optional atom predicate.
 * @param calc_grid Optional pre-computed calculated density.
 * @return DensityScoreResult with RSCC values.
 */
DensityScoreResult RSCC(
    OEChem::OEMolBase& mol,
    const OESystem::OEScalarGrid& grid,
    double resolution,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask = nullptr,
    const OESystem::OEScalarGrid* calc_grid = nullptr,
    const RSCCOptions& options = RSCCOptions());

/**
 * @brief Real-Space R-Factor.
 *
 * R-factor between observed and calculated density at grid points
 * within an adaptive radius of each atom.
 *
 * @param mol Input molecule.
 * @param grid Observed electron density map.
 * @param resolution Resolution in Angstroms.
 * @param mask Optional atom predicate.
 * @param calc_grid Optional pre-computed calculated density.
 * @return DensityScoreResult with RSR values.
 */
DensityScoreResult RSR(
    OEChem::OEMolBase& mol,
    const OESystem::OEScalarGrid& grid,
    double resolution,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask = nullptr,
    const OESystem::OEScalarGrid* calc_grid = nullptr);

/**
 * @brief Q-Score (Pintilie et al., 2020).
 *
 * Compares the observed radial density profile around each atom
 * to a reference Gaussian model.
 *
 * @param mol Input molecule.
 * @param grid Observed electron density map.
 * @param resolution Resolution in Angstroms.
 * @param mask Optional atom predicate.
 * @param options Q-score configuration.
 * @return DensityScoreResult with Q-score values.
 */
DensityScoreResult QScore(
    OEChem::OEMolBase& mol,
    const OESystem::OEScalarGrid& grid,
    double resolution,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask = nullptr,
    const QScoreOptions& options = QScoreOptions());

/**
 * @brief Electron Density Index Averaged, Modified (EDIAm).
 *
 * Scores atoms based on density at atom centers and bond midpoints,
 * normalized by resolution-dependent expected density.
 *
 * @param mol Input molecule.
 * @param grid Observed electron density map.
 * @param resolution Resolution in Angstroms.
 * @param mask Optional atom predicate.
 * @return DensityScoreResult with EDIAm values in [0, 1].
 */
DensityScoreResult EDIAm(
    OEChem::OEMolBase& mol,
    const OESystem::OEScalarGrid& grid,
    double resolution,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask = nullptr);

/**
 * @brief Coverage: fraction of atoms observed in density.
 *
 * Binary metric: an atom has coverage 1 if density at its center
 * exceeds (mean + sigma * std) of the map.
 *
 * @param mol Input molecule.
 * @param grid Observed electron density map.
 * @param sigma Number of standard deviations above mean (default: 1.0).
 * @param mask Optional atom predicate.
 * @return DensityScoreResult with coverage fractions.
 */
DensityScoreResult Coverage(
    OEChem::OEMolBase& mol,
    const OESystem::OEScalarGrid& grid,
    double sigma = 1.0,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask = nullptr);

}  // namespace Maptitude

#endif  // MAPTITUDE_METRIC_H
