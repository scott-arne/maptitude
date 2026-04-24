/**
 * @file maptitude.h
 * @brief Main umbrella header for the Maptitude library.
 *
 * Include this header to access all Maptitude functionality for
 * crystallographic electron density computation and real-space scoring.
 *
 * @code
 * #include <maptitude/maptitude.h>
 *
 * // Compute model density
 * Maptitude::UnitCell cell(a, b, c, alpha, beta, gamma);
 * Maptitude::DensityCalculator calc(cell, symops);
 * auto grid = calc.Calculate(mol, obs_grid, resolution);
 *
 * // Score fit to density
 * auto result = Maptitude::rscc(mol, obs_grid, resolution);
 * @endcode
 */

#ifndef MAPTITUDE_MAPTITUDE_H
#define MAPTITUDE_MAPTITUDE_H

/** @brief Major version number */
#define MAPTITUDE_VERSION_MAJOR 0
/** @brief Minor version number */
#define MAPTITUDE_VERSION_MINOR 2
/** @brief Patch version number */
#define MAPTITUDE_VERSION_PATCH 0

/**
 * @namespace Maptitude
 * @brief Maptitude library namespace containing all density computation
 *        and scoring classes.
 */
namespace Maptitude {

// Forward declarations
class DensityCalculator;
class SpatialIndex;
struct Residue;
struct UnitCell;
struct SymOp;
struct DensityScoreResult;
class QScoreOptions;
class RsccOptions;
class RsrOptions;
class CoverageOptions;

}  // namespace Maptitude

#include "maptitude/Error.h"
#include "maptitude/Residue.h"
#include "maptitude/UnitCell.h"
#include "maptitude/SymOp.h"
#include "maptitude/ScatteringFactors.h"
#include "maptitude/Grid.h"
#include "maptitude/DensityScoreResult.h"
#include "maptitude/QScoreOptions.h"
#include "maptitude/RsccOptions.h"
#include "maptitude/RsrOptions.h"
#include "maptitude/CoverageOptions.h"
#include "maptitude/DensityCalculator.h"
#include "maptitude/Metric.h"
#include "maptitude/GridOps.h"
#include "maptitude/SpatialIndex.h"

#endif  // MAPTITUDE_MAPTITUDE_H
