/**
 * @file DensityScoreResult.h
 * @brief Container for density scoring results.
 */

#ifndef MAPTITUDE_DENSITYSCORERESULT_H
#define MAPTITUDE_DENSITYSCORERESULT_H

#include <map>

#include "maptitude/Residue.h"

namespace Maptitude {

/**
 * @brief Result of a density scoring operation.
 *
 * Contains the overall score plus per-residue and per-atom breakdowns.
 *
 * @code
 * auto result = rscc(mol, grid, resolution);
 * std::cout << "Overall RSCC: " << result.overall << "\n";
 *
 * for (const auto& [res, score] : result.by_residue) {
 *     std::cout << res.ToString() << ": " << score << "\n";
 * }
 * @endcode
 */
struct DensityScoreResult {
    double overall = 0.0;                    ///< Overall score for the structure
    std::map<Residue, double> by_residue;    ///< Per-residue scores
    std::map<unsigned int, double> by_atom;  ///< Per-atom scores (keyed by atom index)

    /**
     * @brief Format as string for display.
     * @return Summary string with overall score and residue/atom counts.
     */
    [[nodiscard]] std::string ToString() const;
};

}  // namespace Maptitude

#endif  // MAPTITUDE_DENSITYSCORERESULT_H
