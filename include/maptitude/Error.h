/**
 * @file Error.h
 * @brief Exception types for Maptitude operations.
 */

#ifndef MAPTITUDE_ERROR_H
#define MAPTITUDE_ERROR_H

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

namespace Maptitude {

/**
 * @brief Exception thrown when structure preparation or validation fails.
 *
 * @code
 * try {
 *     auto result = rscc(mol, grid, resolution);
 * } catch (const StructureError& e) {
 *     std::cerr << "Structure error: " << e.what() << "\n";
 * }
 * @endcode
 */
class StructureError : public std::runtime_error {
public:
    explicit StructureError(const std::string& message)
        : std::runtime_error(message) {}
};

/**
 * @brief Exception thrown when grid operations encounter invalid state.
 */
class GridError : public std::runtime_error {
public:
    explicit GridError(const std::string& message)
        : std::runtime_error(message) {}
};

/**
 * @brief Exception thrown when symmetry operator parsing fails.
 */
class SymOpError : public std::runtime_error {
public:
    explicit SymOpError(const std::string& message)
        : std::runtime_error(message) {}
};

/// Thrown when unit-cell parameters are geometrically invalid or describe a
/// lattice this library does not support.
class CellError : public std::runtime_error {
public:
    explicit CellError(const std::string& message) : std::runtime_error(message) {}
};

/// Reject a resolution that is not a finite positive value.
///
/// Five public entry points take a resolution, and what reaches the argument differs.
/// `ediam` and `DensityCalculator::Calculate` divide by it or by its square. `rscc`,
/// `rsr`, and `qscore` size a radius or a sweep step from it under some radius and
/// sampling models and ignore it entirely under others. Validating here rather than
/// per path keeps a resolution no caller can mean from being honoured in one
/// configuration and silently ignored in the next. `NaN <= 0.0` and `+inf <= 0.0` are
/// both false, so a sign test alone lets them through.
inline void require_usable_resolution(double resolution) {
    if (!std::isfinite(resolution) || resolution <= 0.0) {
        std::ostringstream message;
        message << "Resolution must be a finite positive value (got " << resolution << ")";
        throw GridError(message.str());
    }
}

}  // namespace Maptitude

#endif  // MAPTITUDE_ERROR_H
