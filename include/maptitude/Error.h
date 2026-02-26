/**
 * @file Error.h
 * @brief Exception types for Maptitude operations.
 */

#ifndef MAPTITUDE_ERROR_H
#define MAPTITUDE_ERROR_H

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

}  // namespace Maptitude

#endif  // MAPTITUDE_ERROR_H
