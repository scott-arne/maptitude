/**
 * @file SymOp.h
 * @brief Crystallographic symmetry operator representation and parsing.
 */

#ifndef MAPTITUDE_SYMOP_H
#define MAPTITUDE_SYMOP_H

#include <array>
#include <string>
#include <vector>

namespace Maptitude {

/**
 * @brief A crystallographic symmetry operator consisting of a 3x3 rotation
 *        matrix and a 3-element translation vector.
 *
 * Symmetry operators are expressed in fractional coordinates. The operator
 * transforms a point (x,y,z) as: (x',y',z') = R*(x,y,z) + t
 *
 * @code
 * auto op = SymOp::Parse("x,y,z");        // Identity
 * auto op2 = SymOp::Parse("-x,y+1/2,-z"); // Screw axis
 *
 * auto ops = SymOp::ParseAll("x,y,z\n-x,y+1/2,-z");
 * @endcode
 */
struct SymOp {
    std::array<double, 9> R = {1, 0, 0, 0, 1, 0, 0, 0, 1};  ///< Rotation matrix (row-major)
    std::array<double, 3> t = {0, 0, 0};                       ///< Translation vector

    /// @brief Default constructor creates the identity operator.
    SymOp() = default;

    /**
     * @brief Construct from rotation matrix and translation vector.
     *
     * @param rotation 9-element rotation matrix in row-major order.
     * @param translation 3-element translation vector.
     */
    SymOp(std::array<double, 9> rotation, std::array<double, 3> translation);

    /**
     * @brief Parse a single symmetry operator from a triplet string.
     *
     * Supports standard crystallographic notation such as:
     * - "x,y,z" (identity)
     * - "-x,y+1/2,-z" (fractional translations)
     * - "y,-x,z+1/4"
     *
     * @param triplet Symmetry operator in triplet notation.
     * @return Parsed SymOp.
     * @throws SymOpError if parsing fails.
     */
    static SymOp Parse(const std::string& triplet);

    /**
     * @brief Parse multiple symmetry operators from newline/semicolon-separated text.
     *
     * @param text Block of symmetry operators, one per line or semicolon-separated.
     * @return Vector of parsed SymOp objects.
     * @throws SymOpError if any operator fails to parse.
     */
    static std::vector<SymOp> ParseAll(const std::string& text);

    /**
     * @brief Apply this operator to fractional coordinates.
     *
     * @param u Fractional coordinate along a.
     * @param v Fractional coordinate along b.
     * @param w Fractional coordinate along c.
     * @return Transformed fractional coordinates {u', v', w'}.
     */
    [[nodiscard]] std::array<double, 3> Apply(double u, double v, double w) const;

    /**
     * @brief Format as string for display.
     * @return String representation in triplet notation.
     */
    [[nodiscard]] std::string ToString() const;

    bool operator==(const SymOp& other) const;
    bool operator!=(const SymOp& other) const;
};

}  // namespace Maptitude

#endif  // MAPTITUDE_SYMOP_H
