/**
 * @file UnitCell.h
 * @brief Crystallographic unit cell representation and coordinate transforms.
 */

#ifndef MAPTITUDE_UNITCELL_H
#define MAPTITUDE_UNITCELL_H

#include <array>
#include <string>

namespace Maptitude {

/**
 * @brief Crystallographic unit cell defined by lengths and angles.
 *
 * Provides coordinate transformations between fractional and Cartesian
 * coordinate systems, as well as volume and metric tensor calculations.
 *
 * @code
 * UnitCell cell(50.0, 60.0, 70.0, 90.0, 90.0, 90.0);  // Orthorhombic
 * double vol = cell.Volume();
 *
 * auto frac = cell.CartesianToFractional(10.0, 20.0, 30.0);
 * @endcode
 */
struct UnitCell {
    double a = 0.0;      ///< Cell length a (Angstroms)
    double b = 0.0;      ///< Cell length b (Angstroms)
    double c = 0.0;      ///< Cell length c (Angstroms)
    double alpha = 90.0;  ///< Angle between b and c (degrees)
    double beta = 90.0;   ///< Angle between a and c (degrees)
    double gamma = 90.0;  ///< Angle between a and b (degrees)

    /// @brief Default constructor.
    UnitCell() = default;

    /**
     * @brief Construct unit cell from dimensions.
     *
     * @param a Cell length a in Angstroms.
     * @param b Cell length b in Angstroms.
     * @param c Cell length c in Angstroms.
     * @param alpha Angle alpha in degrees.
     * @param beta Angle beta in degrees.
     * @param gamma Angle gamma in degrees.
     */
    UnitCell(double a, double b, double c,
             double alpha, double beta, double gamma);

    /**
     * @brief Compute the unit cell volume.
     * @return Volume in cubic Angstroms.
     */
    [[nodiscard]] double Volume() const;

    /**
     * @brief Get the orthogonalization matrix (fractional -> Cartesian).
     *
     * The 3x3 matrix M such that (x,y,z)_cart = M * (u,v,w)_frac.
     * Stored in row-major order.
     *
     * @return 9-element array representing the 3x3 matrix.
     */
    [[nodiscard]] std::array<double, 9> OrthogonalizationMatrix() const;

    /**
     * @brief Get the deorthogonalization matrix (Cartesian -> fractional).
     *
     * The 3x3 matrix M^-1 such that (u,v,w)_frac = M^-1 * (x,y,z)_cart.
     * Stored in row-major order.
     *
     * @return 9-element array representing the 3x3 matrix.
     */
    [[nodiscard]] std::array<double, 9> DeorthogonalizationMatrix() const;

    /**
     * @brief Convert Cartesian coordinates to fractional.
     *
     * @param x Cartesian x coordinate.
     * @param y Cartesian y coordinate.
     * @param z Cartesian z coordinate.
     * @return Fractional coordinates {u, v, w}.
     */
    [[nodiscard]] std::array<double, 3> CartesianToFractional(
        double x, double y, double z) const;

    /**
     * @brief Convert fractional coordinates to Cartesian.
     *
     * @param u Fractional coordinate along a.
     * @param v Fractional coordinate along b.
     * @param w Fractional coordinate along c.
     * @return Cartesian coordinates {x, y, z} in Angstroms.
     */
    [[nodiscard]] std::array<double, 3> FractionalToCartesian(
        double u, double v, double w) const;

    /**
     * @brief Format as string for display.
     * @return String like "UnitCell(a=50.0, b=60.0, c=70.0, alpha=90.0, beta=90.0, gamma=90.0)".
     */
    [[nodiscard]] std::string ToString() const;

    bool operator==(const UnitCell& other) const;
    bool operator!=(const UnitCell& other) const;
};

/// Validate that a unit cell is geometrically valid and numerically usable.
///
/// Runs a sequence of checks: (1) each length is finite and positive; (2) each
/// angle is strictly within (0, 180) degrees; (3) the volume radicand exceeds a
/// small positive floor (~1e-9), rejecting both geometrically impossible cells and
/// those too degenerate to compute with; (4) the derived volume is finite and
/// positive; (5) all orthogonalization and deorthogonalization matrix entries are
/// finite and the diagonal is nonzero; (6) the infinity-norm condition number
/// (product of the matrix row-sum norms) stays below a threshold (~1e9), ensuring
/// coordinate conversions are numerically stable for any interior point. Neither
/// the angle-range check nor the radicand check subsumes the other: cos(200°) ==
/// cos(160°), so an out-of-range angle can pass the radicand test, while in-range
/// angles (150, 150, 150) can describe no real lattice.
///
/// Called automatically by the parameterized constructor, by the geometry readers
/// (Volume, OrthogonalizationMatrix, DeorthogonalizationMatrix), and by
/// DensityCalculator's constructor. The default constructor does not validate.
///
/// \param cell The cell to check.
/// \throws CellError If any check fails; the message names the specific cause.
void validate_cell(const UnitCell& cell);

}  // namespace Maptitude

#endif  // MAPTITUDE_UNITCELL_H
