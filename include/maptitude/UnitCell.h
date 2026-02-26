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

}  // namespace Maptitude

#endif  // MAPTITUDE_UNITCELL_H
