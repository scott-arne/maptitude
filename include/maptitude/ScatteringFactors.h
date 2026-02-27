/**
 * @file ScatteringFactors.h
 * @brief Cromer-Mann atomic scattering factor coefficients.
 *
 * Provides compile-time tables of Cromer-Mann parameterization for
 * computing atomic form factors used in structure factor calculation.
 *
 * The scattering factor for an atom is computed as:
 *   f0(s) = c + sum_i{ a_i * exp(-b_i * s^2) }
 * where s = sin(theta)/lambda = |h|/2.
 */

#ifndef MAPTITUDE_SCATTERINGFACTORS_H
#define MAPTITUDE_SCATTERINGFACTORS_H

#include <array>
#include <cstddef>
#include <cstdint>

namespace Maptitude {

/**
 * @brief Cromer-Mann coefficients for a single atomic form factor.
 *
 * f0(s) = c + sum_i{ a[i] * exp(-b[i] * s^2) }  for i = 0..3
 */
struct CromerMannCoeffs {
    std::array<double, 4> a;  ///< Gaussian amplitude coefficients
    std::array<double, 4> b;  ///< Gaussian width coefficients
    double c;                  ///< Constant term

    /**
     * @brief Evaluate the scattering factor at a given s^2.
     *
     * @param s_squared sin^2(theta)/lambda^2 (i.e., stol^2).
     * @return Atomic form factor f0.
     */
    [[nodiscard]] double Evaluate(double s_squared) const;
};

/**
 * @brief Entry in the scattering factor lookup table.
 */
struct ScatteringFactorEntry {
    uint8_t atomic_number;    ///< Atomic number (1-98)
    int8_t formal_charge;     ///< Formal charge (-1 to +6)
    CromerMannCoeffs coeffs;  ///< Cromer-Mann coefficients
};

/**
 * @brief Look up Cromer-Mann coefficients by atomic number and charge.
 *
 * Falls back to neutral atom if the requested charge state is not available.
 *
 * @param atomic_number Atomic number (1-98).
 * @param formal_charge Formal charge (default: 0).
 * @return Pointer to coefficients, or nullptr if not found.
 */
const CromerMannCoeffs* get_scattering_factors(
    unsigned int atomic_number, int formal_charge = 0);

/**
 * @brief Get the full scattering factor table.
 *
 * @param count Output: number of entries in the table.
 * @return Pointer to the first entry.
 */
const ScatteringFactorEntry* get_scattering_factor_table(size_t& count);

}  // namespace Maptitude

#endif  // MAPTITUDE_SCATTERINGFACTORS_H
