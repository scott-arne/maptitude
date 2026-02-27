/**
 * @file Residue.h
 * @brief Residue identifier for per-residue scoring results.
 */

#ifndef MAPTITUDE_RESIDUE_H
#define MAPTITUDE_RESIDUE_H

#include <cstddef>
#include <functional>
#include <string>

namespace OEChem {
class OEAtomBase;
}

namespace Maptitude {

/**
 * @brief Identifies a unique residue by name, number, chain, and insertion code.
 *
 * Used as a key in per-residue scoring result maps. Provides comparison
 * operators for use in ordered containers and a hash specialization for
 * unordered containers.
 *
 * @code
 * auto res = Residue::FromAtom(atom);
 * std::cout << res.ToString() << std::endl;  // "ALA:123: :A"
 * @endcode
 */
struct Residue {
    std::string name;           ///< Residue name (e.g., "ALA", "GLY")
    int number = 0;             ///< PDB residue number
    std::string chain;          ///< Chain identifier (e.g., "A", "B")
    std::string insert_code;    ///< PDB insertion code (usually " ")

    /// @brief Default constructor.
    Residue() = default;

    /**
     * @brief Construct a Residue with all fields.
     *
     * @param name Residue name.
     * @param number Residue number.
     * @param chain Chain identifier.
     * @param insert_code Insertion code (default: " ").
     */
    Residue(std::string name, int number, std::string chain,
            std::string insert_code = " ");

    /**
     * @brief Create a Residue from an OpenEye atom's residue information.
     *
     * @param atom An OEAtomBase from which to extract residue data.
     * @return Residue identifying the atom's residue.
     */
    static Residue FromAtom(const OEChem::OEAtomBase& atom);

    /**
     * @brief Format as "NAME:NUMBER:ICODE:CHAIN" string.
     * @return Formatted string representation.
     */
    [[nodiscard]] std::string ToString() const;

    bool operator==(const Residue& other) const;
    bool operator!=(const Residue& other) const;
    bool operator<(const Residue& other) const;
};

}  // namespace Maptitude

// Hash specialization for use in unordered containers
namespace std {
template <>
struct hash<Maptitude::Residue> {
    size_t operator()(const Maptitude::Residue& r) const;
};
}  // namespace std

#endif  // MAPTITUDE_RESIDUE_H
