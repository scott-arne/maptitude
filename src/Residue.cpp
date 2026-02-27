#include "maptitude/Residue.h"

#include <oechem.h>

namespace Maptitude {

Residue::Residue(std::string name, int number, std::string chain,
                 std::string insert_code)
    : name(std::move(name))
    , number(number)
    , chain(std::move(chain))
    , insert_code(std::move(insert_code)) {}

Residue Residue::FromAtom(const OEChem::OEAtomBase& atom) {
    const OEChem::OEResidue res = OEChem::OEAtomGetResidue(&atom);
    return Residue(
        res.GetName(),
        res.GetResidueNumber(),
        std::string(1, res.GetChainID()),
        std::string(1, res.GetInsertCode()));
}

std::string Residue::ToString() const {
    return name + ":" + std::to_string(number) + ":" + insert_code + ":" + chain;
}

bool Residue::operator==(const Residue& other) const {
    return name == other.name && number == other.number &&
           chain == other.chain && insert_code == other.insert_code;
}

bool Residue::operator!=(const Residue& other) const {
    return !(*this == other);
}

bool Residue::operator<(const Residue& other) const {
    if (chain != other.chain) return chain < other.chain;
    if (number != other.number) return number < other.number;
    if (insert_code != other.insert_code) return insert_code < other.insert_code;
    return name < other.name;
}

}  // namespace Maptitude

namespace std {
size_t hash<Maptitude::Residue>::operator()(const Maptitude::Residue& r) const {
    size_t h = 0;
    h ^= hash<string>()(r.name) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= hash<int>()(r.number) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= hash<string>()(r.chain) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= hash<string>()(r.insert_code) + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
}
}  // namespace std
