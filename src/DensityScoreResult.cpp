#include "maptitude/DensityScoreResult.h"

#include <sstream>

namespace Maptitude {

std::string DensityScoreResult::ToString() const {
    std::ostringstream oss;
    oss << std::fixed;
    oss.precision(3);
    oss << "DensityScoreResult(overall=" << overall
        << ", residues=" << by_residue.size()
        << ", atoms=" << by_atom.size() << ")";
    return oss.str();
}

}  // namespace Maptitude
