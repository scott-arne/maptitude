#include "maptitude/UnitCell.h"

#include <cmath>
#include <sstream>

#include "maptitude/Error.h"

namespace Maptitude {

namespace {
constexpr double PI = 3.14159265358979323846;
constexpr double DEG_TO_RAD = PI / 180.0;
}

UnitCell::UnitCell(const double a, const double b, const double c,
                   const double alpha, const double beta, const double gamma)
    : a(a), b(b), c(c), alpha(alpha), beta(beta), gamma(gamma) {
    validate_cell(*this);
}

double UnitCell::Volume() const {
    const double ca = std::cos(alpha * DEG_TO_RAD);
    const double cb = std::cos(beta * DEG_TO_RAD);
    const double cg = std::cos(gamma * DEG_TO_RAD);
    return a * b * c * std::sqrt(1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg);
}

std::array<double, 9> UnitCell::OrthogonalizationMatrix() const {
    const double ca = std::cos(alpha * DEG_TO_RAD);
    const double cb = std::cos(beta * DEG_TO_RAD);
    const double cg = std::cos(gamma * DEG_TO_RAD);
    const double sg = std::sin(gamma * DEG_TO_RAD);

    const double vol = Volume();

    // Row-major 3x3 matrix: fractional -> Cartesian
    return {
        a,  b * cg,  c * cb,
        0,  b * sg,  c * (ca - cb * cg) / sg,
        0,  0,       vol / (a * b * sg)
    };
}

std::array<double, 9> UnitCell::DeorthogonalizationMatrix() const {
    const double ca = std::cos(alpha * DEG_TO_RAD);
    const double cb = std::cos(beta * DEG_TO_RAD);
    const double cg = std::cos(gamma * DEG_TO_RAD);
    const double sg = std::sin(gamma * DEG_TO_RAD);

    const double vol = Volume();

    // Inverse of orthogonalization matrix: Cartesian -> fractional
    return {
        1.0 / a,  -cg / (a * sg),  (cg * (ca - cb * cg) / sg - cb * sg) * b * c / vol,
        0,        1.0 / (b * sg),  -(ca - cb * cg) * a * c / (vol * sg),
        0,        0,               a * b * sg / vol
    };
}

std::array<double, 3> UnitCell::CartesianToFractional(
    const double x, const double y, const double z) const {
    const auto M = DeorthogonalizationMatrix();
    return {
        M[0] * x + M[1] * y + M[2] * z,
        M[3] * x + M[4] * y + M[5] * z,
        M[6] * x + M[7] * y + M[8] * z
    };
}

std::array<double, 3> UnitCell::FractionalToCartesian(
    const double u, const double v, const double w) const {
    const auto M = OrthogonalizationMatrix();
    return {
        M[0] * u + M[1] * v + M[2] * w,
        M[3] * u + M[4] * v + M[5] * w,
        M[6] * u + M[7] * v + M[8] * w
    };
}

std::string UnitCell::ToString() const {
    std::ostringstream oss;
    oss << "UnitCell(a=" << a << ", b=" << b << ", c=" << c
        << ", alpha=" << alpha << ", beta=" << beta << ", gamma=" << gamma << ")";
    return oss.str();
}

bool UnitCell::operator==(const UnitCell& other) const {
    return a == other.a && b == other.b && c == other.c &&
           alpha == other.alpha && beta == other.beta && gamma == other.gamma;
}

bool UnitCell::operator!=(const UnitCell& other) const {
    return !(*this == other);
}

void validate_cell(const UnitCell& cell) {
    const double lengths[3] = {cell.a, cell.b, cell.c};
    const char* length_names[3] = {"a", "b", "c"};
    for (int i = 0; i < 3; ++i) {
        if (!std::isfinite(lengths[i]) || lengths[i] <= 0.0) {
            std::ostringstream message;
            message << "Unit cell length " << length_names[i] << " must be finite and positive (got "
                    << lengths[i] << ")";
            throw CellError(message.str());
        }
    }

    const double angles[3] = {cell.alpha, cell.beta, cell.gamma};
    const char* angle_names[3] = {"alpha", "beta", "gamma"};
    for (int i = 0; i < 3; ++i) {
        if (!std::isfinite(angles[i]) || angles[i] <= 0.0 || angles[i] >= 180.0) {
            std::ostringstream message;
            message << "Unit cell angle " << angle_names[i]
                    << " must be strictly between 0 and 180 degrees (got " << angles[i] << ")";
            throw CellError(message.str());
        }
    }

    const double ca = std::cos(cell.alpha * DEG_TO_RAD);
    const double cb = std::cos(cell.beta * DEG_TO_RAD);
    const double cg = std::cos(cell.gamma * DEG_TO_RAD);
    const double radicand = 1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg;
    if (radicand <= 0.0) {
        std::ostringstream message;
        message << "Unit cell angles (" << cell.alpha << ", " << cell.beta << ", " << cell.gamma
                << ") describe no real lattice: volume radicand is " << radicand;
        throw CellError(message.str());
    }
}

}  // namespace Maptitude
