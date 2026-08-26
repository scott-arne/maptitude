#include "maptitude/UnitCell.h"

#include <cmath>
#include <sstream>

#include "maptitude/Error.h"

namespace Maptitude {

namespace {
constexpr double PI = 3.14159265358979323846;
constexpr double DEG_TO_RAD = PI / 180.0;

// Minimum radicand for a cell to be considered numerically usable. The radicand
// is a dimensionless function of the angles alone (bounded above by 1 for
// orthogonal cells). The direct polynomial form (1 - ca² - cb² - cg² + 2·ca·cb·cg)
// has absolute error ~1e-15 from rounding (each term is O(1), so it is a few ulps
// of 1.0 regardless of how small the true radicand is), making the computed sign
// unreliable for near-degenerate cells. This floor sits ~6 orders above the noise
// (so the sign question never arises) and ~7 orders below realistic crystallography
// (the most degenerate angle triple in the test suite, 30/30/40 degrees, has
// radicand 0.062).
constexpr double MIN_RADICAND = 1.0e-9;

// Maximum infinity-norm condition number for coordinate conversion. Standard
// rounding analysis for fl(M^-1 fl(M p)) gives round-trip error <= 6u*kappa*||p||_inf,
// with u = 1.11e-16 and ||p||_inf <= 1 for fractional coordinates, so kappa < 1e9
// bounds error at 6.7e-7, comfortably under the 1e-6 target. Measured across
// plausible cells (0.5-5000 A, 20-160 deg), the worst kappa is 15,230 — 4.8 orders
// below this threshold. Named cells: 1 A cube and 3000 A cube both kappa 1; seven
// canonical shapes span 1 to 10.8; 500x500x1200 ribosome 3.79; 2x2x1000 thin plate 500.
constexpr double MAX_CONDITION_NUMBER = 1.0e9;

std::array<double, 9> BuildOrthogonalizationMatrix(double a, double b, double c,
                                                   double ca, double cb, double cg,
                                                   double sg, double vol) {
    // Row-major 3x3 matrix: fractional -> Cartesian
    return {
        a,  b * cg,  c * cb,
        0,  b * sg,  c * (ca - cb * cg) / sg,
        0,  0,       vol / (a * b * sg)
    };
}

std::array<double, 9> BuildDeorthogonalizationMatrix(double a, double b, double c,
                                                     double ca, double cb, double cg,
                                                     double sg, double vol) {
    // Inverse of orthogonalization matrix: Cartesian -> fractional
    return {
        1.0 / a,  -cg / (a * sg),  (cg * (ca - cb * cg) / sg - cb * sg) * b * c / vol,
        0,        1.0 / (b * sg),  -(ca - cb * cg) * a * c / (vol * sg),
        0,        0,               a * b * sg / vol
    };
}
}

UnitCell::UnitCell(const double a, const double b, const double c,
                   const double alpha, const double beta, const double gamma)
    : a(a), b(b), c(c), alpha(alpha), beta(beta), gamma(gamma) {
    validate_cell(*this);
}

double UnitCell::Volume() const {
    validate_cell(*this);
    const double ca = std::cos(alpha * DEG_TO_RAD);
    const double cb = std::cos(beta * DEG_TO_RAD);
    const double cg = std::cos(gamma * DEG_TO_RAD);
    return a * b * c * std::sqrt(1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg);
}

std::array<double, 9> UnitCell::OrthogonalizationMatrix() const {
    validate_cell(*this);
    const double ca = std::cos(alpha * DEG_TO_RAD);
    const double cb = std::cos(beta * DEG_TO_RAD);
    const double cg = std::cos(gamma * DEG_TO_RAD);
    const double sg = std::sin(gamma * DEG_TO_RAD);

    const double vol = Volume();

    return BuildOrthogonalizationMatrix(a, b, c, ca, cb, cg, sg, vol);
}

std::array<double, 9> UnitCell::DeorthogonalizationMatrix() const {
    validate_cell(*this);
    const double ca = std::cos(alpha * DEG_TO_RAD);
    const double cb = std::cos(beta * DEG_TO_RAD);
    const double cg = std::cos(gamma * DEG_TO_RAD);
    const double sg = std::sin(gamma * DEG_TO_RAD);

    const double vol = Volume();

    return BuildDeorthogonalizationMatrix(a, b, c, ca, cb, cg, sg, vol);
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
    if (radicand < MIN_RADICAND) {
        std::ostringstream message;
        message << "Unit cell angles (" << cell.alpha << ", " << cell.beta << ", " << cell.gamma
                << ") are either geometrically impossible or too degenerate: volume radicand is "
                << radicand;
        throw CellError(message.str());
    }

    // Compute the volume locally rather than calling Volume(), which validates
    // its own cell and would recurse back into here.
    //
    // A finite positive volume is not on its own enough to make a cell usable:
    // vol = a*b*c*sqrt(radicand) is left-associative, so a large c can mask an
    // a*b that has already underflowed to zero, leaving denominators such as
    // a*b*sg at zero while vol still looks healthy. Check the derived matrices
    // themselves instead of reasoning about which denominators can degenerate.
    const double sg = std::sin(cell.gamma * DEG_TO_RAD);
    const double vol = cell.a * cell.b * cell.c * std::sqrt(radicand);
    if (!std::isfinite(vol) || vol <= 0.0) {
        std::ostringstream message;
        message << "Unit cell produces unusable volume: " << vol;
        throw CellError(message.str());
    }

    const auto ortho = BuildOrthogonalizationMatrix(cell.a, cell.b, cell.c, ca, cb, cg, sg, vol);
    const auto deortho = BuildDeorthogonalizationMatrix(cell.a, cell.b, cell.c, ca, cb, cg, sg, vol);

    // Verify all matrix entries are finite and the diagonal is nonzero (a zero
    // diagonal entry means a degenerate transform). Indices 0, 4, 8 are the diagonal.
    for (int i = 0; i < 9; ++i) {
        if (!std::isfinite(ortho[i])) {
            throw CellError("Unit cell produces non-finite orthogonalization matrix");
        }
        if (!std::isfinite(deortho[i])) {
            throw CellError("Unit cell produces non-finite deorthogonalization matrix");
        }
    }
    if (ortho[0] == 0.0 || ortho[4] == 0.0 || ortho[8] == 0.0) {
        throw CellError("Unit cell produces degenerate orthogonalization matrix (zero diagonal)");
    }
    if (deortho[0] == 0.0 || deortho[4] == 0.0 || deortho[8] == 0.0) {
        throw CellError("Unit cell produces degenerate deorthogonalization matrix (zero diagonal)");
    }

    // Infinity-norm condition number: bounds the round-trip error for any interior
    // point. Finite matrix entries and nonzero diagonal are not sufficient: cells
    // with extreme length ratios can have well-formed matrices yet produce
    // catastrophically wrong coordinates. A condition-number bound is a *bound*
    // rather than a *sample*, so it is not blind to any failure mode the way
    // short-mantissa probes are. Compute from local arrays; do not call the
    // readers, which call validate_cell, creating infinite recursion.
    const double ortho_row0 = std::abs(ortho[0]) + std::abs(ortho[1]) + std::abs(ortho[2]);
    const double ortho_row1 = std::abs(ortho[3]) + std::abs(ortho[4]) + std::abs(ortho[5]);
    const double ortho_row2 = std::abs(ortho[6]) + std::abs(ortho[7]) + std::abs(ortho[8]);
    const double ortho_norm = std::max({ortho_row0, ortho_row1, ortho_row2});

    const double deortho_row0 = std::abs(deortho[0]) + std::abs(deortho[1]) + std::abs(deortho[2]);
    const double deortho_row1 = std::abs(deortho[3]) + std::abs(deortho[4]) + std::abs(deortho[5]);
    const double deortho_row2 = std::abs(deortho[6]) + std::abs(deortho[7]) + std::abs(deortho[8]);
    const double deortho_norm = std::max({deortho_row0, deortho_row1, deortho_row2});

    const double kappa = ortho_norm * deortho_norm;

    // Require kappa finite before checking threshold: some cells overflow a row
    // sum to inf while every individual matrix entry is finite.
    if (!std::isfinite(kappa) || kappa >= MAX_CONDITION_NUMBER) {
        std::ostringstream message;
        message << "Unit cell is numerically unusable for coordinate conversion (condition number "
                << kappa << " exceeds threshold " << MAX_CONDITION_NUMBER << ")";
        throw CellError(message.str());
    }
}

}  // namespace Maptitude
