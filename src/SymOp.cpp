#include "maptitude/SymOp.h"
#include "maptitude/Error.h"

#include <algorithm>
#include <cmath>
#include <regex>
#include <sstream>

namespace Maptitude {

namespace {

// Parse a single component of a symmetry operator (e.g., "-x", "y+1/2", "z+1/4")
void ParseComponent(const std::string& component, int row, std::array<double, 9>& R,
                    std::array<double, 3>& t) {
    std::string s = component;
    // Remove whitespace
    s.erase(std::remove(s.begin(), s.end(), ' '), s.end());

    double coeff = 0.0;
    double trans = 0.0;
    size_t pos = 0;
    double sign = 1.0;

    while (pos < s.size()) {
        char ch = s[pos];

        if (ch == '+') {
            sign = 1.0;
            ++pos;
        } else if (ch == '-') {
            sign = -1.0;
            ++pos;
        } else if (ch == 'x' || ch == 'X') {
            R[row * 3 + 0] = sign;
            sign = 1.0;
            ++pos;
        } else if (ch == 'y' || ch == 'Y') {
            R[row * 3 + 1] = sign;
            sign = 1.0;
            ++pos;
        } else if (ch == 'z' || ch == 'Z') {
            R[row * 3 + 2] = sign;
            sign = 1.0;
            ++pos;
        } else if (std::isdigit(ch) || ch == '.') {
            // Parse a number - could be a fraction numerator or decimal
            size_t end;
            double num = std::stod(s.substr(pos), &end);
            pos += end;
            if (pos < s.size() && s[pos] == '/') {
                ++pos;
                double denom = std::stod(s.substr(pos), &end);
                pos += end;
                trans += sign * num / denom;
            } else {
                trans += sign * num;
            }
            sign = 1.0;
        } else {
            throw SymOpError("Unexpected character '" + std::string(1, ch) +
                             "' in symmetry operator component: " + component);
        }
    }

    t[row] = trans;
}

}  // namespace

SymOp::SymOp(std::array<double, 9> rotation, std::array<double, 3> translation)
    : R(std::move(rotation)), t(std::move(translation)) {}

SymOp SymOp::Parse(const std::string& triplet) {
    // Split by comma
    std::vector<std::string> parts;
    std::istringstream iss(triplet);
    std::string part;
    while (std::getline(iss, part, ',')) {
        parts.push_back(part);
    }

    if (parts.size() != 3) {
        throw SymOpError("Symmetry operator must have exactly 3 components: " + triplet);
    }

    SymOp op;
    op.R.fill(0.0);
    op.t.fill(0.0);

    for (int i = 0; i < 3; ++i) {
        ParseComponent(parts[i], i, op.R, op.t);
    }

    return op;
}

std::vector<SymOp> SymOp::ParseAll(const std::string& text) {
    std::vector<SymOp> result;
    std::istringstream iss(text);
    std::string line;

    while (std::getline(iss, line)) {
        // Trim whitespace
        size_t start = line.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) continue;
        size_t end = line.find_last_not_of(" \t\r\n");
        line = line.substr(start, end - start + 1);

        if (line.empty()) continue;

        // Handle semicolon-separated operators on a single line
        std::istringstream line_stream(line);
        std::string op_str;
        while (std::getline(line_stream, op_str, ';')) {
            size_t s = op_str.find_first_not_of(" \t");
            if (s == std::string::npos) continue;
            size_t e = op_str.find_last_not_of(" \t");
            op_str = op_str.substr(s, e - s + 1);
            if (!op_str.empty()) {
                result.push_back(Parse(op_str));
            }
        }
    }

    return result;
}

std::array<double, 3> SymOp::Apply(double u, double v, double w) const {
    return {
        R[0] * u + R[1] * v + R[2] * w + t[0],
        R[3] * u + R[4] * v + R[5] * w + t[1],
        R[6] * u + R[7] * v + R[8] * w + t[2]
    };
}

std::string SymOp::ToString() const {
    // Reconstruct triplet notation from matrix + translation
    std::ostringstream oss;
    const char* axes[] = {"x", "y", "z"};

    for (int row = 0; row < 3; ++row) {
        if (row > 0) oss << ",";
        bool first = true;

        for (int col = 0; col < 3; ++col) {
            double val = R[row * 3 + col];
            if (std::abs(val) < 1e-10) continue;

            if (val > 0 && !first) oss << "+";
            if (std::abs(val - 1.0) < 1e-10) {
                oss << axes[col];
            } else if (std::abs(val + 1.0) < 1e-10) {
                oss << "-" << axes[col];
            } else {
                oss << val << "*" << axes[col];
            }
            first = false;
        }

        if (std::abs(t[row]) > 1e-10) {
            if (t[row] > 0 && !first) oss << "+";
            // Try to express as fraction
            double frac = t[row];
            bool found_frac = false;
            for (int denom = 2; denom <= 12; ++denom) {
                double numer = frac * denom;
                if (std::abs(numer - std::round(numer)) < 1e-8) {
                    int n = static_cast<int>(std::round(numer));
                    oss << n << "/" << denom;
                    found_frac = true;
                    break;
                }
            }
            if (!found_frac) {
                oss << frac;
            }
        }
    }

    return oss.str();
}

bool SymOp::operator==(const SymOp& other) const {
    return R == other.R && t == other.t;
}

bool SymOp::operator!=(const SymOp& other) const {
    return !(*this == other);
}

}  // namespace Maptitude
