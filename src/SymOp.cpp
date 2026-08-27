#include "maptitude/SymOp.h"
#include "maptitude/Error.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>

namespace Maptitude {

namespace {

// std::stod signals failure with std::invalid_argument and std::out_of_range, and
// neither is a SymOpError. A symmetry operator is user data out of a CCP4 header or
// a CIF file, so a malformed one has to reach the caller as the library's own typed
// exception rather than through the generic std::exception arm of the SWIG wrapper.
double ParseNumber(const std::string& text, const std::string& component, size_t& consumed) {
    // std::stod accepts more than this grammar intends: a leading sign, the named
    // literals "inf"/"infinity"/"nan", and C99 hex floats such as "0x10". None is a
    // crystallographic translation, and the arithmetic downstream launders some of
    // them into ordinary values -- 1/inf is a finite 0.0, so neither the zero-
    // denominator guard nor the finiteness check on the accumulated translation ever
    // sees it, and "x+1/inf,y,z" parses as the identity. Constrain the token here,
    // where it is still a string, instead of trying to detect it after the fact.
    if (text.empty() || !(std::isdigit(static_cast<unsigned char>(text[0])) || text[0] == '.')) {
        throw SymOpError("Expected a number in symmetry operator component: " + component);
    }

    double value;
    try {
        value = std::stod(text, &consumed);
    } catch (const std::invalid_argument&) {
        throw SymOpError("Expected a number in symmetry operator component: " + component);
    } catch (const std::out_of_range&) {
        throw SymOpError("Number out of range in symmetry operator component: " + component);
    }

    // The leading-character test admits "0x10", whose first character is a digit.
    // Every character std::stod actually consumed has to belong to an unsigned
    // decimal with an optional exponent.
    const std::string token = text.substr(0, consumed);
    for (const char ch : token) {
        if (!std::isdigit(static_cast<unsigned char>(ch)) && ch != '.' && ch != 'e' &&
            ch != 'E' && ch != '+' && ch != '-') {
            throw SymOpError("Malformed number in symmetry operator component: " + component);
        }
    }

    return value;
}

// Parse a single component of a symmetry operator (e.g., "-x", "y+1/2", "z+1/4")
void ParseComponent(const std::string& component, const int row, std::array<double, 9>& R,
                    std::array<double, 3>& t) {
    std::string s = component;
    // Remove whitespace
    s.erase(std::remove(s.begin(), s.end(), ' '), s.end());

    double trans = 0.0;
    size_t pos = 0;
    double sign = 1.0;
    bool axis_seen[3] = {false, false, false};
    // 'sign' alone cannot tell a leading '-' from a stray one: it is just a
    // multiplier the next term happens to pick up. These two flags record the
    // grammar state the multiplier does not -- whether a sign is still owed a term,
    // and whether the component contained any term at all.
    bool pending_sign = false;
    bool saw_term = false;

    while (pos < s.size()) {
        const char ch = s[pos];

        if (ch == '+' || ch == '-') {
            if (pending_sign) {
                throw SymOpError("Consecutive signs in symmetry operator component: " +
                                 component);
            }
            sign = (ch == '-') ? -1.0 : 1.0;
            pending_sign = true;
            ++pos;
        } else if (ch == 'x' || ch == 'X' || ch == 'y' || ch == 'Y' || ch == 'z' || ch == 'Z') {
            // Terms are separated by an explicit '+' or '-'. Without this the loop reads
            // juxtaposition as addition: "xy" becomes x+y, gaining a rotation column, and
            // "2x" becomes x+2, gaining a translation. Neither is a symmetry operator, and
            // both are silently rewritten into one that looks valid.
            if (saw_term && !pending_sign) {
                throw SymOpError(
                    "Missing '+' or '-' between terms in symmetry operator component: " +
                    component);
            }
            // Each axis may appear at most once per component. A repeated axis is not a
            // symmetry operation, and the write below would silently overwrite the earlier
            // coefficient rather than sum it -- "x+x" parsing as 1 and "x-x" as -1. Note
            // this is a per-axis rule, not one non-zero entry per row: trigonal and
            // hexagonal operators such as "x-y,x,z" legitimately populate two columns.
            const int axis = (ch == 'x' || ch == 'X') ? 0 : (ch == 'y' || ch == 'Y') ? 1 : 2;
            if (axis_seen[axis]) {
                throw SymOpError("Axis '" + std::string(1, ch) +
                                 "' appears more than once in symmetry operator component: " +
                                 component);
            }
            axis_seen[axis] = true;
            R[row * 3 + axis] = sign;
            sign = 1.0;
            ++pos;
            pending_sign = false;
            saw_term = true;
        // The <cctype> classifiers are only defined for values representable as
        // unsigned char; a negative char is undefined behavior.
        } else if (std::isdigit(static_cast<unsigned char>(ch)) || ch == '.') {
            // Same term-boundary rule as the axis branch: "x.5" and "x1/2" would otherwise
            // both become x+1/2.
            if (saw_term && !pending_sign) {
                throw SymOpError(
                    "Missing '+' or '-' between terms in symmetry operator component: " +
                    component);
            }
            size_t end;
            const double num = ParseNumber(s.substr(pos), component, end);
            pos += end;
            if (pos < s.size() && s[pos] == '/') {
                ++pos;
                const double denom = ParseNumber(s.substr(pos), component, end);
                pos += end;
                if (denom == 0.0) {
                    throw SymOpError("Zero denominator in symmetry operator component: " +
                                     component);
                }
                trans += sign * num / denom;
            } else {
                trans += sign * num;
            }
            sign = 1.0;
            pending_sign = false;
            saw_term = true;
        } else {
            throw SymOpError("Unexpected character '" + std::string(1, ch) +
                             "' in symmetry operator component: " + component);
        }
    }

    if (pending_sign) {
        throw SymOpError("Symmetry operator component ends with a sign: " + component);
    }
    if (!saw_term) {
        throw SymOpError("Empty symmetry operator component: " + component);
    }
    // A zero denominator is caught above with a specific message, but summed
    // translations can overflow on their own ("1e308+1e308"), so check the result.
    if (!std::isfinite(trans)) {
        throw SymOpError("Non-finite translation in symmetry operator component: " + component);
    }

    t[row] = trans;
}

}  // namespace

SymOp::SymOp(std::array<double, 9> rotation, std::array<double, 3> translation)
    : R(rotation), t(translation) {}

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
        const size_t start = line.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) continue;
        const size_t end = line.find_last_not_of(" \t\r\n");
        line = line.substr(start, end - start + 1);

        if (line.empty()) continue;

        // Handle semicolon-separated operators on a single line
        std::istringstream line_stream(line);
        std::string op_str;
        while (std::getline(line_stream, op_str, ';')) {
            const size_t s = op_str.find_first_not_of(" \t");
            if (s == std::string::npos) continue;
            const size_t e = op_str.find_last_not_of(" \t");
            op_str = op_str.substr(s, e - s + 1);
            if (!op_str.empty()) {
                result.push_back(Parse(op_str));
            }
        }
    }

    return result;
}

std::array<double, 3> SymOp::Apply(const double u, const double v, const double w) const {
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
            const double val = R[row * 3 + col];
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
            const double frac = t[row];
            bool found_frac = false;
            for (int denom = 2; denom <= 12; ++denom) {
                const double numer = frac * denom;
                if (std::abs(numer - std::round(numer)) < 1e-8) {
                    const double rounded = std::round(numer);
                    // Casting a double outside int's range is undefined behavior, and
                    // frac * 2 passes INT_MAX once the translation reaches 2^30. The
                    // saturated cast silently rewrote the operator: "x+1073741824,y,z"
                    // serialized as "x+2147483647/2,y,z" and reparsed as 1073741823.5.
                    // A larger denominator only grows |numer|, so give up on the fraction
                    // form entirely and let the max_digits10 decimal path below render it
                    // exactly.
                    if (rounded < static_cast<double>(std::numeric_limits<int>::min()) ||
                        rounded > static_cast<double>(std::numeric_limits<int>::max())) {
                        break;
                    }
                    oss << static_cast<int>(rounded) << "/" << denom;
                    found_frac = true;
                    break;
                }
            }
            if (!found_frac) {
                // Six significant digits is not enough to name a double: 1/13 writes
                // as 0.0769231 and reads back as a different value. max_digits10 is
                // the shortest precision that round-trips every double exactly.
                const std::streamsize previous =
                    oss.precision(std::numeric_limits<double>::max_digits10);
                oss << frac;
                oss.precision(previous);
            }
            first = false;
        }

        // ParseComponent rejects an empty component, so a row that emitted nothing -- an
        // all-zero rotation and a translation below the suppression threshold -- would make
        // the serializer produce a string its own parser refuses. "0,y,z" round trips;
        // ",y,z" does not parse at all.
        if (first) {
            oss << "0";
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
