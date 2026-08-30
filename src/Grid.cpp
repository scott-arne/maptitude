#include "maptitude/Grid.h"
#include "maptitude/Error.h"

#include <oegrid.h>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>

namespace Maptitude {

namespace {

/// Bound on one double -> float -> double round trip, relative to the magnitude
/// of the value rounded: half an ulp of float.
///
/// Nothing about a grid's geometry is stored in double. Every node coordinate
/// reaches maptitude through ElementToSpatialCoord as a float, so this, scaled by
/// the magnitudes an axis's geometry actually takes, is the noise floor on both
/// the node span and the sampled extent. The scaling is the point: a grid at the
/// Cartesian origin shows only the fractional-to-Cartesian matrix's
/// cos(90 deg) ~ 6.1e-17 term, around 1.5e-15 A, while the same five-node grid
/// moved to 12.3 A reports its first node 1.9e-7 A away -- eight orders larger.
/// A fixed tolerance sized for the first case rejects a caller querying the
/// second grid's own corner.
constexpr double FLOAT_HALF_ULP = 5.9604644775390625e-08;  // 0x1p-24

/// Float roundings that separate a caller's nominal node coordinate from the
/// derived end of the node span.
///
/// Six, each bounded by FLOAT_HALF_ULP times AxisMagnitude: the grid centre and
/// the unit-cell edge as OpenEye stores them, the division of that edge by the
/// dim that gives its internal node interval, the float the endpoint coordinate
/// is returned as, and the two endpoint coordinates a second time, re-entering
/// through the interval get_grid_params derives by walking the full span -- an
/// interval whose relative error is multiplied back by |f|, which at the far
/// face is the n - 1 it was divided by.
///
/// A sweep of every axis this file's tests can build -- dims 2 to 256, 97
/// spacings from 0.15 to 4.1 A, origins from -span to +1000 A, some 270,000
/// grids -- puts the worst node-span error at 2.63 of these units, so the count
/// is a bound with better than a factor of two in hand.
constexpr double NODE_SPAN_ROUNDINGS = 6.0;

/// Float roundings that separate a caller's cell edge from the extent the grid
/// samples, n * spacing.
///
/// Seven are countable, each bounded by FLOAT_HALF_ULP times AxisMagnitude:
/// three at the edge's own magnitude (the edge as OpenEye stores it, the
/// division by the dim, and the caller's edge, itself a float whenever it came
/// from a CCP4 header or GetUnitCell) and four from the two endpoint
/// coordinates, whose roundings are amplified by the n / (n - 1) rescale from
/// the (n - 1)-interval walk to an n-interval extent -- at most 2 each, at the
/// two-node grid this still accepts. Eight is used rather than seven because
/// the count models OpenEye's fractional-to-Cartesian step as a single division
/// where it is in fact a stored skew matrix applied in float.
///
/// The same 270,000-grid sweep puts the worst extent error at 3.95 of these
/// units. Note this cannot borrow same_grid_geometry's tolerance, which is
/// generous for a different reason: that comparison has a derived quantity on
/// both sides, so the float noise is common-mode and mostly cancels. Here one
/// side is the caller's own number and nothing cancels.
constexpr double CELL_EXTENT_ROUNDINGS = 8.0;

/// Read one node's Cartesian coordinate, enforcing derivation checks 2 and 3.
/// @p axis names the axis whose walk needs this node, for the error text.
void ReadNodeCoord(const OESystem::OESkewGrid& grid, const unsigned int element,
                   const char* axis, double out[3]) {
    static const char* const COMPONENT[3] = {"x", "y", "z"};
    float f[3] = {0.0f, 0.0f, 0.0f};

    if (!grid.ElementToSpatialCoord(element, f[0], f[1], f[2])) {
        std::ostringstream message;
        message << "Grid element " << element << " (walking axis " << axis
                << ") has no spatial coordinate; the geometry cannot be derived";
        throw GridError(message.str());
    }

    for (int j = 0; j < 3; ++j) {
        if (!std::isfinite(f[j])) {
            std::ostringstream message;
            message << "Grid element " << element << " (walking axis " << axis
                    << ") has a non-finite " << COMPONENT[j] << " coordinate ("
                    << f[j] << "); the geometry cannot be derived";
            throw GridError(message.str());
        }
        out[j] = f[j];
    }
}

/// Largest magnitude any float in an axis's geometry takes, which sets its noise
/// floor.
///
/// Both endpoints matter, because the grid centre a node coordinate is stored
/// relative to lies between them. So does the cell edge, n * spacing: it is
/// stored as a float too, and for a grid straddling the Cartesian origin it is
/// the largest of the three -- up to 2n / (n - 1) times either endpoint. Scaling
/// by the endpoints alone would then understate the noise fourfold on a two-node
/// axis. The result is never zero: it is at least the cell edge, and
/// get_grid_params has already rejected a non-positive interval.
double AxisMagnitude(const double origin, const unsigned int dim, const double spacing) {
    return std::max(std::max(std::abs(origin), std::abs(origin + (dim - 1u) * spacing)),
                    dim * spacing);
}

/// True when a fractional index lies within its axis's node span, allowing for
/// the float noise on the span's own endpoints.
///
/// The comparison runs in Angstroms rather than in fractional units: the slack is
/// naturally a Cartesian quantity, and multiplying the index by the interval
/// tests the same thing as dividing the slack by it while keeping a division off
/// the library's hottest path. The slack is a property of the grid alone, so a
/// query far from the grid does not widen the span it is tested against.
bool WithinAxisSpan(const double f, const double origin, const unsigned int dim,
                    const double spacing) {
    const double tol =
        NODE_SPAN_ROUNDINGS * FLOAT_HALF_ULP * AxisMagnitude(origin, dim, spacing);
    return f * spacing >= -tol && (f - (dim - 1u)) * spacing <= tol;
}

/// The containment predicate, taking an already-computed fractional index.
///
/// The spec requires `interpolate_density` to share `grid_contains`'s predicate
/// rather than restate it, so the prechecks and the interpolator cannot drift
/// apart. A direct call would not serve: `interpolate_density_at` (Task 2, same
/// translation unit) needs the fractional index itself for the blend, and
/// calling `grid_contains(gp, x, y, z)` would recompute `grid_fractional_index`
/// a second time in the hottest loop in the library. Taking the index instead
/// gives both callers one predicate and one index computation each.
///
/// The range test is written positively, which is what rejects a non-finite
/// index: every comparison against NaN is false, so NaN reports "outside"
/// without help, and so do both infinities. Written the other way round, as
/// `!(f < 0 || f > n - 1)`, NaN would pass.
///
/// The isfinite conjuncts are therefore redundant today -- deleting them leaves
/// every test in this suite passing, which is how that was established rather
/// than assumed. They are kept as a guard on the two ways this could stop being
/// true: a negated rewrite of the range test, and a slack derived from the query
/// instead of from the grid, which would make `f <= (n - 1) + inf` hold. Neither
/// is hypothetical enough to be worth a wrong answer for an infinite coordinate,
/// and the cost is three predictable branches outside the blend.
bool ContainsFractionalIndex(const GridParams& gp,
                             const double fx, const double fy, const double fz) {
    return std::isfinite(fx) && std::isfinite(fy) && std::isfinite(fz) &&
           WithinAxisSpan(fx, gp.x_origin, gp.x_dim, gp.x_spacing) &&
           WithinAxisSpan(fy, gp.y_origin, gp.y_dim, gp.y_spacing) &&
           WithinAxisSpan(fz, gp.z_origin, gp.z_dim, gp.z_spacing);
}

/// Trilinear blend of the eight corners named by per-axis element offsets.
///
/// The offsets arrive pre-multiplied by their axis stride so that one blend
/// serves both interpolators: the periodic path's upper corner may wrap to
/// element 0 of its axis, where the plain path's is always the lower corner plus
/// one stride. @p t holds the weight of the upper corner on each axis.
double BlendTrilinear(const float* values,
                      const unsigned int lo[3], const unsigned int hi[3],
                      const double t[3]) {
    const double c000 = values[lo[2] + lo[1] + lo[0]];
    const double c100 = values[lo[2] + lo[1] + hi[0]];
    const double c010 = values[lo[2] + hi[1] + lo[0]];
    const double c110 = values[lo[2] + hi[1] + hi[0]];
    const double c001 = values[hi[2] + lo[1] + lo[0]];
    const double c101 = values[hi[2] + lo[1] + hi[0]];
    const double c011 = values[hi[2] + hi[1] + lo[0]];
    const double c111 = values[hi[2] + hi[1] + hi[0]];

    const double c00 = c000 * (1.0 - t[0]) + c100 * t[0];
    const double c10 = c010 * (1.0 - t[0]) + c110 * t[0];
    const double c01 = c001 * (1.0 - t[0]) + c101 * t[0];
    const double c11 = c011 * (1.0 - t[0]) + c111 * t[0];

    const double c0 = c00 * (1.0 - t[1]) + c10 * t[1];
    const double c1 = c01 * (1.0 - t[1]) + c11 * t[1];

    return c0 * (1.0 - t[2]) + c1 * t[2];
}

/// True when two quantities agree to within @p tol scaled by their magnitude.
///
/// The operands descend from OpenEye's float grid coordinates, where one ulp
/// exceeds an absolute 1e-6 above roughly 8.4 Angstroms — so an absolute
/// comparison is exact float equality on every real crystallographic cell.
/// The unit floor keeps the tolerance absolute for sub-Angstrom quantities
/// such as node spacings, and gives the degree-valued cell angles a scale
/// their own magnitude supplies.
bool NearlyEqual(const double p, const double q, const double tol) {
    return std::abs(p - q) <= tol * std::max({1.0, std::abs(p), std::abs(q)});
}

}  // namespace

GridParams get_grid_params(const OESystem::OESkewGrid& grid) {
    static const char* const AXIS[3] = {"x", "y", "z"};
    // Off-axis drift above this over an axis's full span means the sampling is
    // not axis-aligned, so no per-axis spacing describes it.
    constexpr double MAX_OFF_AXIS_LEAK = 1e-4;

    const unsigned int n[3] = {grid.GetXDim(), grid.GetYDim(), grid.GetZDim()};

    // Check 1 runs over all three axes before any walk, so a grid that fails
    // both check 1 and check 5 reports check 1.
    for (int i = 0; i < 3; ++i) {
        if (n[i] < 2) {
            std::ostringstream message;
            message << "Grid axis " << AXIS[i] << " has dimension " << n[i]
                    << "; deriving a node interval needs at least 2 nodes on every axis";
            throw GridError(message.str());
        }
    }

    // Elements linearize x-fastest: el = iz*nx*ny + iy*nx + ix.
    const unsigned int step[3] = {1u, n[0], n[0] * n[1]};

    double origin[3];
    ReadNodeCoord(grid, 0u, AXIS[0], origin);

    double spacing[3];
    for (int i = 0; i < 3; ++i) {
        double node[3];
        ReadNodeCoord(grid, (n[i] - 1u) * step[i], AXIS[i], node);

        for (int j = 0; j < 3; ++j) {
            if (j == i) continue;
            const double leak = std::abs(node[j] - origin[j]);
            if (leak > MAX_OFF_AXIS_LEAK) {
                std::ostringstream message;
                message << "Grid axis " << AXIS[i] << " leaks " << leak
                        << " A into " << AXIS[j] << " over its full span (limit "
                        << MAX_OFF_AXIS_LEAK << " A); maptitude requires axis-aligned sampling";
                throw CellError(message.str());
            }
        }

        // Averaging over every interval on the axis, rather than measuring one
        // adjacent step, divides ElementToSpatialCoord's float noise by n_i - 1.
        // On 1d26 that is the difference between 0.902082443 (single step) and
        // 0.902083317 (full walk) against a true 0.902083333.
        spacing[i] = (node[i] - origin[i]) / (n[i] - 1u);

        if (!std::isfinite(spacing[i]) || spacing[i] <= 0.0) {
            std::ostringstream message;
            message << "Grid axis " << AXIS[i] << " has a derived node interval of "
                    << spacing[i] << " A; the interval must be finite and positive";
            throw GridError(message.str());
        }
    }

    return GridParams{origin[0], origin[1], origin[2],
                      n[0],      n[1],      n[2],
                      spacing[0], spacing[1], spacing[2]};
}

UnitCellParams get_unit_cell(const OESystem::OESkewGrid& grid) {
    float a = 0.0f, b = 0.0f, c = 0.0f, alpha = 0.0f, beta = 0.0f, gamma = 0.0f;
    if (!grid.HasUnitCell() || !grid.GetUnitCell(a, b, c, alpha, beta, gamma)) {
        throw CellError("Grid has no unit cell; its cell parameters cannot be read");
    }
    return UnitCellParams{a, b, c, alpha, beta, gamma};
}

void grid_node_origin(const GridParams& gp, double& x, double& y, double& z) {
    x = gp.x_origin;
    y = gp.y_origin;
    z = gp.z_origin;
}

void grid_fractional_index(const GridParams& gp,
                           const double x, const double y, const double z,
                           double& fx, double& fy, double& fz) {
    fx = (x - gp.x_origin) / gp.x_spacing;
    fy = (y - gp.y_origin) / gp.y_spacing;
    fz = (z - gp.z_origin) / gp.z_spacing;
}

bool grid_contains(const GridParams& gp,
                   const double x, const double y, const double z) {
    double fx = 0.0, fy = 0.0, fz = 0.0;
    grid_fractional_index(gp, x, y, z, fx, fy, fz);
    return ContainsFractionalIndex(gp, fx, fy, fz);
}

double interpolate_density_at(const GridParams& gp, const float* values,
                              const double x, const double y, const double z,
                              const double default_value) {
    double f[3] = {0.0, 0.0, 0.0};
    grid_fractional_index(gp, x, y, z, f[0], f[1], f[2]);

    // Steps 2 and 3 of the algorithm, through the one predicate grid_contains
    // also calls, so the containment prechecks and the interpolator cannot drift
    // apart. Passing the index rather than the coordinate keeps this to one
    // grid_fractional_index per point in the batch loops.
    if (!ContainsFractionalIndex(gp, f[0], f[1], f[2])) {
        return default_value;
    }

    // Clamping the base index keeps a point exactly on the far face inside the
    // last cell with t == 1.0, rather than indexing one node past the end. Both
    // clamps carry load now that the containment check admits a boundary
    // tolerance: floor() puts the corner one cell below the array just below
    // f == 0 and one cell above the last full cell just above f == n - 1.
    // Clamping the weight to match stops the blend extrapolating that tolerance
    // back out past the edge node.
    const unsigned int n[3] = {gp.x_dim, gp.y_dim, gp.z_dim};
    const unsigned int stride[3] = {1u, gp.x_dim, gp.x_dim * gp.y_dim};
    unsigned int lo[3], hi[3];
    double t[3];
    for (int i = 0; i < 3; ++i) {
        const double base = std::floor(f[i]);
        const double clamped = std::min(std::max(base, 0.0),
                                        static_cast<double>(n[i] - 2u));
        lo[i] = static_cast<unsigned int>(clamped) * stride[i];
        hi[i] = lo[i] + stride[i];
        t[i] = std::min(std::max(f[i] - clamped, 0.0), 1.0);
    }

    return BlendTrilinear(values, lo, hi, t);
}

bool same_grid_geometry(const OESystem::OESkewGrid& lhs,
                        const OESystem::OESkewGrid& rhs,
                        const double tol) {
    const GridParams a = get_grid_params(lhs);
    const GridParams b = get_grid_params(rhs);

    if (a.x_dim != b.x_dim || a.y_dim != b.y_dim || a.z_dim != b.z_dim) return false;
    if (!NearlyEqual(a.x_spacing, b.x_spacing, tol)) return false;
    if (!NearlyEqual(a.y_spacing, b.y_spacing, tol)) return false;
    if (!NearlyEqual(a.z_spacing, b.z_spacing, tol)) return false;
    if (!NearlyEqual(a.x_origin, b.x_origin, tol)) return false;
    if (!NearlyEqual(a.y_origin, b.y_origin, tol)) return false;
    if (!NearlyEqual(a.z_origin, b.z_origin, tol)) return false;

    if (lhs.HasUnitCell() != rhs.HasUnitCell()) return false;
    if (lhs.HasUnitCell()) {
        const UnitCellParams ca = get_unit_cell(lhs);
        const UnitCellParams cb = get_unit_cell(rhs);
        if (!NearlyEqual(ca.a, cb.a, tol)) return false;
        if (!NearlyEqual(ca.b, cb.b, tol)) return false;
        if (!NearlyEqual(ca.c, cb.c, tol)) return false;
        if (!NearlyEqual(ca.alpha, cb.alpha, tol)) return false;
        if (!NearlyEqual(ca.beta, cb.beta, tol)) return false;
        if (!NearlyEqual(ca.gamma, cb.gamma, tol)) return false;
    }
    return true;
}

std::vector<double> grid_to_vector(const OESystem::OESkewGrid& grid) {
    const unsigned int size = grid.GetSize();
    const float* values = grid.GetValues();
    std::vector<double> result(size);
    for (unsigned int i = 0; i < size; ++i) {
        result[i] = values[i];
    }
    return result;
}

void vector_to_grid(const std::vector<double>& values, OESystem::OESkewGrid& grid) {
    const unsigned int size = grid.GetSize();
    if (values.size() != size) {
        std::ostringstream message;
        message << "vector_to_grid needs exactly " << size << " values for this grid, got "
                << values.size() << "; a short vector previously left the tail of the grid "
                   "holding stale density";
        throw GridError(message.str());
    }
    float* out = grid.GetValues();
    for (unsigned int i = 0; i < size; ++i) {
        out[i] = static_cast<float>(values[i]);
    }
}

double interpolate_density(const OESystem::OESkewGrid& grid,
                           const double x, const double y, const double z,
                           const double default_value) {
    const GridParams gp = get_grid_params(grid);
    return interpolate_density_at(gp, grid.GetValues(), x, y, z, default_value);
}

std::vector<double> interpolate_density_batch(
    const OESystem::OESkewGrid& grid,
    const std::vector<double>& points,
    const size_t num_points,
    const double default_value) {
    const GridParams gp = get_grid_params(grid);
    const float* values = grid.GetValues();
    std::vector<double> result(num_points);
    for (size_t i = 0; i < num_points; ++i) {
        result[i] = interpolate_density_at(
            gp, values, points[i * 3], points[i * 3 + 1], points[i * 3 + 2],
            default_value);
    }
    return result;
}

void require_commensurate_cell(const GridParams& gp, const double cell_a,
                               const double cell_b, const double cell_c) {
    static const char* const EDGE[3] = {"a", "b", "c"};
    static const char* const AXIS[3] = {"x", "y", "z"};
    const double given[3] = {cell_a, cell_b, cell_c};
    const unsigned int n[3] = {gp.x_dim, gp.y_dim, gp.z_dim};
    const double spacing[3] = {gp.x_spacing, gp.y_spacing, gp.z_spacing};
    const double origin[3] = {gp.x_origin, gp.y_origin, gp.z_origin};

    for (int i = 0; i < 3; ++i) {
        const double extent = n[i] * spacing[i];
        const double tol = CELL_EXTENT_ROUNDINGS * FLOAT_HALF_ULP *
                           AxisMagnitude(origin[i], n[i], spacing[i]);
        // Negated so a non-finite edge lands here rather than passing a
        // comparison it cannot satisfy either way.
        if (!(std::abs(given[i] - extent) <= tol)) {
            std::ostringstream message;
            message << "Periodic interpolation needs cell edge " << EDGE[i] << " to equal the "
                       "extent the grid samples along " << AXIS[i] << ": got " << given[i]
                    << " A against " << extent << " A (" << n[i] << " nodes at "
                    << spacing[i] << " A), a difference of " << std::abs(given[i] - extent)
                    << " A against a " << tol << " A allowance for float node coordinates; "
                       "the grid does not tile that cell";
            throw CellError(message.str());
        }
    }
}

namespace {

/// The periodic blend, on a cell already checked commensurate.
///
/// Split out so the batch entry point can check the cell once and then run this
/// per point.
double WrapAndBlendPeriodic(const GridParams& gp, const float* values,
                            const double x, const double y, const double z,
                            const double default_value) {
    double f[3] = {0.0, 0.0, 0.0};
    grid_fractional_index(gp, x, y, z, f[0], f[1], f[2]);

    // The only remaining way out: fmod of a non-finite index is NaN, which would
    // index the array with garbage. What is tested is the fractional index, not
    // the coordinate -- a finite coordinate large enough that dividing it by the
    // spacing overflows lands here too. Every point that survives this has a
    // value, so it is the sole use of default_value on the periodic path.
    if (!std::isfinite(f[0]) || !std::isfinite(f[1]) || !std::isfinite(f[2])) {
        return default_value;
    }

    const unsigned int n[3] = {gp.x_dim, gp.y_dim, gp.z_dim};
    const unsigned int stride[3] = {1u, gp.x_dim, gp.x_dim * gp.y_dim};
    unsigned int lo[3], hi[3];
    double t[3];
    for (int i = 0; i < 3; ++i) {
        // Wrapping the fractional index modulo the node count, rather than the
        // Cartesian coordinate modulo the cell edge, makes the reduction itself
        // exact: fmod is exact by IEEE 754, and the period is the integer n_i.
        // The wrapped index therefore carries only the error already in f[i], and
        // a point a hundred cells out is placed no less accurately than one just
        // past the edge. Reducing the coordinate instead would subtract a rounded
        // multiple of a rounded cell edge, and that error would grow with the
        // number of cells crossed. Adding the period back for a negative remainder
        // can round the sum up to exactly n_i, which the modulus below absorbs.
        const double period = static_cast<double>(n[i]);
        double w = std::fmod(f[i], period);
        if (w < 0.0) w += period;

        const double base = std::floor(w);
        const unsigned int i0 = static_cast<unsigned int>(base) % n[i];
        lo[i] = i0 * stride[i];
        hi[i] = ((i0 + 1u) % n[i]) * stride[i];
        t[i] = std::min(std::max(w - base, 0.0), 1.0);
    }

    return BlendTrilinear(values, lo, hi, t);
}

}  // namespace

double interpolate_density_periodic_at(
    const GridParams& gp, const float* values,
    const double x, const double y, const double z,
    const double cell_a, const double cell_b, const double cell_c,
    const double default_value) {
    require_commensurate_cell(gp, cell_a, cell_b, cell_c);
    return WrapAndBlendPeriodic(gp, values, x, y, z, default_value);
}

double interpolate_density_periodic(
    const OESystem::OESkewGrid& grid,
    const double x, const double y, const double z,
    const double cell_a, const double cell_b, const double cell_c,
    const double default_value) {
    const GridParams gp = get_grid_params(grid);
    return interpolate_density_periodic_at(gp, grid.GetValues(), x, y, z,
                                           cell_a, cell_b, cell_c, default_value);
}

std::vector<double> interpolate_density_periodic_batch(
    const OESystem::OESkewGrid& grid,
    const std::vector<double>& points,
    const size_t num_points,
    const double cell_a, const double cell_b, const double cell_c,
    const double default_value) {
    const GridParams gp = get_grid_params(grid);
    // The cell is a property of the batch, not of a point in it, so the check
    // runs once here rather than num_points times inside the loop.
    require_commensurate_cell(gp, cell_a, cell_b, cell_c);

    const float* values = grid.GetValues();
    std::vector<double> result(num_points);
    for (size_t i = 0; i < num_points; ++i) {
        result[i] = WrapAndBlendPeriodic(
            gp, values, points[i * 3], points[i * 3 + 1], points[i * 3 + 2],
            default_value);
    }
    return result;
}

std::vector<unsigned int> get_atom_grid_points(
    const OESystem::OESkewGrid& grid,
    const double x, const double y, const double z, const double radius) {
    std::vector<unsigned int> result;

    const GridParams gp = get_grid_params(grid);
    const double r2 = radius * radius;
    const double origin[3] = {gp.x_origin, gp.y_origin, gp.z_origin};
    const double spacing[3] = {gp.x_spacing, gp.y_spacing, gp.z_spacing};
    const unsigned int dim[3] = {gp.x_dim, gp.y_dim, gp.z_dim};
    const double centre[3] = {x, y, z};

    int lo[3], hi[3];
    for (int i = 0; i < 3; ++i) {
        lo[i] = std::max(0, static_cast<int>(
            std::floor((centre[i] - radius - origin[i]) / spacing[i])));
        hi[i] = std::min(static_cast<int>(dim[i]) - 1, static_cast<int>(
            std::ceil((centre[i] + radius - origin[i]) / spacing[i])));
    }

    for (int ix = lo[0]; ix <= hi[0]; ++ix) {
        const double dx = origin[0] + ix * spacing[0] - x;
        for (int iy = lo[1]; iy <= hi[1]; ++iy) {
            const double dy = origin[1] + iy * spacing[1] - y;
            for (int iz = lo[2]; iz <= hi[2]; ++iz) {
                const double dz = origin[2] + iz * spacing[2] - z;
                if (dx * dx + dy * dy + dz * dz <= r2) {
                    result.push_back(static_cast<unsigned int>(
                        iz * static_cast<int>(dim[0]) * static_cast<int>(dim[1]) +
                        iy * static_cast<int>(dim[0]) + ix));
                }
            }
        }
    }

    return result;
}

}  // namespace Maptitude
