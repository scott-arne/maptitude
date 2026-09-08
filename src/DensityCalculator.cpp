#include "maptitude/DensityCalculator.h"
#include "maptitude/detail/ScoringHelpers.h"
#include "maptitude/Error.h"
#include "maptitude/Grid.h"
#include "maptitude/ScatteringFactors.h"
#include "maptitude/UnitCell.h"

#include <oechem.h>
#include <oegrid.h>
#include <pocketfft_hdronly.h>

#include <cmath>
#include <complex>
#include <cstddef>
#include <memory>
#include <sstream>
#include <vector>

#ifdef MAPTITUDE_USE_OPENMP
#include <omp.h>
#endif

namespace Maptitude {

static constexpr double TWO_PI = 6.283185307179586;
static constexpr double DEFAULT_BFACTOR = 20.0;
static constexpr double PROBE_RADIUS = 1.4;
static constexpr double DEG_TO_RAD = 3.14159265358979323846 / 180.0;

/// The structure-factor pipeline computes 1/d^2 as
/// (h/a)^2 + (k/b)^2 + (l/c)^2 in four places, which is only correct for an
/// orthorhombic lattice. Non-orthorhombic cells would return a plausible wrong
/// answer; reject them until general lattice support lands.
static void RequireOrthorhombic(const UnitCell& cell) {
    constexpr double COSINE_TOLERANCE = 1e-9;
    const double cosines[3] = {std::cos(cell.alpha * DEG_TO_RAD), std::cos(cell.beta * DEG_TO_RAD),
                               std::cos(cell.gamma * DEG_TO_RAD)};
    const char* names[3] = {"alpha", "beta", "gamma"};
    for (int i = 0; i < 3; ++i) {
        if (std::fabs(cosines[i]) >= COSINE_TOLERANCE) {
            std::ostringstream message;
            message << "DensityCalculator supports orthorhombic cells only; angle " << names[i]
                    << " is " << (i == 0 ? cell.alpha : (i == 1 ? cell.beta : cell.gamma))
                    << " degrees";
            throw CellError(message.str());
        }
    }
}

// ---- Impl ----

struct DensityCalculator::Impl {
    UnitCell cell;
    std::vector<SymOp> symops;

    Impl(const UnitCell& cell, const std::vector<SymOp>& symops)
        : cell(cell), symops(symops) {}
};

DensityCalculator::DensityCalculator(const UnitCell& cell,
                                     const std::vector<SymOp>& symops)
    : pimpl_(std::make_unique<Impl>(cell, symops)) {
    // The members of UnitCell are public and mutable, so a cell can be
    // invalidated after construction. Re-check at the consumption point.
    validate_cell(pimpl_->cell);
    RequireOrthorhombic(pimpl_->cell);
}

DensityCalculator::~DensityCalculator() = default;

DensityCalculator::DensityCalculator(DensityCalculator&&) noexcept = default;
DensityCalculator& DensityCalculator::operator=(DensityCalculator&&) noexcept = default;

// ---- Internal structures ----

struct AtomData {
    double frac_x, frac_y, frac_z;
    double bfactor;
    int type_index;  // index into unique scattering factor types
};

using detail::MillerIndex;

// ---- Miller index generation ----

std::vector<MillerIndex> detail::GenerateMillerIndices(
    const double a, const double b, const double c, const double resolution) {
    const double s_max = 1.0 / resolution;
    const double s_max2 = s_max * s_max;
    const double h_extent = std::ceil(a * s_max);
    const double k_extent = std::ceil(b * s_max);
    const double l_extent = std::ceil(c * s_max);

    // Bound the loop volume before narrowing the extents to int. Both steps need this
    // guard: the triple loop below is O((2a/resolution)^3) and does not finish for a
    // small enough resolution, and an extent past INT_MAX makes the narrowing itself
    // undefined behavior. The product is formed in double, which saturates to infinity
    // instead of wrapping, so the comparison holds however extreme the request is.
    // NaN cannot arise -- the resolution is checked finite and positive on entry to
    // Calculate and the cell edges are validated at construction.
    const double box_points =
        (2.0 * h_extent + 1.0) * (2.0 * k_extent + 1.0) * (2.0 * l_extent + 1.0);
    if (box_points > MAX_MILLER_BOX_POINTS) {
        std::ostringstream message;
        message << "Resolution " << resolution << " A over cell edges a = " << a << " A, b = "
                << b << " A, c = " << c << " A needs a Miller-index box of " << box_points
                << " points, over the " << MAX_MILLER_BOX_POINTS
                << " limit (MAX_MILLER_BOX_POINTS); raise the resolution or use a smaller cell";
        throw GridError(message.str());
    }

    const int h_max = static_cast<int>(h_extent);
    const int k_max = static_cast<int>(k_extent);
    const int l_max = static_cast<int>(l_extent);

    std::vector<MillerIndex> indices;
    for (int h = -h_max; h <= h_max; ++h) {
        for (int k = -k_max; k <= k_max; ++k) {
            for (int l = -l_max; l <= l_max; ++l) {
                if (h == 0 && k == 0 && l == 0) continue;
                // Square in double, not int. The box bound above is on the product of
                // the three extents, so a cell with one long edge and two short ones
                // reaches a single extent past floor(sqrt(INT_MAX)) = 46340 while the
                // box stays far under the limit. `h * h` in int then overflowed, and the
                // wrapped negative s2 passed the test below, admitting reflections from
                // outside the requested shell. Widening the arithmetic rather than
                // bounding each axis keeps the anisotropic cells that are legitimate.
                //
                // Exact for every input the box bound admits, so no reflection that was
                // already in the shell moves: the extents are under 1.2e7, whose squares
                // are well inside the 2^53 range where double represents every integer,
                // and the int product was converted to this same double before dividing.
                const double hd = h, kd = k, ld = l;
                const double s2 = (hd * hd) / (a * a) +
                                   (kd * kd) / (b * b) +
                                   (ld * ld) / (c * c);
                if (s2 <= s_max2) {
                    indices.push_back({h, k, l, s2 / 4.0});
                }
            }
        }
    }
    return indices;
}

// ---- Build solvent mask ----

static std::vector<double> BuildSolventMask(
    OEChem::OEMolBase& mol,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask,
    const double a, const double b, const double c,
    const int nx, const int ny, const int nz,
    const std::vector<SymOp>& symops) {
    std::vector<double> sol_mask(nx * ny * nz, 1.0);

    // Fractional grid coordinates
    std::vector<double> fi(nx), fj(ny), fk(nz);
    for (int i = 0; i < nx; ++i) fi[i] = static_cast<double>(i) / nx;
    for (int j = 0; j < ny; ++j) fj[j] = static_cast<double>(j) / ny;
    for (int k = 0; k < nz; ++k) fk[k] = static_cast<double>(k) / nz;

    float coords[3];
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() == 1) continue;
        if (mask && !(*mask)(*atom)) continue;

        mol.GetCoords(&(*atom), coords);
        const double xf = coords[0] / a;
        const double yf = coords[1] / b;
        const double zf = coords[2] / c;

        double vdw = atom->GetRadius();
        if (vdw <= 0.0) vdw = 1.7;
        const double r_excl = vdw + PROBE_RADIUS;
        const double r_excl2 = r_excl * r_excl;

        for (const auto& op : symops) {
            // Apply symmetry: r_sym = R * frac + t
            const double sx = op.R[0] * xf + op.R[1] * yf + op.R[2] * zf + op.t[0];
            const double sy = op.R[3] * xf + op.R[4] * yf + op.R[5] * zf + op.t[1];
            const double sz = op.R[6] * xf + op.R[7] * yf + op.R[8] * zf + op.t[2];

            for (int i = 0; i < nx; ++i) {
                double dx = fi[i] - sx;
                dx -= std::round(dx);
                const double dx_cart = dx * a;

                for (int j = 0; j < ny; ++j) {
                    double dy = fj[j] - sy;
                    dy -= std::round(dy);
                    const double dy_cart = dy * b;

                    const double dxy2 = dx_cart * dx_cart + dy_cart * dy_cart;
                    if (dxy2 > r_excl2) continue;

                    for (int k = 0; k < nz; ++k) {
                        double dz = fk[k] - sz;
                        dz -= std::round(dz);
                        const double dz_cart = dz * c;

                        const double dist2 = dxy2 + dz_cart * dz_cart;
                        if (dist2 <= r_excl2) {
                            sol_mask[i * ny * nz + j * nz + k] = 0.0;
                        }
                    }
                }
            }
        }
    }

    return sol_mask;
}

// ---- Trilinear interpolation from UC grid to output grid ----

static void InterpolateUCToGrid(
    const double* rho_3d, const int nx, const int ny, const int nz,
    const double a, const double b, const double c,
    const OESystem::OESkewGrid& out_template,
    OESystem::OESkewGrid& out_grid) {
    const GridParams gp = get_grid_params(out_template);
    float* out_values = out_grid.GetValues();

    for (unsigned int iz = 0; iz < gp.z_dim; ++iz) {
        const double z = gp.z_origin + iz * gp.z_spacing;
        double fz = std::fmod(z / c, 1.0);
        if (fz < 0.0) fz += 1.0;
        const double gk = fz * nz;
        const int k0 = static_cast<int>(gk) % nz;
        const int k1 = (k0 + 1) % nz;
        const double dk = gk - static_cast<int>(gk);

        for (unsigned int iy = 0; iy < gp.y_dim; ++iy) {
            const double y = gp.y_origin + iy * gp.y_spacing;
            double fy = std::fmod(y / b, 1.0);
            if (fy < 0.0) fy += 1.0;
            const double gj = fy * ny;
            const int j0 = static_cast<int>(gj) % ny;
            const int j1 = (j0 + 1) % ny;
            const double dj = gj - static_cast<int>(gj);

            for (unsigned int ix = 0; ix < gp.x_dim; ++ix) {
                const double x = gp.x_origin + ix * gp.x_spacing;
                double fx = std::fmod(x / a, 1.0);
                if (fx < 0.0) fx += 1.0;
                const double gi = fx * nx;
                const int i0 = static_cast<int>(gi) % nx;
                const int i1 = (i0 + 1) % nx;
                const double di = gi - static_cast<int>(gi);

                // Index: [i][j][k] = i * ny * nz + j * nz + k
                auto idx = [&](const int ii, const int jj, const int kk) {
                    return ii * ny * nz + jj * nz + kk;
                };

                const double val =
                    rho_3d[idx(i0, j0, k0)] * (1 - di) * (1 - dj) * (1 - dk) +
                    rho_3d[idx(i1, j0, k0)] * di * (1 - dj) * (1 - dk) +
                    rho_3d[idx(i0, j1, k0)] * (1 - di) * dj * (1 - dk) +
                    rho_3d[idx(i0, j0, k1)] * (1 - di) * (1 - dj) * dk +
                    rho_3d[idx(i1, j1, k0)] * di * dj * (1 - dk) +
                    rho_3d[idx(i1, j0, k1)] * di * (1 - dj) * dk +
                    rho_3d[idx(i0, j1, k1)] * (1 - di) * dj * dk +
                    rho_3d[idx(i1, j1, k1)] * di * dj * dk;

                const unsigned int elem = iz * gp.x_dim * gp.y_dim +
                                    iy * gp.x_dim + ix;
                out_values[elem] = static_cast<float>(val);
            }
        }
    }
}

// ==== FFT helpers ====

namespace {

/// The transforms run over std::complex<double> volumes in C order, so an
/// element's position matches how every loop in Calculate addresses it:
/// flat = (i*ny + j)*nz + k.
using ComplexVolume = std::vector<std::complex<double>>;

/// Allocate a zero-filled FFT volume, reporting exhaustion as a GridError.
/// MAX_FFT_GRID_POINTS admits grids whose working set runs to gigabytes, so a
/// failed allocation is an input-driven outcome the caller can act on rather
/// than a programming error, and it belongs in this library's error type
/// instead of escaping as std::bad_alloc.
ComplexVolume MakeVolume(size_t grid_size, const char* what) {
    try {
        return ComplexVolume(grid_size);
    } catch (const std::bad_alloc&) {
        std::ostringstream msg;
        msg << "allocation failed for " << what << " (" << grid_size
            << " complex elements)";
        throw GridError(msg.str());
    }
}

/// Return a volume's storage to the allocator. The peak-memory figure
/// documented on MAX_FFT_GRID_POINTS counts only the volumes live at once, so
/// releasing each as soon as it is consumed is load-bearing rather than tidy.
/// clear() would keep the allocation; only the swap gives it back.
void ReleaseVolume(ComplexVolume& volume) {
    ComplexVolume().swap(volume);
}

/// Run an out-of-place 3D complex-to-complex transform over a whole volume.
///
/// pocketfft::FORWARD is the exp(-2*pi*i*h*x) direction, density to structure
/// factors, and BACKWARD its exp(+2*pi*i*h*x) inverse. Neither is normalized:
/// the scale factor is 1, and each caller divides by the point count itself.
///
/// PocketFFT carries no planner state between calls (POCKETFFT_CACHE_SIZE is 0
/// by default), so there is nothing global here to serialize and concurrent
/// Calculate calls need no lock.
void Transform3d(int nx, int ny, int nz, const ComplexVolume& in, ComplexVolume& out,
                 bool forward) {
    // pocketfft strides are byte counts, not element counts. ptrdiff_t
    // arithmetic throughout: MAX_FFT_GRID_POINTS admits an outer stride past
    // 3e9 bytes, which int cannot hold.
    constexpr std::ptrdiff_t ELEMENT_BYTES = sizeof(ComplexVolume::value_type);
    const pocketfft::shape_t shape{static_cast<size_t>(nx), static_cast<size_t>(ny),
                                   static_cast<size_t>(nz)};
    const pocketfft::stride_t stride{static_cast<std::ptrdiff_t>(ny) * nz * ELEMENT_BYTES,
                                     static_cast<std::ptrdiff_t>(nz) * ELEMENT_BYTES,
                                     ELEMENT_BYTES};
    const pocketfft::shape_t axes{0, 1, 2};
    pocketfft::c2c(shape, stride, stride, axes, forward, in.data(), out.data(), 1.0);
}

}  // namespace

// ==== Main Calculate method ====

OESystem::OESkewGrid* DensityCalculator::Calculate(
    OEChem::OEMolBase& mol,
    const OESystem::OESkewGrid& obs_grid,
    double resolution,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask,
    double k_sol,
    double b_sol,
    bool include_h,
    unsigned int n_scale_shells) const {
    require_usable_resolution(resolution);
    if (n_scale_shells < 1 || n_scale_shells > MAX_SCALE_SHELLS) {
        std::ostringstream message;
        message << "n_scale_shells must be between 1 and " << MAX_SCALE_SHELLS
                << " (MAX_SCALE_SHELLS), got " << n_scale_shells
                << "; lower the bin count. The bins partition the resolution range and even a "
                   "0.5 A dataset has far fewer independent shells than the limit. The "
                   "shell-edge table holds n_scale_shells + 1 doubles, so a value near "
                   "UINT_MAX asks for tens of gigabytes, and at UINT_MAX itself the unsigned "
                   "addition wraps to zero";
        throw GridError(message.str());
    }

    // The bulk-solvent mask reads each atom's radius. A molecule fresh from a
    // file carries none, and the scorers assign Bondi radii in place, so
    // without this the same molecule gave one Fc before rscc had seen it and
    // another after; the mask's 1.7 A fallback is carbon's radius, wrong for
    // every other element.
    detail::assign_missing_radii(mol);

    const UnitCell& cell = pimpl_->cell;
    const auto& symops = pimpl_->symops;
    const double a = cell.a, b = cell.b, c = cell.c;

    // ----------------------------------------------------------------
    // Step 1: Prepare structure - extract atom data
    // ----------------------------------------------------------------
    float coords_buf[3];
    std::vector<AtomData> atoms_data;
    std::vector<std::pair<unsigned int, int>> unique_keys;  // (Z, charge)

    auto find_or_add_type = [&](unsigned int z, int charge) -> int {
        auto key = std::make_pair(z, charge);
        for (size_t i = 0; i < unique_keys.size(); ++i) {
            if (unique_keys[i] == key) return static_cast<int>(i);
        }
        unique_keys.push_back(key);
        return static_cast<int>(unique_keys.size() - 1);
    };

    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (!include_h && atom->GetAtomicNum() == 1) continue;
        if (mask && !(*mask)(*atom)) continue;

        const unsigned int z_num = atom->GetAtomicNum();
        const int charge = atom->GetFormalCharge();

        // Look up scattering factors (fall back to neutral)
        const CromerMannCoeffs* cm = get_scattering_factors(z_num, charge);
        if (!cm) continue;

        // Get effective key for type indexing
        const unsigned int eff_z = z_num;
        int eff_charge = charge;
        if (!get_scattering_factors(z_num, charge)) {
            eff_charge = 0;
        }

        mol.GetCoords(&(*atom), coords_buf);

        // Convert to fractional using unit cell
        const auto frac = cell.CartesianToFractional(
            coords_buf[0], coords_buf[1], coords_buf[2]);

        const OEChem::OEResidue res = OEChem::OEAtomGetResidue(&(*atom));
        double bfac = res.GetBFactor();
        if (bfac <= 0.0) bfac = DEFAULT_BFACTOR;

        const int type_idx = find_or_add_type(eff_z, eff_charge);
        atoms_data.push_back({frac[0], frac[1], frac[2], bfac, type_idx});
    }

    if (atoms_data.empty()) {
        throw StructureError("No scorable atoms with known scattering factors");
    }

    // ----------------------------------------------------------------
    // Step 2: Generate Miller indices within resolution sphere
    // ----------------------------------------------------------------
    auto miller = detail::GenerateMillerIndices(a, b, c, resolution);
    const size_t n_refl = miller.size();

    // ----------------------------------------------------------------
    // Step 3: Precompute scattering factor table per type per reflection
    // ----------------------------------------------------------------
    const size_t n_types = unique_keys.size();
    // f_s_table[refl_idx * n_types + type_idx]
    std::vector<double> f_s_table(n_refl * n_types);
    for (size_t ti = 0; ti < n_types; ++ti) {
        const auto* cm = get_scattering_factors(
            unique_keys[ti].first, unique_keys[ti].second);
        for (size_t ri = 0; ri < n_refl; ++ri) {
            f_s_table[ri * n_types + ti] = cm->Evaluate(miller[ri].stol2);
        }
    }

    // ----------------------------------------------------------------
    // Step 4: Expand atoms by symmetry operators
    // ----------------------------------------------------------------
    const size_t n_atoms_asu = atoms_data.size();
    size_t n_ops = symops.size();
    if (n_ops == 0) n_ops = 1;  // identity

    const size_t n_expanded = n_atoms_asu * n_ops;
    std::vector<double> frac_x(n_expanded), frac_y(n_expanded), frac_z(n_expanded);
    std::vector<double> bfacs(n_expanded);
    std::vector<int> type_indices(n_expanded);

    size_t idx = 0;
    for (const auto& ad : atoms_data) {
        if (symops.empty()) {
            frac_x[idx] = ad.frac_x;
            frac_y[idx] = ad.frac_y;
            frac_z[idx] = ad.frac_z;
            bfacs[idx] = ad.bfactor;
            type_indices[idx] = ad.type_index;
            ++idx;
        } else {
            for (const auto& op : symops) {
                const auto sym = op.Apply(ad.frac_x, ad.frac_y, ad.frac_z);
                frac_x[idx] = sym[0];
                frac_y[idx] = sym[1];
                frac_z[idx] = sym[2];
                bfacs[idx] = ad.bfactor;
                type_indices[idx] = ad.type_index;
                ++idx;
            }
        }
    }

    // ----------------------------------------------------------------
    // Step 5: Accumulate structure factors Fc (OpenMP parallel)
    // ----------------------------------------------------------------
    std::vector<double> Fc_real(n_refl, 0.0);
    std::vector<double> Fc_imag(n_refl, 0.0);

#ifdef MAPTITUDE_USE_OPENMP
    #pragma omp parallel for schedule(dynamic, 64)
#endif
    for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(n_refl); ++i) {
        const int h = miller[i].h;
        const int k = miller[i].k;
        const int l = miller[i].l;
        const double s2 = miller[i].stol2;
        double re = 0.0, im = 0.0;

        for (size_t j = 0; j < n_expanded; ++j) {
            const double f_s = f_s_table[static_cast<size_t>(i) * n_types + type_indices[j]];
            const double dw = std::exp(-bfacs[j] * s2);
            const double f_dw = f_s * dw;
            const double phase = -TWO_PI * (
                h * frac_x[j] + k * frac_y[j] + l * frac_z[j]);
            re += f_dw * std::cos(phase);
            im += f_dw * std::sin(phase);
        }

        Fc_real[i] = re;
        Fc_imag[i] = im;
    }

    // ----------------------------------------------------------------
    // Step 6: Scatter Fc into 3D FFT array
    // ----------------------------------------------------------------
    // The FFT grid samples the unit cell, so its counts come from the cell
    // edges and the map's per-axis intervals. Deriving them from one scalar
    // spacing resampled one of the three axes: on 1d26 the header says
    // 48x48x24 and a scalar spacing gives 48x48x27.

    // How far the map's sampling and the cell's may disagree before they are
    // taken to describe different samplings rather than the same one recorded
    // with rounding. Measured on 1d26, the disagreement is 2.5e-16 on a and
    // 1.2e-16 on b but 7.4e-08 on c; the worst of the three still sits four
    // orders of magnitude inside this limit.
    constexpr double SAMPLING_AGREEMENT_TOLERANCE = 1e-3;

    const GridParams obs_gp = get_grid_params(obs_grid);
    const double edges[3] = {a, b, c};
    const double intervals[3] = {obs_gp.x_spacing, obs_gp.y_spacing, obs_gp.z_spacing};
    const char* const EDGE_NAMES[3] = {"a", "b", "c"};
    int counts[3];
    for (int i = 0; i < 3; ++i) {
        // Bound the sample count before narrowing it. The cell edge is validated only as
        // finite and positive and the node interval only as finite and positive, so
        // nothing caps their quotient and a cast past INT_MAX is undefined behavior. The
        // quotient is formed in double, which saturates to infinity instead of wrapping,
        // so the comparison holds however extreme the pair is.
        const double samples = std::round(edges[i] / intervals[i]);
        if (!std::isfinite(samples) || samples > MAX_FFT_GRID_POINTS) {
            std::ostringstream message;
            message << "Cell edge " << EDGE_NAMES[i] << " = " << edges[i]
                    << " A at the map's node interval of " << intervals[i] << " A needs "
                    << samples << " FFT samples along that axis, over the "
                    << MAX_FFT_GRID_POINTS
                    << " limit (MAX_FFT_GRID_POINTS); use a map with coarser sampling or a "
                       "smaller cell";
            throw GridError(message.str());
        }
        counts[i] = static_cast<int>(samples);
        if (counts[i] < 1) continue;  // reported by the existing check below
        const double implied = edges[i] / counts[i];
        const double disagreement = std::abs(implied - intervals[i]) / intervals[i];
        if (disagreement > SAMPLING_AGREEMENT_TOLERANCE) {
            std::ostringstream message;
            message << "Cell edge " << EDGE_NAMES[i] << " = " << edges[i]
                    << " A does not divide into the map's node interval of "
                    << intervals[i] << " A: " << counts[i] << " samples imply "
                    << implied << " A, a relative disagreement of " << disagreement
                    << " (limit " << SAMPLING_AGREEMENT_TOLERANCE
                    << "). The map and the cell describe different samplings";
            throw GridError(message.str());
        }
    }
    const int nx = counts[0];
    const int ny = counts[1];
    const int nz = counts[2];

    // A node interval above twice a cell edge rounds that dimension to zero,
    // which both sizes the FFT allocation at zero and makes the Miller-index wrap
    // below a division by zero -- undefined behavior, and SIGFPE on x86-64. At
    // exactly twice the edge the ratio is 0.5 and std::round(0.5) is 1, so this
    // branch does not fire; that case is rejected by the divisibility check
    // above, whose relative disagreement of 0.5 is far past its tolerance.
    if (nx < 1 || ny < 1 || nz < 1) {
        const int failing = (nx < 1) ? 0 : (ny < 1) ? 1 : 2;
        std::ostringstream message;
        message << "Node interval " << intervals[failing] << " A is too coarse for cell edge "
                << EDGE_NAMES[failing] << " = " << edges[failing]
                << " A: the FFT grid would be " << counts[failing]
                << " points along that axis. Use a spacing below half the shortest cell edge";
        throw GridError(message.str());
    }

    // Three counts that each pass the per-axis bound can still multiply past it, and the
    // product is what indexes the FFT array: `hi * ny * nz` in the scatter loop below is
    // int arithmetic, so it overflows before the size_t element count does. Formed in
    // double for the same reason as the per-axis check, and evaluated after the per-axis
    // loop rather than inside it, because the three counts are only all known once that
    // loop has finished.
    const double total_samples = static_cast<double>(nx) * ny * nz;
    if (total_samples > MAX_FFT_GRID_POINTS) {
        std::ostringstream message;
        message << "Cell edges a = " << a << " A, b = " << b << " A and c = " << c
                << " A at the map's node intervals of " << intervals[0] << " A, "
                << intervals[1] << " A and " << intervals[2] << " A need an FFT grid of "
                << nx << " x " << ny << " x " << nz << " = " << total_samples
                << " points, over the " << MAX_FFT_GRID_POINTS
                << " limit (MAX_FFT_GRID_POINTS); use a map with coarser sampling or a "
                   "smaller cell";
        throw GridError(message.str());
    }

    const size_t grid_size = static_cast<size_t>(nx) * ny * nz;
    // The scatter below accumulates into this volume, so it has to start at
    // zero; MakeVolume value-initializes.
    ComplexVolume Fc_3d = MakeVolume(grid_size, "the calculated structure factors");

    for (size_t i = 0; i < n_refl; ++i) {
        const int hi = ((miller[i].h % nx) + nx) % nx;
        const int ki = ((miller[i].k % ny) + ny) % ny;
        const int li = ((miller[i].l % nz) + nz) % nz;
        const size_t flat = hi * ny * nz + ki * nz + li;
        Fc_3d[flat] += std::complex<double>(Fc_real[i], Fc_imag[i]);
    }

    // ----------------------------------------------------------------
    // Step 7: Flat bulk solvent correction (optional)
    // ----------------------------------------------------------------
    if (k_sol != 0.0) {
        auto sol_mask = BuildSolventMask(mol, mask, a, b, c, nx, ny, nz,
            symops.empty() ? std::vector<SymOp>{SymOp()} : symops);

        // FFT the solvent mask
        ComplexVolume mask_fft = MakeVolume(grid_size, "the solvent mask FFT");
        ComplexVolume mask_in = MakeVolume(grid_size, "the solvent mask input");
        for (size_t i = 0; i < grid_size; ++i) {
            mask_in[i] = sol_mask[i];
        }

        Transform3d(nx, ny, nz, mask_in, mask_fft, pocketfft::FORWARD);
        ReleaseVolume(mask_in);

        // Compute S^2 for each FFT grid point and apply correction
        for (int i = 0; i < nx; ++i) {
            const double hi = (i <= nx / 2) ? i : i - nx;
            for (int j = 0; j < ny; ++j) {
                const double ki = (j <= ny / 2) ? j : j - ny;
                for (int k = 0; k < nz; ++k) {
                    const double li = (k <= nz / 2) ? k : k - nz;
                    const double s2 = (hi / a) * (hi / a) +
                                      (ki / b) * (ki / b) +
                                      (li / c) * (li / c);
                    const double correction = k_sol * std::exp(-b_sol * s2 / 4.0);
                    const size_t flat = i * ny * nz + j * nz + k;
                    Fc_3d[flat] += correction * mask_fft[flat];
                }
            }
        }
    }

    // ----------------------------------------------------------------
    // Step 8: Inverse FFT -> real-space density
    // ----------------------------------------------------------------
    ComplexVolume rho_complex = MakeVolume(grid_size, "the density map");
    Transform3d(nx, ny, nz, Fc_3d, rho_complex, pocketfft::BACKWARD);
    ReleaseVolume(Fc_3d);

    const double V = a * b * c;
    std::vector<double> rho_3d(grid_size);
    const double scale_factor = static_cast<double>(nx * ny * nz) / V;
    for (size_t i = 0; i < grid_size; ++i) {
        rho_3d[i] = rho_complex[i].real() * scale_factor /
                     static_cast<double>(grid_size);
    }

    // ----------------------------------------------------------------
    // Step 9: Per-shell amplitude scaling (optional)
    // ----------------------------------------------------------------
    if (n_scale_shells > 1) {
        // FFT the calculated density
        ComplexVolume rho_in = MakeVolume(grid_size, "the calculated density input");
        ComplexVolume F_calc =
            MakeVolume(grid_size, "the calculated structure factor FFT");
        for (size_t i = 0; i < grid_size; ++i) {
            rho_in[i] = rho_3d[i];
        }
        Transform3d(nx, ny, nz, rho_in, F_calc, pocketfft::FORWARD);
        ReleaseVolume(rho_in);

        // Sample observed density onto UC grid and FFT
        ComplexVolume obs_in = MakeVolume(grid_size, "the observed density input");
        ComplexVolume Fobs_3d =
            MakeVolume(grid_size, "the observed structure factor FFT");

        const float* obs_values = obs_grid.GetValues();
        for (int i = 0; i < nx; ++i) {
            const double x = obs_gp.x_origin + (static_cast<double>(i) / nx) * a;
            for (int j = 0; j < ny; ++j) {
                const double y = obs_gp.y_origin +
                                 (static_cast<double>(j) / ny) * b;
                for (int k = 0; k < nz; ++k) {
                    const double z = obs_gp.z_origin +
                                     (static_cast<double>(k) / nz) * c;
                    const size_t flat = i * ny * nz + j * nz + k;
                    obs_in[flat] =
                        interpolate_density_at(obs_gp, obs_values, x, y, z, 0.0);
                }
            }
        }

        Transform3d(nx, ny, nz, obs_in, Fobs_3d, pocketfft::FORWARD);
        ReleaseVolume(obs_in);

        // Per-shell scaling
        const double s2_max = 1.0 / (resolution * resolution);
        std::vector<double> shell_edges(n_scale_shells + 1);
        for (unsigned int i = 0; i <= n_scale_shells; ++i) {
            shell_edges[i] = s2_max * i / n_scale_shells;
        }

        for (unsigned int shell = 0; shell < n_scale_shells; ++shell) {
            const double s2_lo = shell_edges[shell];
            const double s2_hi = shell_edges[shell + 1];

            double sum_fobs_fc = 0.0;
            double sum_fc2 = 0.0;

            for (int i = 0; i < nx; ++i) {
                const double hi = (i <= nx / 2) ? i : i - nx;
                for (int j = 0; j < ny; ++j) {
                    const double ki = (j <= ny / 2) ? j : j - ny;
                    for (int k = 0; k < nz; ++k) {
                        const double li = (k <= nz / 2) ? k : k - nz;
                        const double s2 = (hi / a) * (hi / a) +
                                          (ki / b) * (ki / b) +
                                          (li / c) * (li / c);

                        bool in_shell;
                        if (shell == 0) {
                            in_shell = (s2 > 0) && (s2 <= s2_hi);
                        } else {
                            in_shell = (s2 > s2_lo) && (s2 <= s2_hi);
                        }

                        if (in_shell) {
                            const size_t flat = i * ny * nz + j * nz + k;
                            // Spelled out rather than std::abs, which is
                            // hypot: the pinned shell values were measured on
                            // this arithmetic.
                            const double fc_amp = std::sqrt(
                                F_calc[flat].real() * F_calc[flat].real() +
                                F_calc[flat].imag() * F_calc[flat].imag());
                            const double fobs_amp = std::sqrt(
                                Fobs_3d[flat].real() * Fobs_3d[flat].real() +
                                Fobs_3d[flat].imag() * Fobs_3d[flat].imag());
                            sum_fobs_fc += fobs_amp * fc_amp;
                            sum_fc2 += fc_amp * fc_amp;
                        }
                    }
                }
            }

            if (sum_fc2 > 0.0) {
                const double k_shell = sum_fobs_fc / sum_fc2;
                for (int i = 0; i < nx; ++i) {
                    const double hi = (i <= nx / 2) ? i : i - nx;
                    for (int j = 0; j < ny; ++j) {
                        const double ki = (j <= ny / 2) ? j : j - ny;
                        for (int k = 0; k < nz; ++k) {
                            const double li = (k <= nz / 2) ? k : k - nz;
                            const double s2 = (hi / a) * (hi / a) +
                                              (ki / b) * (ki / b) +
                                              (li / c) * (li / c);
                            bool in_shell;
                            if (shell == 0) {
                                in_shell = (s2 > 0) && (s2 <= s2_hi);
                            } else {
                                in_shell = (s2 > s2_lo) && (s2 <= s2_hi);
                            }
                            if (in_shell) {
                                const size_t flat = i * ny * nz + j * nz + k;
                                F_calc[flat] *= k_shell;
                            }
                        }
                    }
                }
            }
        }

        ReleaseVolume(Fobs_3d);

        // Inverse FFT scaled Fc back to real space
        ComplexVolume rho_scaled = MakeVolume(grid_size, "the scaled density map");
        Transform3d(nx, ny, nz, F_calc, rho_scaled, pocketfft::BACKWARD);
        ReleaseVolume(F_calc);

        for (size_t i = 0; i < grid_size; ++i) {
            rho_3d[i] = rho_scaled[i].real() / static_cast<double>(grid_size);
        }
    }

    ReleaseVolume(rho_complex);

    // ----------------------------------------------------------------
    // Step 10: Trilinear interpolation onto output grid
    // ----------------------------------------------------------------
    auto* calc_grid = new OESystem::OESkewGrid(obs_grid);
    InterpolateUCToGrid(rho_3d.data(), nx, ny, nz, a, b, c,
                        obs_grid, *calc_grid);

    // ----------------------------------------------------------------
    // Step 11: Global linear scaling
    // ----------------------------------------------------------------
    float scale_coords[3];
    double sum_obs_calc = 0.0, sum_calc2 = 0.0;

    const GridParams calc_gp = get_grid_params(*calc_grid);
    const float* obs_values = obs_grid.GetValues();
    const float* calc_values = calc_grid->GetValues();

    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() == 1) continue;
        if (mask && !(*mask)(*atom)) continue;

        mol.GetCoords(&(*atom), scale_coords);
        const float fx = scale_coords[0], fy = scale_coords[1], fz = scale_coords[2];

        if (grid_contains(obs_gp, fx, fy, fz)) {
            const double obs_val = interpolate_density_at(obs_gp, obs_values, fx, fy, fz, 0.0);
            const double calc_val = interpolate_density_at(calc_gp, calc_values, fx, fy, fz, 0.0);
            sum_obs_calc += obs_val * calc_val;
            sum_calc2 += calc_val * calc_val;
        }
    }

    if (sum_calc2 > 0.0) {
        const double k_scale = sum_obs_calc / sum_calc2;
        float* out = calc_grid->GetValues();
        const unsigned int grid_sz = calc_grid->GetSize();
        for (unsigned int i = 0; i < grid_sz; ++i) {
            out[i] = static_cast<float>(out[i] * k_scale);
        }
    }

    return calc_grid;
}

}  // namespace Maptitude
