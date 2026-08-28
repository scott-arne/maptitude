#include "maptitude/DensityCalculator.h"
#include "maptitude/Error.h"
#include "maptitude/Grid.h"
#include "maptitude/ScatteringFactors.h"
#include "maptitude/UnitCell.h"

#include <oechem.h>
#include <oegrid.h>
#include <fftw3.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <memory>
#include <mutex>
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

struct MillerIndex {
    int h, k, l;
    double stol2;  // sin^2(theta)/lambda^2 = s^2/4
};

// ---- Miller index generation ----

static std::vector<MillerIndex> GenerateMillerIndices(
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
                const double s2 = (h * h) / (a * a) +
                                   (k * k) / (b * b) +
                                   (l * l) / (c * c);
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
    const OESystem::OEScalarGrid& out_template,
    OESystem::OEScalarGrid& out_grid) {
    const GridParams gp = get_grid_params(out_template);

    for (unsigned int iz = 0; iz < gp.z_dim; ++iz) {
        const double z = gp.z_origin + iz * gp.spacing;
        double fz = std::fmod(z / c, 1.0);
        if (fz < 0.0) fz += 1.0;
        const double gk = fz * nz;
        const int k0 = static_cast<int>(gk) % nz;
        const int k1 = (k0 + 1) % nz;
        const double dk = gk - static_cast<int>(gk);

        for (unsigned int iy = 0; iy < gp.y_dim; ++iy) {
            const double y = gp.y_origin + iy * gp.spacing;
            double fy = std::fmod(y / b, 1.0);
            if (fy < 0.0) fy += 1.0;
            const double gj = fy * ny;
            const int j0 = static_cast<int>(gj) % ny;
            const int j1 = (j0 + 1) % ny;
            const double dj = gj - static_cast<int>(gj);

            for (unsigned int ix = 0; ix < gp.x_dim; ++ix) {
                const double x = gp.x_origin + ix * gp.spacing;
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
                out_grid[elem] = static_cast<float>(val);
            }
        }
    }
}

// ==== FFTW RAII wrappers ====

namespace {
// This mutex has internal linkage and is sufficient only while this file is
// the sole FFTW translation unit. A second FFTW TU must share this mutex
// through an internal header rather than copying this block.

/// FFTW buffers come from fftw_malloc, not operator new, so they need their own
/// deleter. RAII here guards the throw paths: the nine null checks added during
/// this conversion are themselves the throw sites, and they are reachable in
/// practice — under planner contention fftw_plan_dft_3d returns NULL.
struct FftwComplexDeleter {
    void operator()(fftw_complex* p) const noexcept {
        if (p != nullptr) {
            fftw_free(p);
        }
    }
};

using FftwComplexPtr = std::unique_ptr<fftw_complex, FftwComplexDeleter>;

/// fftw_plan is a pointer to an opaque fftw_plan_s.
struct FftwPlanDeleter {
    void operator()(fftw_plan_s* p) const noexcept;
};

using FftwPlanPtr = std::unique_ptr<fftw_plan_s, FftwPlanDeleter>;

/// fftw_plan_dft_3d and fftw_destroy_plan mutate global planner state and are
/// not thread-safe. fftw_execute on an already-created plan is, so this guards
/// only creation and destruction. A C++ caller invoking Calculate from multiple
/// threads reaches this directly; Python callers are currently serialized by the
/// GIL since the module is not built with SWIG threading.
std::mutex& FftwPlannerMutex() {
    static std::mutex mutex;
    return mutex;
}

void FftwPlanDeleter::operator()(fftw_plan_s* p) const noexcept {
    if (p != nullptr) {
        std::lock_guard<std::mutex> lock(FftwPlannerMutex());
        fftw_destroy_plan(p);
    }
}

/// Create a plan under the planner lock.
FftwPlanPtr MakePlan3d(int n0, int n1, int n2, fftw_complex* in, fftw_complex* out, int sign,
                       unsigned int flags) {
    std::lock_guard<std::mutex> lock(FftwPlannerMutex());
    return FftwPlanPtr(fftw_plan_dft_3d(n0, n1, n2, in, out, sign, flags));
}

}  // namespace

// ==== Main Calculate method ====

OESystem::OEScalarGrid* DensityCalculator::Calculate(
    OEChem::OEMolBase& mol,
    const OESystem::OEScalarGrid& obs_grid,
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
                << "; lower the bin count. Above the limit `n_scale_shells + 1` wraps the "
                   "unsigned addition that sizes the shell-edge table";
        throw GridError(message.str());
    }

    const UnitCell& cell = pimpl_->cell;
    const auto& symops = pimpl_->symops;
    const double a = cell.a, b = cell.b, c = cell.c;
    const double sp = obs_grid.GetSpacing();

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
    auto miller = GenerateMillerIndices(a, b, c, resolution);
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
    const int nx = static_cast<int>(std::round(a / sp));
    const int ny = static_cast<int>(std::round(b / sp));
    const int nz = static_cast<int>(std::round(c / sp));

    // A spacing at or above twice a cell edge rounds that dimension to zero, which
    // both sizes the FFT allocation at zero and makes the Miller-index wrap below a
    // division by zero -- undefined behavior, and SIGFPE on x86-64.
    if (nx < 1 || ny < 1 || nz < 1) {
        const double edge = (nx < 1) ? a : (ny < 1) ? b : c;
        const char* axis = (nx < 1) ? "a" : (ny < 1) ? "b" : "c";
        const int dim = (nx < 1) ? nx : (ny < 1) ? ny : nz;
        std::ostringstream message;
        message << "Grid spacing " << sp << " A is too coarse for cell edge " << axis << " = "
                << edge << " A: the FFT grid would be " << dim
                << " points along that axis. Use a spacing below half the shortest cell edge";
        throw GridError(message.str());
    }

    const size_t grid_size = static_cast<size_t>(nx) * ny * nz;
    FftwComplexPtr Fc_3d_owner(fftw_alloc_complex(grid_size));
    if (!Fc_3d_owner) {
        std::ostringstream msg;
        msg << "FFTW allocation failed for the calculated structure factors ("
            << grid_size << " fftw_complex elements)";
        throw GridError(msg.str());
    }
    fftw_complex* Fc_3d = Fc_3d_owner.get();
    std::fill(reinterpret_cast<double*>(Fc_3d),
              reinterpret_cast<double*>(Fc_3d) + 2 * grid_size, 0.0);

    for (size_t i = 0; i < n_refl; ++i) {
        const int hi = ((miller[i].h % nx) + nx) % nx;
        const int ki = ((miller[i].k % ny) + ny) % ny;
        const int li = ((miller[i].l % nz) + nz) % nz;
        const size_t flat = hi * ny * nz + ki * nz + li;
        Fc_3d[flat][0] += Fc_real[i];
        Fc_3d[flat][1] += Fc_imag[i];
    }

    // ----------------------------------------------------------------
    // Step 7: Flat bulk solvent correction (optional)
    // ----------------------------------------------------------------
    if (k_sol != 0.0) {
        auto sol_mask = BuildSolventMask(mol, mask, a, b, c, nx, ny, nz,
            symops.empty() ? std::vector<SymOp>{SymOp()} : symops);

        // FFT the solvent mask
        FftwComplexPtr mask_fft_owner(fftw_alloc_complex(grid_size));
        if (!mask_fft_owner) {
            std::ostringstream msg;
            msg << "FFTW allocation failed for the solvent mask FFT ("
                << grid_size << " fftw_complex elements)";
            throw GridError(msg.str());
        }
        fftw_complex* mask_fft = mask_fft_owner.get();
        FftwComplexPtr mask_in_owner(fftw_alloc_complex(grid_size));
        if (!mask_in_owner) {
            std::ostringstream msg;
            msg << "FFTW allocation failed for the solvent mask input ("
                << grid_size << " fftw_complex elements)";
            throw GridError(msg.str());
        }
        fftw_complex* mask_in = mask_in_owner.get();
        for (size_t i = 0; i < grid_size; ++i) {
            mask_in[i][0] = sol_mask[i];
            mask_in[i][1] = 0.0;
        }

        FftwPlanPtr mask_plan = MakePlan3d(
            nx, ny, nz, mask_in, mask_fft, FFTW_FORWARD, FFTW_ESTIMATE);
        if (!mask_plan) {
            throw GridError("FFTW planning failed for the solvent mask");
        }
        fftw_execute(mask_plan.get());
        mask_plan.reset();
        mask_in_owner.reset();

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
                    Fc_3d[flat][0] += correction * mask_fft[flat][0];
                    Fc_3d[flat][1] += correction * mask_fft[flat][1];
                }
            }
        }

        mask_fft_owner.reset();
    }

    // ----------------------------------------------------------------
    // Step 8: Inverse FFT -> real-space density
    // ----------------------------------------------------------------
    FftwComplexPtr rho_complex_owner(fftw_alloc_complex(grid_size));
    if (!rho_complex_owner) {
        std::ostringstream msg;
        msg << "FFTW allocation failed for the density map ("
            << grid_size << " fftw_complex elements)";
        throw GridError(msg.str());
    }
    fftw_complex* rho_complex = rho_complex_owner.get();
    FftwPlanPtr ifft_plan = MakePlan3d(
        nx, ny, nz, Fc_3d, rho_complex, FFTW_BACKWARD, FFTW_ESTIMATE);
    if (!ifft_plan) {
        throw GridError("FFTW planning failed for the structure-factor inverse FFT");
    }
    fftw_execute(ifft_plan.get());
    ifft_plan.reset();
    Fc_3d_owner.reset();

    const double V = a * b * c;
    std::vector<double> rho_3d(grid_size);
    const double scale_factor = static_cast<double>(nx * ny * nz) / V;
    for (size_t i = 0; i < grid_size; ++i) {
        rho_3d[i] = rho_complex[i][0] * scale_factor /
                     static_cast<double>(grid_size);
    }

    // ----------------------------------------------------------------
    // Step 9: Per-shell amplitude scaling (optional)
    // ----------------------------------------------------------------
    if (n_scale_shells > 1) {
        // FFT the calculated density
        FftwComplexPtr rho_in_owner(fftw_alloc_complex(grid_size));
        if (!rho_in_owner) {
            std::ostringstream msg;
            msg << "FFTW allocation failed for the calculated density input ("
                << grid_size << " fftw_complex elements)";
            throw GridError(msg.str());
        }
        fftw_complex* rho_in = rho_in_owner.get();
        FftwComplexPtr F_calc_owner(fftw_alloc_complex(grid_size));
        if (!F_calc_owner) {
            std::ostringstream msg;
            msg << "FFTW allocation failed for the calculated structure factor FFT ("
                << grid_size << " fftw_complex elements)";
            throw GridError(msg.str());
        }
        fftw_complex* F_calc = F_calc_owner.get();
        for (size_t i = 0; i < grid_size; ++i) {
            rho_in[i][0] = rho_3d[i];
            rho_in[i][1] = 0.0;
        }
        FftwPlanPtr fwd_plan = MakePlan3d(
            nx, ny, nz, rho_in, F_calc, FFTW_FORWARD, FFTW_ESTIMATE);
        if (!fwd_plan) {
            throw GridError("FFTW planning failed for the calculated-density forward FFT");
        }
        fftw_execute(fwd_plan.get());
        fwd_plan.reset();
        rho_in_owner.reset();

        // Sample observed density onto UC grid and FFT
        FftwComplexPtr obs_in_owner(fftw_alloc_complex(grid_size));
        if (!obs_in_owner) {
            std::ostringstream msg;
            msg << "FFTW allocation failed for the observed density input ("
                << grid_size << " fftw_complex elements)";
            throw GridError(msg.str());
        }
        fftw_complex* obs_in = obs_in_owner.get();
        FftwComplexPtr Fobs_3d_owner(fftw_alloc_complex(grid_size));
        if (!Fobs_3d_owner) {
            std::ostringstream msg;
            msg << "FFTW allocation failed for the observed structure factor FFT ("
                << grid_size << " fftw_complex elements)";
            throw GridError(msg.str());
        }
        fftw_complex* Fobs_3d = Fobs_3d_owner.get();

        for (int i = 0; i < nx; ++i) {
            const double x = obs_grid.GetXMin() + (static_cast<double>(i) / nx) * a;
            for (int j = 0; j < ny; ++j) {
                const double y = obs_grid.GetYMin() +
                                 (static_cast<double>(j) / ny) * b;
                for (int k = 0; k < nz; ++k) {
                    const double z = obs_grid.GetZMin() +
                                     (static_cast<double>(k) / nz) * c;
                    const size_t flat = i * ny * nz + j * nz + k;
                    obs_in[flat][0] = interpolate_density(
                        obs_grid, x, y, z, 0.0);
                    obs_in[flat][1] = 0.0;
                }
            }
        }

        FftwPlanPtr obs_fwd = MakePlan3d(
            nx, ny, nz, obs_in, Fobs_3d, FFTW_FORWARD, FFTW_ESTIMATE);
        if (!obs_fwd) {
            throw GridError("FFTW planning failed for the observed-density forward FFT");
        }
        fftw_execute(obs_fwd.get());
        obs_fwd.reset();
        obs_in_owner.reset();

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
                            const double fc_amp = std::sqrt(
                                F_calc[flat][0] * F_calc[flat][0] +
                                F_calc[flat][1] * F_calc[flat][1]);
                            const double fobs_amp = std::sqrt(
                                Fobs_3d[flat][0] * Fobs_3d[flat][0] +
                                Fobs_3d[flat][1] * Fobs_3d[flat][1]);
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
                                F_calc[flat][0] *= k_shell;
                                F_calc[flat][1] *= k_shell;
                            }
                        }
                    }
                }
            }
        }

        Fobs_3d_owner.reset();

        // Inverse FFT scaled Fc back to real space
        FftwComplexPtr rho_scaled_owner(fftw_alloc_complex(grid_size));
        if (!rho_scaled_owner) {
            std::ostringstream msg;
            msg << "FFTW allocation failed for the scaled density map ("
                << grid_size << " fftw_complex elements)";
            throw GridError(msg.str());
        }
        fftw_complex* rho_scaled = rho_scaled_owner.get();
        FftwPlanPtr scale_ifft = MakePlan3d(
            nx, ny, nz, F_calc, rho_scaled, FFTW_BACKWARD, FFTW_ESTIMATE);
        if (!scale_ifft) {
            throw GridError("FFTW planning failed for the scaled structure-factor inverse FFT");
        }
        fftw_execute(scale_ifft.get());
        scale_ifft.reset();
        F_calc_owner.reset();

        for (size_t i = 0; i < grid_size; ++i) {
            rho_3d[i] = rho_scaled[i][0] / static_cast<double>(grid_size);
        }
        rho_scaled_owner.reset();
    }

    rho_complex_owner.reset();

    // ----------------------------------------------------------------
    // Step 10: Trilinear interpolation onto output grid
    // ----------------------------------------------------------------
    auto* calc_grid = new OESystem::OEScalarGrid(obs_grid);
    InterpolateUCToGrid(rho_3d.data(), nx, ny, nz, a, b, c,
                        obs_grid, *calc_grid);

    // ----------------------------------------------------------------
    // Step 11: Global linear scaling
    // ----------------------------------------------------------------
    float scale_coords[3];
    double sum_obs_calc = 0.0, sum_calc2 = 0.0;

    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() == 1) continue;
        if (mask && !(*mask)(*atom)) continue;

        mol.GetCoords(&(*atom), scale_coords);
        const float fx = scale_coords[0], fy = scale_coords[1], fz = scale_coords[2];

        if (obs_grid.IsInGrid(fx, fy, fz)) {
            const double obs_val = OESystem::OEFloatGridLinearInterpolate(
                obs_grid, fx, fy, fz, 0.0f);
            const double calc_val = OESystem::OEFloatGridLinearInterpolate(
                *calc_grid, fx, fy, fz, 0.0f);
            sum_obs_calc += obs_val * calc_val;
            sum_calc2 += calc_val * calc_val;
        }
    }

    if (sum_calc2 > 0.0) {
        const double k_scale = sum_obs_calc / sum_calc2;
        const unsigned int grid_sz = calc_grid->GetSize();
        for (unsigned int i = 0; i < grid_sz; ++i) {
            (*calc_grid)[i] = static_cast<float>((*calc_grid)[i] * k_scale);
        }
    }

    return calc_grid;
}

}  // namespace Maptitude
