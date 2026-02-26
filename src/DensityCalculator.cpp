#include "maptitude/DensityCalculator.h"
#include "maptitude/Error.h"
#include "maptitude/Grid.h"
#include "maptitude/ScatteringFactors.h"

#include <oechem.h>
#include <oegrid.h>
#include <fftw3.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

#ifdef MAPTITUDE_USE_OPENMP
#include <omp.h>
#endif

namespace Maptitude {

static constexpr double TWO_PI = 6.283185307179586;
static constexpr double DEFAULT_BFACTOR = 20.0;
static constexpr double PROBE_RADIUS = 1.4;

// ---- Impl ----

struct DensityCalculator::Impl {
    UnitCell cell;
    std::vector<SymOp> symops;

    Impl(const UnitCell& cell, const std::vector<SymOp>& symops)
        : cell(cell), symops(symops) {}
};

DensityCalculator::DensityCalculator(const UnitCell& cell,
                                     const std::vector<SymOp>& symops)
    : pimpl_(std::make_unique<Impl>(cell, symops)) {}

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
    double a, double b, double c, double resolution) {
    double s_max = 1.0 / resolution;
    double s_max2 = s_max * s_max;
    int h_max = static_cast<int>(std::ceil(a * s_max));
    int k_max = static_cast<int>(std::ceil(b * s_max));
    int l_max = static_cast<int>(std::ceil(c * s_max));

    std::vector<MillerIndex> indices;
    for (int h = -h_max; h <= h_max; ++h) {
        for (int k = -k_max; k <= k_max; ++k) {
            for (int l = -l_max; l <= l_max; ++l) {
                if (h == 0 && k == 0 && l == 0) continue;
                double s2 = (h * h) / (a * a) +
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
    double a, double b, double c,
    int nx, int ny, int nz,
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
        double xf = coords[0] / a;
        double yf = coords[1] / b;
        double zf = coords[2] / c;

        double vdw = atom->GetRadius();
        if (vdw <= 0.0) vdw = 1.7;
        double r_excl = vdw + PROBE_RADIUS;
        double r_excl2 = r_excl * r_excl;

        for (const auto& op : symops) {
            // Apply symmetry: r_sym = R * frac + t
            double sx = op.R[0] * xf + op.R[1] * yf + op.R[2] * zf + op.t[0];
            double sy = op.R[3] * xf + op.R[4] * yf + op.R[5] * zf + op.t[1];
            double sz = op.R[6] * xf + op.R[7] * yf + op.R[8] * zf + op.t[2];

            for (int i = 0; i < nx; ++i) {
                double dx = fi[i] - sx;
                dx -= std::round(dx);
                double dx_cart = dx * a;

                for (int j = 0; j < ny; ++j) {
                    double dy = fj[j] - sy;
                    dy -= std::round(dy);
                    double dy_cart = dy * b;

                    double dxy2 = dx_cart * dx_cart + dy_cart * dy_cart;
                    if (dxy2 > r_excl2) continue;

                    for (int k = 0; k < nz; ++k) {
                        double dz = fk[k] - sz;
                        dz -= std::round(dz);
                        double dz_cart = dz * c;

                        double dist2 = dxy2 + dz_cart * dz_cart;
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
    const double* rho_3d, int nx, int ny, int nz,
    double a, double b, double c,
    const OESystem::OEScalarGrid& out_template,
    OESystem::OEScalarGrid& out_grid) {
    GridParams gp = GetGridParams(out_template);

    for (unsigned int iz = 0; iz < gp.z_dim; ++iz) {
        double z = gp.z_origin + iz * gp.spacing;
        double fz = std::fmod(z / c, 1.0);
        if (fz < 0.0) fz += 1.0;
        double gk = fz * nz;
        int k0 = static_cast<int>(gk) % nz;
        int k1 = (k0 + 1) % nz;
        double dk = gk - static_cast<int>(gk);

        for (unsigned int iy = 0; iy < gp.y_dim; ++iy) {
            double y = gp.y_origin + iy * gp.spacing;
            double fy = std::fmod(y / b, 1.0);
            if (fy < 0.0) fy += 1.0;
            double gj = fy * ny;
            int j0 = static_cast<int>(gj) % ny;
            int j1 = (j0 + 1) % ny;
            double dj = gj - static_cast<int>(gj);

            for (unsigned int ix = 0; ix < gp.x_dim; ++ix) {
                double x = gp.x_origin + ix * gp.spacing;
                double fx = std::fmod(x / a, 1.0);
                if (fx < 0.0) fx += 1.0;
                double gi = fx * nx;
                int i0 = static_cast<int>(gi) % nx;
                int i1 = (i0 + 1) % nx;
                double di = gi - static_cast<int>(gi);

                // Index: [i][j][k] = i * ny * nz + j * nz + k
                auto idx = [&](int ii, int jj, int kk) {
                    return ii * ny * nz + jj * nz + kk;
                };

                double val =
                    rho_3d[idx(i0, j0, k0)] * (1 - di) * (1 - dj) * (1 - dk) +
                    rho_3d[idx(i1, j0, k0)] * di * (1 - dj) * (1 - dk) +
                    rho_3d[idx(i0, j1, k0)] * (1 - di) * dj * (1 - dk) +
                    rho_3d[idx(i0, j0, k1)] * (1 - di) * (1 - dj) * dk +
                    rho_3d[idx(i1, j1, k0)] * di * dj * (1 - dk) +
                    rho_3d[idx(i1, j0, k1)] * di * (1 - dj) * dk +
                    rho_3d[idx(i0, j1, k1)] * (1 - di) * dj * dk +
                    rho_3d[idx(i1, j1, k1)] * di * dj * dk;

                unsigned int elem = iz * gp.x_dim * gp.y_dim +
                                    iy * gp.x_dim + ix;
                out_grid[elem] = static_cast<float>(val);
            }
        }
    }
}

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
    if (resolution <= 0.0) {
        throw GridError("Resolution must be positive");
    }
    if (n_scale_shells < 1) {
        throw GridError("n_scale_shells must be >= 1");
    }

    const UnitCell& cell = pimpl_->cell;
    const auto& symops = pimpl_->symops;
    double a = cell.a, b = cell.b, c = cell.c;
    double sp = obs_grid.GetSpacing();

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

        unsigned int z_num = atom->GetAtomicNum();
        int charge = atom->GetFormalCharge();

        // Look up scattering factors (fall back to neutral)
        const CromerMannCoeffs* cm = GetScatteringFactors(z_num, charge);
        if (!cm) continue;

        // Get effective key for type indexing
        unsigned int eff_z = z_num;
        int eff_charge = charge;
        if (!GetScatteringFactors(z_num, charge)) {
            eff_charge = 0;
        }

        mol.GetCoords(&(*atom), coords_buf);

        // Convert to fractional using unit cell
        auto frac = cell.CartesianToFractional(
            coords_buf[0], coords_buf[1], coords_buf[2]);

        OEChem::OEResidue res = OEChem::OEAtomGetResidue(&(*atom));
        double bfac = res.GetBFactor();
        if (bfac <= 0.0) bfac = DEFAULT_BFACTOR;

        int type_idx = find_or_add_type(eff_z, eff_charge);
        atoms_data.push_back({frac[0], frac[1], frac[2], bfac, type_idx});
    }

    if (atoms_data.empty()) {
        throw StructureError("No scorable atoms with known scattering factors");
    }

    // ----------------------------------------------------------------
    // Step 2: Generate Miller indices within resolution sphere
    // ----------------------------------------------------------------
    auto miller = GenerateMillerIndices(a, b, c, resolution);
    size_t n_refl = miller.size();

    // ----------------------------------------------------------------
    // Step 3: Precompute scattering factor table per type per reflection
    // ----------------------------------------------------------------
    size_t n_types = unique_keys.size();
    // f_s_table[refl_idx * n_types + type_idx]
    std::vector<double> f_s_table(n_refl * n_types);
    for (size_t ti = 0; ti < n_types; ++ti) {
        const auto* cm = GetScatteringFactors(
            unique_keys[ti].first, unique_keys[ti].second);
        for (size_t ri = 0; ri < n_refl; ++ri) {
            f_s_table[ri * n_types + ti] = cm->Evaluate(miller[ri].stol2);
        }
    }

    // ----------------------------------------------------------------
    // Step 4: Expand atoms by symmetry operators
    // ----------------------------------------------------------------
    size_t n_atoms_asu = atoms_data.size();
    size_t n_ops = symops.size();
    if (n_ops == 0) n_ops = 1;  // identity

    size_t n_expanded = n_atoms_asu * n_ops;
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
                auto sym = op.Apply(ad.frac_x, ad.frac_y, ad.frac_z);
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
    for (size_t i = 0; i < n_refl; ++i) {
        int h = miller[i].h;
        int k = miller[i].k;
        int l = miller[i].l;
        double s2 = miller[i].stol2;
        double re = 0.0, im = 0.0;

        for (size_t j = 0; j < n_expanded; ++j) {
            double f_s = f_s_table[i * n_types + type_indices[j]];
            double dw = std::exp(-bfacs[j] * s2);
            double f_dw = f_s * dw;
            double phase = -TWO_PI * (
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
    int nx = static_cast<int>(std::round(a / sp));
    int ny = static_cast<int>(std::round(b / sp));
    int nz = static_cast<int>(std::round(c / sp));

    size_t grid_size = static_cast<size_t>(nx) * ny * nz;
    fftw_complex* Fc_3d = fftw_alloc_complex(grid_size);
    std::fill(reinterpret_cast<double*>(Fc_3d),
              reinterpret_cast<double*>(Fc_3d) + 2 * grid_size, 0.0);

    for (size_t i = 0; i < n_refl; ++i) {
        int hi = ((miller[i].h % nx) + nx) % nx;
        int ki = ((miller[i].k % ny) + ny) % ny;
        int li = ((miller[i].l % nz) + nz) % nz;
        size_t flat = hi * ny * nz + ki * nz + li;
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
        fftw_complex* mask_fft = fftw_alloc_complex(grid_size);
        fftw_complex* mask_in = fftw_alloc_complex(grid_size);
        for (size_t i = 0; i < grid_size; ++i) {
            mask_in[i][0] = sol_mask[i];
            mask_in[i][1] = 0.0;
        }

        fftw_plan mask_plan = fftw_plan_dft_3d(
            nx, ny, nz, mask_in, mask_fft, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(mask_plan);
        fftw_destroy_plan(mask_plan);
        fftw_free(mask_in);

        // Compute S^2 for each FFT grid point and apply correction
        for (int i = 0; i < nx; ++i) {
            double hi = (i <= nx / 2) ? i : i - nx;
            for (int j = 0; j < ny; ++j) {
                double ki = (j <= ny / 2) ? j : j - ny;
                for (int k = 0; k < nz; ++k) {
                    double li = (k <= nz / 2) ? k : k - nz;
                    double s2 = (hi / a) * (hi / a) +
                                (ki / b) * (ki / b) +
                                (li / c) * (li / c);
                    double correction = k_sol * std::exp(-b_sol * s2 / 4.0);
                    size_t flat = i * ny * nz + j * nz + k;
                    Fc_3d[flat][0] += correction * mask_fft[flat][0];
                    Fc_3d[flat][1] += correction * mask_fft[flat][1];
                }
            }
        }

        fftw_free(mask_fft);
    }

    // ----------------------------------------------------------------
    // Step 8: Inverse FFT -> real-space density
    // ----------------------------------------------------------------
    fftw_complex* rho_complex = fftw_alloc_complex(grid_size);
    fftw_plan ifft_plan = fftw_plan_dft_3d(
        nx, ny, nz, Fc_3d, rho_complex, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(ifft_plan);
    fftw_destroy_plan(ifft_plan);
    fftw_free(Fc_3d);

    double V = a * b * c;
    std::vector<double> rho_3d(grid_size);
    double scale_factor = static_cast<double>(nx * ny * nz) / V;
    for (size_t i = 0; i < grid_size; ++i) {
        rho_3d[i] = rho_complex[i][0] * scale_factor /
                     static_cast<double>(grid_size);
    }

    // ----------------------------------------------------------------
    // Step 9: Per-shell amplitude scaling (optional)
    // ----------------------------------------------------------------
    if (n_scale_shells > 1) {
        // FFT the calculated density
        fftw_complex* rho_in = fftw_alloc_complex(grid_size);
        fftw_complex* F_calc = fftw_alloc_complex(grid_size);
        for (size_t i = 0; i < grid_size; ++i) {
            rho_in[i][0] = rho_3d[i];
            rho_in[i][1] = 0.0;
        }
        fftw_plan fwd_plan = fftw_plan_dft_3d(
            nx, ny, nz, rho_in, F_calc, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(fwd_plan);
        fftw_destroy_plan(fwd_plan);
        fftw_free(rho_in);

        // Sample observed density onto UC grid and FFT
        fftw_complex* obs_in = fftw_alloc_complex(grid_size);
        fftw_complex* Fobs_3d = fftw_alloc_complex(grid_size);

        GridParams obs_gp = GetGridParams(obs_grid);
        for (int i = 0; i < nx; ++i) {
            double x = obs_grid.GetXMin() + (static_cast<double>(i) / nx) * a;
            for (int j = 0; j < ny; ++j) {
                double y = obs_grid.GetYMin() +
                           (static_cast<double>(j) / ny) * b;
                for (int k = 0; k < nz; ++k) {
                    double z = obs_grid.GetZMin() +
                               (static_cast<double>(k) / nz) * c;
                    size_t flat = i * ny * nz + j * nz + k;
                    obs_in[flat][0] = InterpolateDensity(
                        obs_grid, x, y, z, 0.0);
                    obs_in[flat][1] = 0.0;
                }
            }
        }

        fftw_plan obs_fwd = fftw_plan_dft_3d(
            nx, ny, nz, obs_in, Fobs_3d, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(obs_fwd);
        fftw_destroy_plan(obs_fwd);
        fftw_free(obs_in);

        // Per-shell scaling
        double s2_max = 1.0 / (resolution * resolution);
        std::vector<double> shell_edges(n_scale_shells + 1);
        for (unsigned int i = 0; i <= n_scale_shells; ++i) {
            shell_edges[i] = s2_max * i / n_scale_shells;
        }

        for (unsigned int shell = 0; shell < n_scale_shells; ++shell) {
            double s2_lo = shell_edges[shell];
            double s2_hi = shell_edges[shell + 1];

            double sum_fobs_fc = 0.0;
            double sum_fc2 = 0.0;

            for (int i = 0; i < nx; ++i) {
                double hi = (i <= nx / 2) ? i : i - nx;
                for (int j = 0; j < ny; ++j) {
                    double ki = (j <= ny / 2) ? j : j - ny;
                    for (int k = 0; k < nz; ++k) {
                        double li = (k <= nz / 2) ? k : k - nz;
                        double s2 = (hi / a) * (hi / a) +
                                    (ki / b) * (ki / b) +
                                    (li / c) * (li / c);

                        bool in_shell;
                        if (shell == 0) {
                            in_shell = (s2 > 0) && (s2 <= s2_hi);
                        } else {
                            in_shell = (s2 > s2_lo) && (s2 <= s2_hi);
                        }

                        if (in_shell) {
                            size_t flat = i * ny * nz + j * nz + k;
                            double fc_amp = std::sqrt(
                                F_calc[flat][0] * F_calc[flat][0] +
                                F_calc[flat][1] * F_calc[flat][1]);
                            double fobs_amp = std::sqrt(
                                Fobs_3d[flat][0] * Fobs_3d[flat][0] +
                                Fobs_3d[flat][1] * Fobs_3d[flat][1]);
                            sum_fobs_fc += fobs_amp * fc_amp;
                            sum_fc2 += fc_amp * fc_amp;
                        }
                    }
                }
            }

            if (sum_fc2 > 0.0) {
                double k_shell = sum_fobs_fc / sum_fc2;
                for (int i = 0; i < nx; ++i) {
                    double hi = (i <= nx / 2) ? i : i - nx;
                    for (int j = 0; j < ny; ++j) {
                        double ki = (j <= ny / 2) ? j : j - ny;
                        for (int k = 0; k < nz; ++k) {
                            double li = (k <= nz / 2) ? k : k - nz;
                            double s2 = (hi / a) * (hi / a) +
                                        (ki / b) * (ki / b) +
                                        (li / c) * (li / c);
                            bool in_shell;
                            if (shell == 0) {
                                in_shell = (s2 > 0) && (s2 <= s2_hi);
                            } else {
                                in_shell = (s2 > s2_lo) && (s2 <= s2_hi);
                            }
                            if (in_shell) {
                                size_t flat = i * ny * nz + j * nz + k;
                                F_calc[flat][0] *= k_shell;
                                F_calc[flat][1] *= k_shell;
                            }
                        }
                    }
                }
            }
        }

        fftw_free(Fobs_3d);

        // Inverse FFT scaled Fc back to real space
        fftw_complex* rho_scaled = fftw_alloc_complex(grid_size);
        fftw_plan scale_ifft = fftw_plan_dft_3d(
            nx, ny, nz, F_calc, rho_scaled, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(scale_ifft);
        fftw_destroy_plan(scale_ifft);
        fftw_free(F_calc);

        for (size_t i = 0; i < grid_size; ++i) {
            rho_3d[i] = rho_scaled[i][0] / static_cast<double>(grid_size);
        }
        fftw_free(rho_scaled);
    }

    fftw_free(rho_complex);

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
        float fx = scale_coords[0], fy = scale_coords[1], fz = scale_coords[2];

        if (obs_grid.IsInGrid(fx, fy, fz)) {
            double obs_val = OESystem::OEFloatGridLinearInterpolate(
                obs_grid, fx, fy, fz, 0.0f);
            double calc_val = OESystem::OEFloatGridLinearInterpolate(
                *calc_grid, fx, fy, fz, 0.0f);
            sum_obs_calc += obs_val * calc_val;
            sum_calc2 += calc_val * calc_val;
        }
    }

    if (sum_calc2 > 0.0) {
        double k_scale = sum_obs_calc / sum_calc2;
        unsigned int grid_sz = calc_grid->GetSize();
        for (unsigned int i = 0; i < grid_sz; ++i) {
            (*calc_grid)[i] = static_cast<float>((*calc_grid)[i] * k_scale);
        }
    }

    return calc_grid;
}

}  // namespace Maptitude
