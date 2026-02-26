"""Benchmark: Scoring metrics — speed + agreement between maptitude and bms-bio.

Benchmarks all 5 scoring metrics (RSCC, RSR, Q-score, EDIAm, Coverage) on
the EM dataset (7CEC) and a representative X-ray dataset (340d). Each package
uses its own fc_density for RSCC/RSR (self-contained pipeline).

Usage:
    python benchmarks/bench_scoring.py
"""

import math

import numpy as np
from openeye import oechem

import maptitude as mpt
from bms_bio.grid import create as bms_create
from bms_bio.grid import score as bms_score

from helpers import (
    ASSET_DIR, RESOLUTIONS, load_mol, load_mrc_grid, load_ccp4_grid,
    wrap_and_pad, bench,
)


def _atom_correlation(by_atom_a: dict, by_atom_b: dict) -> float:
    """Pearson correlation between two per-atom score dicts."""
    common = set(by_atom_a.keys()) & set(by_atom_b.keys())
    vals_a, vals_b = [], []
    for k in sorted(common):
        va, vb = by_atom_a[k], by_atom_b[k]
        if not (math.isnan(va) or math.isnan(vb)):
            vals_a.append(va)
            vals_b.append(vb)
    if len(vals_a) < 3:
        return float("nan")
    return float(np.corrcoef(vals_a, vals_b)[0, 1])


def _run_dataset(label: str, mol, obs_grid, resolution, cell_dims=None,
                 symops_text=None):
    """Benchmark all 5 metrics for one dataset."""
    print(f"\nDataset: {label}")
    print("-" * 100)

    # Pre-compute calc grids for RSCC/RSR (each package uses its own)
    is_xray = cell_dims is not None

    if is_xray:
        a, b, c = cell_dims
        bms_fc = bms_create.fc_density(
            mol, obs_grid, resolution, (a, b, c), k_sol=0.35)

        cell = mpt.UnitCell(a, b, c, 90.0, 90.0, 90.0)
        symops = (mpt.parse_symops(symops_text) if symops_text else None)
        mpt_fc = mpt.fc_density(
            mol, obs_grid, resolution, cell, k_sol=0.35, symops=symops)
    else:
        # EM — use simple pseudo-cell from grid dimensions
        sp = obs_grid.GetSpacing()
        pseudo_cell = (
            obs_grid.GetXDim() * sp,
            obs_grid.GetYDim() * sp,
            obs_grid.GetZDim() * sp,
        )
        bms_fc = bms_create.fc_density(
            mol, obs_grid, resolution, pseudo_cell, k_sol=0.0)
        mpt_fc = None  # EM RSCC/RSR will use bms fc for both

    # Configure maptitude RSCC to use Scaled@1.5 (matches bms-bio default)
    rscc_opts = mpt.RSCCOptions()
    rscc_opts.atom_radius_method = mpt.AtomRadius.Scaled
    rscc_opts.atom_radius_scaling = 1.5

    # Define metric runners
    # Each entry: (name, bms_func, mpt_func)
    metrics = []

    if is_xray:
        metrics.append(("rscc",
            lambda: bms_score.rscc(mol, obs_grid, resolution,
                                   calc_grid=bms_fc),
            lambda: mpt.rscc(mol, obs_grid, resolution,
                             calc_grid=mpt_fc, options=rscc_opts)))
        metrics.append(("rsr",
            lambda: bms_score.rsr(mol, obs_grid, resolution,
                                  calc_grid=bms_fc),
            lambda: mpt.rsr(mol, obs_grid, resolution,
                            calc_grid=mpt_fc)))
    else:
        metrics.append(("rscc",
            lambda: bms_score.rscc(mol, obs_grid, resolution,
                                   calc_grid=bms_fc),
            lambda: mpt.rscc(mol, obs_grid, resolution,
                             calc_grid=bms_fc, options=rscc_opts)))
        metrics.append(("rsr",
            lambda: bms_score.rsr(mol, obs_grid, resolution,
                                  calc_grid=bms_fc),
            lambda: mpt.rsr(mol, obs_grid, resolution,
                            calc_grid=bms_fc)))

    metrics.extend([
        ("qscore",
            lambda: bms_score.qscore(mol, obs_grid, resolution),
            lambda: mpt.qscore(mol, obs_grid, resolution)),
        ("ediam",
            lambda: bms_score.ediam(mol, obs_grid, resolution),
            lambda: mpt.ediam(mol, obs_grid, resolution)),
        ("coverage",
            lambda: bms_score.coverage(mol, obs_grid, 1.0),
            lambda: mpt.coverage(mol, obs_grid, 1.0)),
    ])

    # Run benchmarks
    print(f"{'Metric':10s} {'bms-bio (ms)':>13s} {'maptitude (ms)':>15s} "
          f"{'Speedup':>8s} {'bms-bio':>9s} {'maptitude':>10s} "
          f"{'Delta':>7s} {'Atom r':>7s}")
    print("-" * 85)

    for name, bms_fn, mpt_fn in metrics:
        t_bms, r_bms = bench(bms_fn, n=5, warmup=1)
        t_mpt, r_mpt = bench(mpt_fn, n=5, warmup=1)

        speedup = t_bms / t_mpt if t_mpt > 0 else float("inf")
        delta = abs(r_bms.overall - r_mpt.overall)

        bms_atoms = {k: v for k, v in r_bms.by_atom.items()}
        mpt_atoms = {k: v for k, v in r_mpt.by_atom.items()}
        corr = _atom_correlation(bms_atoms, mpt_atoms)

        print(f"{name:10s} {t_bms*1000:13.2f} {t_mpt*1000:15.2f} "
              f"{speedup:7.1f}x {r_bms.overall:9.4f} {r_mpt.overall:10.4f} "
              f"{delta:7.4f} {corr:7.3f}")

    print()


def run():
    print("=" * 72)
    print("Scoring Metrics: Speed + Agreement (maptitude vs bms-bio)")
    print("=" * 72)

    # EM dataset (7CEC)
    mol_em = load_mol(ASSET_DIR / "390_7cec_A100.cif")
    grid_em = load_mrc_grid(ASSET_DIR / "390_emd_30342_A_z4.mrc")
    n_em = sum(1 for _ in mol_em.GetAtoms(oechem.OEIsHeavy()))
    print(f"\nEM: 7CEC, {n_em} atoms, res={RESOLUTIONS['7cec']} A")

    _run_dataset(
        f"EM (7CEC, {n_em} atoms, {RESOLUTIONS['7cec']} A)",
        mol_em, grid_em, RESOLUTIONS["7cec"])

    # X-ray dataset (340d)
    mol_xr = load_mol(ASSET_DIR / "340d.cif")
    grid_xr, cell_xr, sym_xr = load_ccp4_grid(ASSET_DIR / "340d_2fofc.ccp4")
    grid_xr = wrap_and_pad(grid_xr, mol_xr, cell_xr)
    n_xr = sum(1 for _ in mol_xr.GetAtoms(oechem.OEIsHeavy()))
    print(f"X-ray: 340d, {n_xr} atoms, res={RESOLUTIONS['340d']} A")

    _run_dataset(
        f"X-ray (340d, {n_xr} atoms, {RESOLUTIONS['340d']} A)",
        mol_xr, grid_xr, RESOLUTIONS["340d"],
        cell_dims=cell_xr, symops_text=sym_xr)


if __name__ == "__main__":
    run()
