"""Benchmark: Fc density generation — maptitude vs bms-bio.

Compares timing and grid-level Pearson correlation between the two
implementations on all 3 X-ray structures (1d26, 3q9g, 340d).

Usage:
    python benchmarks/bench_fc_density.py
"""

import numpy as np
from openeye import oechem

import maptitude as mpt
from bms_bio.grid import create as bms_create

from helpers import (
    ASSET_DIR, RESOLUTIONS, load_mol, load_ccp4_grid, wrap_and_pad, bench,
)


def _grid_to_array(grid):
    """Extract grid values as a 1-D numpy array."""
    return np.array([grid.GetValue(i) for i in range(grid.GetSize())])


def _pearson(a, b):
    """Pearson correlation between two arrays."""
    if len(a) < 2:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def run():
    print("=" * 72)
    print("Fc Density Generation: maptitude vs bms-bio")
    print("=" * 72)

    # Load all X-ray datasets
    datasets = {}
    for pdb in ["1d26", "3q9g", "340d"]:
        mol = load_mol(ASSET_DIR / f"{pdb}.cif")
        grid, cell_dims, symops_text = load_ccp4_grid(
            ASSET_DIR / f"{pdb}_2fofc.ccp4")
        grid = wrap_and_pad(grid, mol, cell_dims)
        n_heavy = sum(1 for _ in mol.GetAtoms(oechem.OEIsHeavy()))
        datasets[pdb] = {
            "mol": mol, "grid": grid, "cell": cell_dims,
            "symops_text": symops_text, "n_heavy": n_heavy,
        }
        print(f"  {pdb}: {n_heavy} heavy atoms, "
              f"grid {grid.GetXDim()}x{grid.GetYDim()}x{grid.GetZDim()}, "
              f"res={RESOLUTIONS[pdb]} A")

    print()

    # Benchmark
    rows = []
    for pdb in ["1d26", "3q9g", "340d"]:
        d = datasets[pdb]
        a, b, c = d["cell"]
        res = RESOLUTIONS[pdb]

        # bms-bio
        t_bms, fc_bms = bench(
            lambda: bms_create.fc_density(
                d["mol"], d["grid"], res, (a, b, c), k_sol=0.35),
            n=3, warmup=1)

        # maptitude
        cell = mpt.UnitCell(a, b, c, 90.0, 90.0, 90.0)
        symops = (mpt.parse_symops(d["symops_text"])
                  if d["symops_text"] else None)
        t_mpt, fc_mpt = bench(
            lambda: mpt.fc_density(
                d["mol"], d["grid"], res, cell, k_sol=0.35, symops=symops),
            n=3, warmup=1)

        # Grid correlation
        arr_bms = _grid_to_array(fc_bms)
        arr_mpt = _grid_to_array(fc_mpt)
        r = _pearson(arr_bms, arr_mpt)

        speedup = t_bms / t_mpt if t_mpt > 0 else float("inf")
        rows.append({
            "pdb": pdb, "atoms": d["n_heavy"],
            "bms_ms": t_bms * 1000, "mpt_ms": t_mpt * 1000,
            "speedup": speedup, "grid_r": r,
        })

    # Report
    print(f"{'PDB':5s} {'Atoms':>6s} {'bms-bio (ms)':>13s} "
          f"{'maptitude (ms)':>15s} {'Speedup':>8s} {'Grid r':>8s}")
    print("-" * 60)
    for r in rows:
        print(f"{r['pdb']:5s} {r['atoms']:6d} {r['bms_ms']:13.1f} "
              f"{r['mpt_ms']:15.1f} {r['speedup']:7.1f}x {r['grid_r']:8.4f}")
    print()


if __name__ == "__main__":
    run()
