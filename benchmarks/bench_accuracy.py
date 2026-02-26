"""Benchmark: Accuracy vs published reference values.

Part A: RSCC/RSR vs wwPDB validation pipeline (Smart et al., 2018) for 3 X-ray
structures, comparing maptitude (Scaled@1.5 and Binned strategies) and bms-bio
against published per-ligand values.

Part B: Q-score vs mapq reference (Pintilie et al., 2020) for the 7CEC EM
dataset, comparing per-residue Q-scores.

Usage:
    python benchmarks/bench_accuracy.py
"""

import math

import numpy as np
from openeye import oechem

import maptitude as mpt
from bms_bio.grid import create as bms_create
from bms_bio.grid import score as bms_score

from helpers import (
    ASSET_DIR, RESOLUTIONS, PUBLISHED_RSCC_RSR, MAPQ_QSCORES,
    load_mol, load_mrc_grid, load_ccp4_grid, wrap_and_pad, ligand_mask,
)


# ---------------------------------------------------------------------------
# Part A: RSCC/RSR vs published wwPDB values
# ---------------------------------------------------------------------------

def _part_a():
    print("=" * 72)
    print("Part A: RSCC/RSR vs Published wwPDB Values (Smart et al., 2018)")
    print("=" * 72)

    # Load all 3 X-ray structures
    data = {}
    for pdb in ["1d26", "3q9g", "340d"]:
        mol = load_mol(ASSET_DIR / f"{pdb}.cif")
        grid, cell_dims, symops_text = load_ccp4_grid(
            ASSET_DIR / f"{pdb}_2fofc.ccp4")
        grid = wrap_and_pad(grid, mol, cell_dims)
        data[pdb] = {
            "mol": mol, "grid": grid, "cell": cell_dims,
            "symops_text": symops_text,
        }

    # Compute Fc density and per-ligand RSCC/RSR for each implementation
    bms_results = {}  # pdb -> {(resname, resnum): {"rscc": ..., "rsr": ...}}
    mpt_scaled = {}   # Scaled@1.5
    mpt_binned = {}   # Binned

    for pdb in data:
        d = data[pdb]
        a, b, c = d["cell"]
        res = RESOLUTIONS[pdb]

        # bms-bio fc_density + scoring
        bms_fc = bms_create.fc_density(
            d["mol"], d["grid"], res, (a, b, c), k_sol=0.35)

        # maptitude fc_density
        cell = mpt.UnitCell(a, b, c, 90.0, 90.0, 90.0)
        symops = (mpt.parse_symops(d["symops_text"])
                  if d["symops_text"] else None)
        mpt_fc = mpt.fc_density(
            d["mol"], d["grid"], res, cell, k_sol=0.35, symops=symops)

        # RSCC options
        opts_scaled = mpt.RSCCOptions()
        opts_scaled.atom_radius_method = mpt.AtomRadius.Scaled
        opts_scaled.atom_radius_scaling = 1.5

        opts_binned = mpt.RSCCOptions()
        opts_binned.atom_radius_method = mpt.AtomRadius.Binned

        bms_results[pdb] = {}
        mpt_scaled[pdb] = {}
        mpt_binned[pdb] = {}

        for entry in PUBLISHED_RSCC_RSR[pdb]:
            mask = ligand_mask(entry["resname"], entry["chain"],
                               entry["resnum"])
            key = (entry["resname"], entry["resnum"])

            # bms-bio
            r_rscc = bms_score.rscc(
                d["mol"], d["grid"], res, mask=mask, calc_grid=bms_fc)
            r_rsr = bms_score.rsr(
                d["mol"], d["grid"], res, mask=mask, calc_grid=bms_fc)
            bms_results[pdb][key] = {
                "rscc": r_rscc.overall, "rsr": r_rsr.overall}

            # maptitude Scaled@1.5
            r_rscc_s = mpt.RSCC(d["mol"], d["grid"], res,
                                mask, mpt_fc, opts_scaled)
            r_rsr_s = mpt.rsr(d["mol"], d["grid"], res,
                              mask=mask, calc_grid=mpt_fc)
            mpt_scaled[pdb][key] = {
                "rscc": r_rscc_s.overall, "rsr": r_rsr_s.overall}

            # maptitude Binned
            r_rscc_b = mpt.RSCC(d["mol"], d["grid"], res,
                                mask, mpt_fc, opts_binned)
            mpt_binned[pdb][key] = {"rscc": r_rscc_b.overall}

    # --- RSCC table ---
    print("\nPer-Ligand RSCC vs Published")
    print(f"{'PDB':5s} {'Ligand':10s} {'Published':>10s} {'bms-bio':>9s} "
          f"{'Scaled@1.5':>11s} {'Binned':>8s} "
          f"{'d(Scaled)':>10s} {'d(Binned)':>10s}")
    print("-" * 80)

    for pdb in ["1d26", "3q9g", "340d"]:
        for entry in PUBLISHED_RSCC_RSR[pdb]:
            key = (entry["resname"], entry["resnum"])
            pub = entry["rscc"]
            bms = bms_results[pdb][key]["rscc"]
            s15 = mpt_scaled[pdb][key]["rscc"]
            bn = mpt_binned[pdb][key]["rscc"]
            label = f"{entry['resname']}:{entry['resnum']}"

            print(f"{pdb:5s} {label:10s} {pub:10.3f} {bms:9.3f} "
                  f"{s15:11.3f} {bn:8.3f} "
                  f"{s15 - pub:+10.3f} {bn - pub:+10.3f}")

    # --- RSR table ---
    print("\nPer-Ligand RSR vs Published")
    print(f"{'PDB':5s} {'Ligand':10s} {'Published':>10s} {'bms-bio':>9s} "
          f"{'maptitude':>10s} {'d(maptitude)':>13s}")
    print("-" * 65)

    for pdb in ["1d26", "3q9g", "340d"]:
        for entry in PUBLISHED_RSCC_RSR[pdb]:
            key = (entry["resname"], entry["resnum"])
            pub = entry["rsr"]
            bms = bms_results[pdb][key]["rsr"]
            mpt_val = mpt_scaled[pdb][key]["rsr"]
            label = f"{entry['resname']}:{entry['resnum']}"

            print(f"{pdb:5s} {label:10s} {pub:10.3f} {bms:9.3f} "
                  f"{mpt_val:10.3f} {mpt_val - pub:+13.3f}")

    print()


# ---------------------------------------------------------------------------
# Part B: Q-score vs mapq reference
# ---------------------------------------------------------------------------

def _part_b():
    print("=" * 72)
    print("Part B: Q-Score vs mapq Reference (7CEC / EMD-30342)")
    print("=" * 72)

    mol = load_mol(ASSET_DIR / "390_7cec_A100.cif")
    grid = load_mrc_grid(ASSET_DIR / "390_emd_30342_A_z4.mrc")
    resolution = RESOLUTIONS["7cec"]

    result = mpt.qscore(mol, grid, resolution)

    # Build residue-number -> score mapping from maptitude result
    mpt_by_res = {}
    for res, score in result.by_residue.items():
        mpt_by_res[res.number] = score

    print(f"\nPer-Residue Q-Score vs mapq Reference (7CEC)")
    print(f"{'Res#':>5s} {'Name':5s} {'mapq':>8s} {'maptitude':>10s} "
          f"{'Delta':>8s}")
    print("-" * 40)

    mapq_vals = []
    mpt_vals = []

    for resnum in sorted(MAPQ_QSCORES.keys()):
        name, mapq_q = MAPQ_QSCORES[resnum]
        mpt_q = mpt_by_res.get(resnum, float("nan"))
        delta = mpt_q - mapq_q if not math.isnan(mpt_q) else float("nan")

        if not math.isnan(mpt_q):
            mapq_vals.append(mapq_q)
            mpt_vals.append(mpt_q)

        delta_str = f"{delta:+8.4f}" if not math.isnan(delta) else "     NaN"
        mpt_str = f"{mpt_q:10.4f}" if not math.isnan(mpt_q) else "       NaN"
        print(f"{resnum:5d} {name:5s} {mapq_q:8.4f} {mpt_str} {delta_str}")

    # Overall
    mapq_overall = np.mean([v for _, (_, v) in
                            sorted(MAPQ_QSCORES.items())])
    print(f"{'':5s} {'Overall':5s} {mapq_overall:8.4f} "
          f"{result.overall:10.4f} {result.overall - mapq_overall:+8.4f}")

    # Correlation
    if len(mapq_vals) >= 3:
        corr = float(np.corrcoef(mapq_vals, mpt_vals)[0, 1])
        print(f"\nPearson r: {corr:.3f}")
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run():
    _part_a()
    _part_b()


if __name__ == "__main__":
    run()
