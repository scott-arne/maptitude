"""Cover the benchmark map-loading path, which the benchmarks cannot cover.

All three benchmarks import ``bms_bio`` at module scope, so none of them can be
imported without it, let alone run. Extracting their shared loading step into
``helpers.load_xray_dataset`` -- which imports no ``bms_bio`` -- is what makes
that path reachable from this suite.
"""

import sys
from pathlib import Path

import pytest
from maptitude import get_unit_cell, parse_symops, read_map

oechem = pytest.importorskip("openeye.oechem", reason="OpenEye toolkit not installed")

_BENCH_DIR = Path(__file__).resolve().parents[2] / "benchmarks"
if str(_BENCH_DIR) not in sys.path:
    sys.path.insert(0, str(_BENCH_DIR))

from helpers import ASSET_DIR, load_xray_dataset


def test_load_xray_dataset_threads_the_map_into_every_return_value():
    """Each value must come from the map named by the argument.

    The failure this guards is silent: the loader hands a cell triple and a
    grid to the padding helper, and a stale or misordered value there shifts
    the scores computed from that dataset without raising. 340d's a and b
    edges are equal, so an a/b transposition is invisible here -- what this
    pins is that each value tracks the map, not that every permutation of the
    triple is distinguishable.
    """
    mol, grid, cell_dims, symops_text = load_xray_dataset("340d")

    reference = read_map(ASSET_DIR / "340d_2fofc.ccp4")
    reference_cell = get_unit_cell(reference.grid)

    assert cell_dims == (reference_cell.a, reference_cell.b, reference_cell.c)
    assert symops_text == reference.symops
    assert parse_symops(symops_text)

    assert sum(1 for _ in mol.GetAtoms(oechem.OEIsHeavy())) > 0

    # Skipping the padding step would hand back the map's own grid, so the
    # inequality has to be strict to catch it -- `<=` admits that exact case.
    # For 340d the box cut around the molecule is 36x53x35 against the map's
    # 61x61x37; this pins that the step ran, not the size it produced.
    assert grid.GetSize() < reference.grid.GetSize()
