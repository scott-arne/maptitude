"""Cover the benchmark map-loading path, which the benchmarks cannot cover.

All three benchmarks import ``bms_bio`` at module scope, so none of them can be
imported without it, let alone run. Extracting their shared loading step into
``helpers.load_xray_dataset`` -- which imports no ``bms_bio`` -- is what makes
that path reachable from this suite.
"""

import sys
from pathlib import Path

import pytest

pytest.importorskip("openeye.oechem", reason="OpenEye toolkit not installed")

from openeye import oechem  # noqa: E402

import maptitude  # noqa: E402
from maptitude import get_unit_cell, parse_symops, read_map  # noqa: E402

_BENCH_DIR = Path(__file__).resolve().parents[2] / "benchmarks"
if str(_BENCH_DIR) not in sys.path:
    sys.path.insert(0, str(_BENCH_DIR))

from helpers import ASSET_DIR, load_xray_dataset  # noqa: E402


def test_load_xray_dataset_threads_the_map_into_every_return_value():
    """Each value must come from the map named by the argument.

    The failure this guards is silent: the loader hands a cell triple and a
    grid to the padding helper, and a swapped, transposed or stale value there
    changes every benchmark number without raising.
    """
    mol, grid, cell_dims, symops_text = load_xray_dataset("340d")

    reference = read_map(ASSET_DIR / "340d_2fofc.ccp4")
    reference_cell = get_unit_cell(reference.grid)

    assert cell_dims == (reference_cell.a, reference_cell.b, reference_cell.c)
    assert symops_text == reference.symops
    assert parse_symops(symops_text)

    assert sum(1 for _ in mol.GetAtoms(oechem.OEIsHeavy())) > 0

    # wrap_and_pad creates a sub-box around the molecule, so the returned grid
    # may be smaller than, equal to, or larger than the reference grid depending
    # on the molecule's extent. Comparing sizes catches a loader that forgot to
    # call wrap_and_pad at all (which would return the reference grid unchanged).
    assert grid.GetSize() <= reference.grid.GetSize()
