"""Regenerate the committed CCP4 test map.

The output is committed so tests do not depend on a Python toolchain. Rerun only
if the fixture needs to change::

    .venv/bin/python tests/data/make_test_map.py

Note: regeneration is not byte-identical. The density region is stable, but 2
bytes in the CCP4 EXTRA field change each run due to pointer-derived metadata.
"""

from __future__ import annotations

import math
from pathlib import Path

from openeye import oechem, oegrid


def build_map(path: Path) -> None:
    """Write a small Gaussian density map in CCP4 format.

    :param path: Destination file. Overwritten if it exists.
    """
    half_width = 5.0
    spacing = 0.5
    minmax = oechem.OEDoubleArray([-half_width] * 3 + [half_width] * 3)
    grid = oegrid.OEScalarGrid(minmax, spacing)

    sigma = 1.2
    two_sigma_sq = 2.0 * sigma * sigma
    for i in range(grid.GetSize()):
        x, y, z = grid.ElementToSpatialCoord(i)
        grid.SetValue(i, math.exp(-(x * x + y * y + z * z) / two_sigma_sq))

    if not oegrid.OEWriteGrid(str(path), grid):
        raise RuntimeError(f"OEWriteGrid failed for {path}")


if __name__ == "__main__":
    build_map(Path(__file__).parent / "test_map.ccp4")
