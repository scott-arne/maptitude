"""Regenerate the committed CCP4 test map.

The output is committed so tests do not depend on a Python toolchain. Rerun only
if the fixture needs to change::

    .venv/bin/python tests/data/make_test_map.py

Note: rerunning this script does not reproduce the committed file. That file was
written by the earlier scalar-carrier path, which departs from the analytic
Gaussian by up to 8.4e-06; the skew path used here departs by up to 2.0e-08, so
the density region changes by up to 8.4e-06.
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
    dim = int(2 * half_width / spacing) + 1

    grid = oegrid.OESkewGrid()
    assert grid.SetDim(dim, dim, dim)
    # The cell is deliberately twice the sampled extent, at twice the division
    # count, and the space group is set explicitly. The CCP4 writer needs both,
    # the grid needs neither: with no space group the writer emits an origin one
    # node too high, and with one it emits the grid's own cell, which the reader
    # then re-expands by a node per axis on any map whose sampled count equals
    # its cell division count. Doubling keeps the two apart. Spacing, midpoint
    # and node coordinates are identical under either cell, so use the doubled
    # form only where a grid is about to be written.
    cell_edge = 2 * dim * spacing
    assert grid.SetUnitCell(cell_edge, cell_edge, cell_edge, 90.0, 90.0, 90.0,
                            2 * dim, 2 * dim, 2 * dim)
    assert grid.SetMid(0.0, 0.0, 0.0)
    assert grid.SetSpaceGroup(1)

    sigma = 1.2
    two_sigma_sq = 2.0 * sigma * sigma
    values = oechem.OEFloatArray(grid.GetSize())
    for i in range(grid.GetSize()):
        x, y, z = grid.ElementToSpatialCoord(i)
        values[i] = math.exp(-(x * x + y * y + z * z) / two_sigma_sq)
    assert grid.SetValues(values, grid.GetSize())

    if not oegrid.OEWriteGrid(str(path), grid):
        raise RuntimeError(f"OEWriteGrid failed for {path}")


if __name__ == "__main__":
    build_map(Path(__file__).parent / "test_map.ccp4")
