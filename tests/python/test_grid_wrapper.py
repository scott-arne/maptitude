"""Grids returned to Python must be independently owned copies.

The wrapper used to delete OpenEye's allocation and swap in maptitude's
pointer, so the two runtimes freed each other's memory. These tests exercise
the return path hard enough that a heap problem shows up as a crash.
"""

from __future__ import annotations

import gc

import pytest

import maptitude

oechem = pytest.importorskip("openeye.oechem")
oegrid = pytest.importorskip("openeye.oegrid")


def _make_grid(value: float = 1.0):
    grid = oegrid.OEScalarGrid(oechem.OEDoubleArray([-4.0, -4.0, -4.0, 4.0, 4.0, 4.0]), 1.0)
    for i in range(grid.GetSize()):
        grid.SetValue(i, value)
    return grid


def test_combined_grid_has_the_expected_values() -> None:
    a = _make_grid(1.0)
    b = _make_grid(2.0)
    result = maptitude.combine_maps(a, b, maptitude.MapOp.ADD)
    assert result.GetSize() == a.GetSize()
    for i in range(result.GetSize()):
        assert result.GetValue(i) == pytest.approx(3.0)


def test_returned_grid_survives_its_inputs() -> None:
    """The result must not alias memory owned by the arguments."""
    result = maptitude.combine_maps(_make_grid(1.0), _make_grid(2.0), maptitude.MapOp.ADD)
    gc.collect()
    assert result.GetValue(0) == pytest.approx(3.0)


def test_returned_grid_geometry_matches_the_input() -> None:
    a = _make_grid()
    result = maptitude.combine_maps(a, _make_grid(), maptitude.MapOp.ADD)
    assert result.GetXDim() == a.GetXDim()
    assert result.GetYDim() == a.GetYDim()
    assert result.GetZDim() == a.GetZDim()
    assert result.GetSpacing() == pytest.approx(a.GetSpacing())
    assert result.GetXMid() == pytest.approx(a.GetXMid())
    assert result.GetYMid() == pytest.approx(a.GetYMid())
    assert result.GetZMid() == pytest.approx(a.GetZMid())


def test_repeated_allocation_does_not_corrupt_the_heap() -> None:
    """A single mismatched free may not fault; a few hundred usually will."""
    for _ in range(300):
        result = maptitude.combine_maps(_make_grid(1.0), _make_grid(2.0), maptitude.MapOp.ADD)
        assert result.GetValue(0) == pytest.approx(3.0)
        del result
    gc.collect()
