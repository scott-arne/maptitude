"""Grids returned to Python must be independently owned copies.

The wrapper used to delete OpenEye's allocation and swap in maptitude's
pointer, so the two runtimes freed each other's memory. These tests verify
that the returned grid is a correct, independently owned copy with the right
values and geometry.
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
    """The result must still be readable and correct after its inputs are gone."""
    result = maptitude.combine_maps(_make_grid(1.0), _make_grid(2.0), maptitude.MapOp.ADD)
    # Belt and braces: refcounting has already freed the temporaries above.
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


def test_repeated_allocation_completes_without_error() -> None:
    """The return path must not throw or crash under repeated use."""
    for _ in range(300):
        result = maptitude.combine_maps(_make_grid(1.0), _make_grid(2.0), maptitude.MapOp.ADD)
        assert result.GetValue(0) == pytest.approx(3.0)
        del result
    gc.collect()


def test_combine_maps_rejects_zero_spacing() -> None:
    """Grids with zero spacing produce degenerate geometry that operator= silently loses.

    The typemap helper must detect and reject the silent geometry corruption rather
    than returning a plausible-looking 1x1x1 grid holding 0.0.
    """
    grid_a = oegrid.OEScalarGrid(oechem.OEDoubleArray([-4.0, -4.0, -4.0, 4.0, 4.0, 4.0]), 0.0)
    for i in range(grid_a.GetSize()):
        grid_a.SetValue(i, 1.0)
    grid_b = oegrid.OEScalarGrid(oechem.OEDoubleArray([-4.0, -4.0, -4.0, 4.0, 4.0, 4.0]), 0.0)
    for i in range(grid_b.GetSize()):
        grid_b.SetValue(i, 2.0)
    with pytest.raises(RuntimeError):
        maptitude.combine_maps(grid_a, grid_b, maptitude.MapOp.ADD)


def test_combine_maps_rejects_infinite_spacing() -> None:
    """Grids with infinite spacing trigger a setter failure that operator= ignores.

    OpenEye prints 'Warning: SetSpacing unable to handle NaN: inf' and the
    destination stays at default geometry. The typemap helper must reject that.
    """
    grid_a = oegrid.OEScalarGrid(oechem.OEDoubleArray([-4.0, -4.0, -4.0, 4.0, 4.0, 4.0]), float("inf"))
    for i in range(grid_a.GetSize()):
        grid_a.SetValue(i, 1.0)
    grid_b = oegrid.OEScalarGrid(oechem.OEDoubleArray([-4.0, -4.0, -4.0, 4.0, 4.0, 4.0]), float("inf"))
    for i in range(grid_b.GetSize()):
        grid_b.SetValue(i, 2.0)
    with pytest.raises(RuntimeError):
        maptitude.combine_maps(grid_a, grid_b, maptitude.MapOp.ADD)
