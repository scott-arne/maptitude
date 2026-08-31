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


def _fill(grid, value: float):
    """Set every node of *grid* to *value*.

    OESkewGrid has no per-element setter, so the whole array goes across in one
    SetValues call. GetValues hands back a copy, so writing to that would be
    discarded.
    """
    values = oechem.OEFloatArray(grid.GetSize())
    for i in range(grid.GetSize()):
        values[i] = value
    assert grid.SetValues(values, grid.GetSize())
    return grid


def _make_grid(value: float = 1.0):
    # The 9x9x9 dims and the (0, 0, 0) midpoint are what the extents-box
    # constructor derived from [-4, -4, -4, 4, 4, 4] at a 1.0 A interval.
    grid = oegrid.OESkewGrid()
    assert grid.SetDim(9, 9, 9)
    assert grid.SetUnitCell(9.0, 9.0, 9.0, 90.0, 90.0, 90.0, 9, 9, 9)
    assert grid.SetMid(0.0, 0.0, 0.0)
    return _fill(grid, value)


def _make_degenerate_grid(cell_edge: float, value: float):
    """Build a 9x9x9 grid whose node coordinates are all NaN.

    A zero or infinite cell edge leaves the sampling undefined, so every node
    coordinate comes back NaN and no geometry can be derived from the grid.
    """
    grid = oegrid.OESkewGrid()
    assert grid.SetDim(9, 9, 9)
    assert grid.SetUnitCell(cell_edge, cell_edge, cell_edge, 90.0, 90.0, 90.0, 9, 9, 9)
    return _fill(grid, value)


def test_combined_grid_has_the_expected_values() -> None:
    a = _make_grid(1.0)
    b = _make_grid(2.0)
    result = maptitude.combine_maps(a, b, maptitude.MapOp.ADD)
    assert result.GetSize() == a.GetSize()
    values = result.GetValues()
    for i in range(result.GetSize()):
        assert values[i] == pytest.approx(3.0)


def test_returned_grid_survives_its_inputs() -> None:
    """The result must still be readable and correct after its inputs are gone."""
    result = maptitude.combine_maps(_make_grid(1.0), _make_grid(2.0), maptitude.MapOp.ADD)
    # Belt and braces: refcounting has already freed the temporaries above.
    gc.collect()
    assert result.GetValues()[0] == pytest.approx(3.0)


def test_returned_grid_geometry_matches_the_input() -> None:
    a = _make_grid()
    result = maptitude.combine_maps(a, _make_grid(), maptitude.MapOp.ADD)
    assert result.GetXDim() == a.GetXDim()
    assert result.GetYDim() == a.GetYDim()
    assert result.GetZDim() == a.GetZDim()
    assert maptitude.same_grid_geometry(result, a)
    assert result.GetXMid() == pytest.approx(a.GetXMid())
    assert result.GetYMid() == pytest.approx(a.GetYMid())
    assert result.GetZMid() == pytest.approx(a.GetZMid())


def test_returned_grid_is_an_independent_copy() -> None:
    """Writing to the returned grid must not reach back into either input.

    The out typemap copy-assigns into a Python-owned grid rather than swapping
    pointers, so the result and the inputs must not share storage. Every other
    test in this file only reads the result, so an aliased buffer would satisfy
    all of them; it shows up here, as a wrong value in an input rather than as a
    crash.

    This file also used to reach the typemap's copy validation on its failing
    branch, through the two degenerate-grid tests below. Those now throw inside
    combine_maps, before it returns, and an out typemap only runs on a return.
    """
    a = _make_grid(1.0)
    b = _make_grid(2.0)
    result = maptitude.combine_maps(a, b, maptitude.MapOp.ADD)

    assert result.GetSize() == a.GetSize()
    assert maptitude.same_grid_geometry(result, a)
    assert result.GetValues()[0] == pytest.approx(3.0)

    _fill(result, 99.0)

    assert result.GetValues()[0] == pytest.approx(99.0)
    assert a.GetValues()[0] == pytest.approx(1.0)
    assert b.GetValues()[0] == pytest.approx(2.0)


def test_repeated_allocation_completes_without_error() -> None:
    """The return path must not throw or crash under repeated use."""
    for _ in range(300):
        result = maptitude.combine_maps(_make_grid(1.0), _make_grid(2.0), maptitude.MapOp.ADD)
        assert result.GetValues()[0] == pytest.approx(3.0)
        del result
    gc.collect()


def test_combine_maps_rejects_zero_cell_edge() -> None:
    """A zero cell edge leaves every node coordinate NaN.

    combine_maps gates on same_grid_geometry, which routes through
    get_grid_params and refuses a grid whose sampling cannot be measured, so
    the call is rejected before it computes anything.
    """
    grid_a = _make_degenerate_grid(0.0, 1.0)
    grid_b = _make_degenerate_grid(0.0, 2.0)
    with pytest.raises(maptitude.GridError):
        maptitude.combine_maps(grid_a, grid_b, maptitude.MapOp.ADD)


def test_combine_maps_rejects_infinite_cell_edge() -> None:
    """An infinite cell edge is refused on the same path as a zero one.

    OpenEye warns on stderr that SetSpacing cannot handle the value and leaves
    the node coordinates NaN, which get_grid_params then rejects.
    """
    grid_a = _make_degenerate_grid(float("inf"), 1.0)
    grid_b = _make_degenerate_grid(float("inf"), 2.0)
    with pytest.raises(maptitude.GridError):
        maptitude.combine_maps(grid_a, grid_b, maptitude.MapOp.ADD)


def _make_anisotropic_grid(value: float = 1.0):
    """Build a grid with all three axes distinct in both dim and spacing.

    Uses dims (4, 5, 6) and spacing (0.5, 1.0, 2.0), giving edges (2.0, 5.0, 12.0).
    An x/y spacing or dim swap is immediately visible through get_grid_params.
    """
    grid = oegrid.OESkewGrid()
    assert grid.SetDim(4, 5, 6)
    assert grid.SetUnitCell(2.0, 5.0, 12.0, 90.0, 90.0, 90.0, 4, 5, 6)
    assert grid.SetMid(0.0, 0.0, 0.0)
    gp = maptitude.get_grid_params(grid)
    assert gp.x_dim == 4
    assert gp.y_dim == 5
    assert gp.z_dim == 6
    assert gp.x_spacing == pytest.approx(0.5)
    assert gp.y_spacing == pytest.approx(1.0)
    assert gp.z_spacing == pytest.approx(2.0)
    return _fill(grid, value)


def test_diff_to_calc_preserves_anisotropic_geometry() -> None:
    """diff_to_calc must preserve per-axis dims and spacings in order.

    This test would fail if x and y dims or spacings were transposed by the
    out typemap, because the fixture has all three axes distinct.
    """
    obs = _make_anisotropic_grid(2.0)
    diff = _make_anisotropic_grid(1.0)
    calc = maptitude.diff_to_calc(obs, diff)
    assert calc is not None
    gp = maptitude.get_grid_params(calc)
    assert gp.x_dim == 4
    assert gp.y_dim == 5
    assert gp.z_dim == 6
    assert gp.x_spacing == pytest.approx(0.5)
    assert gp.y_spacing == pytest.approx(1.0)
    assert gp.z_spacing == pytest.approx(2.0)
    assert calc.GetValues()[0] == pytest.approx(0.0)


def test_scale_map_preserves_anisotropic_geometry() -> None:
    """scale_map must preserve per-axis dims and spacings after in-place mutation.

    This test would fail if x and y dims or spacings were transposed by the
    non-const in typemap, because the fixture has all three axes distinct.
    """
    grid = _make_anisotropic_grid(2.0)
    maptitude.scale_map(grid, 3.0)
    gp = maptitude.get_grid_params(grid)
    assert gp.x_dim == 4
    assert gp.y_dim == 5
    assert gp.z_dim == 6
    assert gp.x_spacing == pytest.approx(0.5)
    assert gp.y_spacing == pytest.approx(1.0)
    assert gp.z_spacing == pytest.approx(2.0)
    assert grid.GetValues()[0] == pytest.approx(6.0)


def test_wrap_and_pad_grid_padding_branch_preserves_anisotropic_geometry() -> None:
    """wrap_and_pad_grid must preserve per-axis dims and spacings when padding.

    Places an atom far outside the grid extent to force the padding branch.
    The test would fail if x and y dims or spacings were transposed by the
    out typemap, because the fixture has all three axes distinct.

    The discriminating assertion is 'padded is not obs', which proves the
    padding branch was taken. Without it, the test would pass on the nullptr
    branch (which returns the original grid) and prove nothing about the
    returned-grid path.
    """
    grid = _make_anisotropic_grid(1.0)
    mol = oechem.OEGraphMol()
    atom = mol.NewAtom(6)
    mol.SetCoords(atom, oechem.OEFloatArray([50.0, 50.0, 50.0]))
    padded = maptitude.wrap_and_pad_grid(grid, mol, 2.0, 5.0, 12.0, padding=3.0)
    assert padded is not grid
    gp = maptitude.get_grid_params(padded)
    assert gp.x_spacing == pytest.approx(0.5)
    assert gp.y_spacing == pytest.approx(1.0)
    assert gp.z_spacing == pytest.approx(2.0)
