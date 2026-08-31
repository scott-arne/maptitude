"""The atom-mask argument must be type-checked before it is dereferenced.

The cross-SWIG-runtime workaround reads an object's `this` attribute and casts
it to a C++ pointer. Without a type check first, a wrong-typed mask is a
segfault, not an exception.
"""

from __future__ import annotations

import maptitude
import pytest

oechem = pytest.importorskip("openeye.oechem")
oegrid = pytest.importorskip("openeye.oegrid")


def _make_grid(dim: int, spacing: float):
    """Build a cubic grid centred on the origin.

    OESkewGrid has no extents-box constructor, so the dims and midpoint such a
    box derived are set explicitly.
    """
    grid = oegrid.OESkewGrid()
    assert grid.SetDim(dim, dim, dim)
    edge = dim * spacing
    assert grid.SetUnitCell(edge, edge, edge, 90.0, 90.0, 90.0, dim, dim, dim)
    assert grid.SetMid(0.0, 0.0, 0.0)
    return grid


def _make_uniform_grid():
    """A 17x17x17 grid at a 0.5 A node interval, every node 1.0."""
    grid = _make_grid(17, 0.5)
    ones = oechem.OEFloatArray(grid.GetSize())
    for i in range(grid.GetSize()):
        ones[i] = 1.0
    assert grid.SetValues(ones, grid.GetSize())
    return grid


@pytest.fixture()
def mol_and_grids():
    mol = oechem.OEGraphMol()
    atom = mol.NewAtom(6)
    mol.SetCoords(atom, (0.0, 0.0, 0.0))
    return mol, _make_uniform_grid(), _make_uniform_grid()


@pytest.mark.parametrize("bad_mask", ["not a predicate", 42, [1, 2, 3], {"a": 1}, object()])
def test_wrong_typed_builtin_mask_raises_rather_than_crashing(mol_and_grids, bad_mask) -> None:
    """Objects with no `this` attribute. These are the easy half."""
    mol, obs, calc = mol_and_grids
    with pytest.raises(TypeError):
        maptitude.rscc(mol, obs, 2.0, bad_mask, calc)


@pytest.mark.parametrize(
    "factory",
    [
        pytest.param(lambda: oechem.OEGraphMol(), id="OEGraphMol"),
        pytest.param(lambda: _make_grid(3, 1.0), id="OESkewGrid"),
        pytest.param(lambda: oechem.OEFloatArray([0.0, 0.0, 0.0]), id="OEFloatArray"),
        pytest.param(lambda: oechem.OEIsHeavy, id="OEIsHeavy-class-not-instance"),
    ],
)
def test_wrong_typed_swig_mask_raises_rather_than_crashing(mol_and_grids, factory) -> None:
    """The dangerous half: real SWIG proxies that DO have a `this` attribute.

    A builtin fails safely by accident — `_maptitude_extract_swig_ptr` finds no
    `this` attribute and returns NULL before it can pun anything. These objects
    have one, so the extractor happily reinterpret_casts an unrelated OpenEye
    C++ object to a predicate pointer. This is the case the isinstance check in
    the typemap exists to stop, and the only case that distinguishes the fixed
    build from the broken one.
    """
    mol, obs, calc = mol_and_grids
    with pytest.raises(TypeError):
        maptitude.rscc(mol, obs, 2.0, factory(), calc)


def test_none_mask_is_still_accepted(mol_and_grids) -> None:
    mol, obs, calc = mol_and_grids
    result = maptitude.rscc(mol, obs, 2.0, None, calc)
    assert result is not None


def test_a_real_predicate_is_still_accepted(mol_and_grids) -> None:
    mol, obs, calc = mol_and_grids
    mask = oechem.OEIsHeavy()
    result = maptitude.rscc(mol, obs, 2.0, mask, calc)
    assert result is not None


def test_forged_this_without_class_trickery_raises(mol_and_grids) -> None:
    """A plain Python object with a forged `this` borrowed from a real SWIG proxy.

    This is the pre-fix unsafe path: an object that has no `__class__` trickery
    but carries a `this` attribute stolen from a legitimate OpenEye object. The
    extractor will read it, but the type check must reject it before the cast.
    """
    mol, obs, calc = mol_and_grids
    victim = oechem.OEGraphMol()

    class ForgedThis:
        def __init__(self, donor):
            self.this = donor.this

    with pytest.raises(TypeError):
        maptitude.rscc(mol, obs, 2.0, ForgedThis(victim), calc)


@pytest.fixture()
def two_atom_mol_and_grids():
    """Two heavy atoms: carbon at origin, oxygen at (1.5, 0, 0)."""
    mol = oechem.OEGraphMol()
    carbon = mol.NewAtom(6)
    oxygen = mol.NewAtom(8)
    mol.SetCoords(carbon, (0.0, 0.0, 0.0))
    mol.SetCoords(oxygen, (1.5, 0.0, 0.0))
    return mol, _make_uniform_grid(), _make_uniform_grid()


def test_mask_is_actually_used(two_atom_mol_and_grids) -> None:
    """Verify the mask parameter filters atoms, not just whether the call succeeds.

    On a two-atom fixture (carbon and oxygen), mask=None includes both atoms,
    mask=OEHasAtomicNum(8) includes only oxygen, mask=OEHasAtomicNum(6) includes
    only carbon. Assert on the by_atom keys to prove the right atoms were selected.
    """
    mol, obs, calc = two_atom_mol_and_grids

    # mask=None should include both atoms (indices 0 and 1)
    result_none = maptitude.rscc(mol, obs, 2.0, None, calc)
    assert set(result_none.by_atom.keys()) == {0, 1}

    # mask=OEHasAtomicNum(8) should include only oxygen (index 1)
    result_oxygen = maptitude.rscc(mol, obs, 2.0, oechem.OEHasAtomicNum(8), calc)
    assert set(result_oxygen.by_atom.keys()) == {1}

    # mask=OEHasAtomicNum(6) should include only carbon (index 0)
    result_carbon = maptitude.rscc(mol, obs, 2.0, oechem.OEHasAtomicNum(6), calc)
    assert set(result_carbon.by_atom.keys()) == {0}
