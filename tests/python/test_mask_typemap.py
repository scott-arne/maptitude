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


@pytest.fixture()
def mol_and_grids():
    mol = oechem.OEGraphMol()
    atom = mol.NewAtom(6)
    mol.SetCoords(atom, (0.0, 0.0, 0.0))
    obs = oegrid.OEScalarGrid(oechem.OEDoubleArray([-4.0, -4.0, -4.0, 4.0, 4.0, 4.0]), 0.5)
    calc = oegrid.OEScalarGrid(oechem.OEDoubleArray([-4.0, -4.0, -4.0, 4.0, 4.0, 4.0]), 0.5)
    for i in range(obs.GetSize()):
        obs.SetValue(i, 1.0)
        calc.SetValue(i, 1.0)
    return mol, obs, calc


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
        pytest.param(
            lambda: oegrid.OEScalarGrid(
                oechem.OEDoubleArray([-1.0, -1.0, -1.0, 1.0, 1.0, 1.0]), 1.0
            ),
            id="OEScalarGrid",
        ),
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
