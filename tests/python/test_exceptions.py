"""Assert the typed exception hierarchy survives SWIG regeneration.

Each name has to traverse five layers: a C global, a %init binding on the
extension module, a %pythoncode alias onto the SWIG proxy, the %exception
raise site, and the package re-export. A break in any one of them shows up
here.
"""

from __future__ import annotations

import subprocess
import sys

import maptitude
import pytest

EXCEPTION_NAMES = ["StructureError", "GridError", "SymOpError", "CellError"]
ALL_EXCEPTION_NAMES = EXCEPTION_NAMES + ["MaptitudeError"]


def test_base_exception_is_exported() -> None:
    assert issubclass(maptitude.MaptitudeError, Exception)


@pytest.mark.parametrize("name", EXCEPTION_NAMES)
def test_exception_is_exported_and_parented(name: str) -> None:
    cls = getattr(maptitude, name)
    assert issubclass(cls, maptitude.MaptitudeError)
    assert issubclass(cls, Exception)
    assert cls.__base__ is maptitude.MaptitudeError


@pytest.mark.parametrize("name", ALL_EXCEPTION_NAMES)
def test_exception_is_in_dunder_all(name: str) -> None:
    assert name in maptitude.__all__


@pytest.mark.parametrize("name", ALL_EXCEPTION_NAMES)
def test_exception_is_not_a_value_error(name: str) -> None:
    """The point of the task: these no longer collapse into ValueError.

    Pinning this stops anyone reintroducing the old behavior as a
    backwards-compatibility shim, which the project has ruled out.
    """
    assert not issubclass(getattr(maptitude, name), ValueError)


@pytest.mark.parametrize("name", ALL_EXCEPTION_NAMES)
def test_all_three_layers_expose_the_same_class(name: str) -> None:
    """The class is created on the extension module and must survive two hops.

    %pythoncode aliases it onto the SWIG proxy and __init__.py re-exports from
    the proxy. Comparing only the package against the proxy is tautological --
    the package imports from the proxy. Pinning the extension too is what
    catches an alias typo such as ``CellError = _maptitude.GridError``, which
    would otherwise pass every other test in this file until Task 8 gives
    CellError a live throw.
    """
    from maptitude import maptitude as proxy

    extension = getattr(proxy, name)
    assert getattr(maptitude, name) is extension
    assert getattr(proxy._maptitude, name) is extension


@pytest.mark.parametrize("name", ALL_EXCEPTION_NAMES)
def test_exception_module_is_the_public_package(name: str) -> None:
    """The dotted name given to PyErr_NewException is what tracebacks show."""
    assert getattr(maptitude, name).__module__ == "maptitude"


def test_concrete_exceptions_are_distinct_classes() -> None:
    """Two names bound to one class would silently restore the ambiguity."""
    classes = {getattr(maptitude, name) for name in EXCEPTION_NAMES}
    assert len(classes) == len(EXCEPTION_NAMES)


def test_structure_error_is_raised_with_its_own_type() -> None:
    """A molecule with no scorable atoms is a StructureError, not a GridError.

    The explicit calc_grid isolates this from validation-order changes -- without
    it, the test would fail as GridError if the missing-grid check ran first.
    """
    oechem = pytest.importorskip("openeye.oechem")
    oegrid = pytest.importorskip("openeye.oegrid")

    grid = oegrid.OEScalarGrid(oechem.OEDoubleArray([-4.0, -4.0, -4.0, 4.0, 4.0, 4.0]), 0.5)

    with pytest.raises(maptitude.StructureError) as excinfo:
        maptitude.rscc(oechem.OEGraphMol(), grid, 2.0, calc_grid=grid)
    assert type(excinfo.value) is maptitude.StructureError
    assert "No scorable heavy atoms" in str(excinfo.value)


def test_grid_error_is_raised_with_its_own_type() -> None:
    """rscc without a calc_grid raises GridError, not a bare ValueError."""
    oechem = pytest.importorskip("openeye.oechem")
    oegrid = pytest.importorskip("openeye.oegrid")

    mol = oechem.OEGraphMol()
    atom = mol.NewAtom(6)
    mol.SetCoords(atom, (0.0, 0.0, 0.0))
    grid = oegrid.OEScalarGrid(oechem.OEDoubleArray([-4.0, -4.0, -4.0, 4.0, 4.0, 4.0]), 0.5)

    with pytest.raises(maptitude.GridError) as excinfo:
        maptitude.rscc(mol, grid, 2.0)
    assert type(excinfo.value) is maptitude.GridError
    assert "calc_grid is required" in str(excinfo.value)


def test_symop_error_is_raised_with_its_own_type() -> None:
    """Each typed arm of %exception needs its own live raise.

    Without this, transposing two arms so the wrong global reaches
    PyErr_SetString would leave the whole file green.
    """
    with pytest.raises(maptitude.SymOpError) as excinfo:
        maptitude.parse_symop("not a symop")
    assert type(excinfo.value) is maptitude.SymOpError
    assert "exactly 3 components" in str(excinfo.value)


# Run in a child interpreter: re-importing in this process would leave the
# pytest session holding a stale maptitude and corrupt every later test.
_REIMPORT_SCRIPT = """
import sys

import maptitude

held = maptitude.SymOpError
held_base = maptitude.MaptitudeError

for name in [key for key in sys.modules if "maptitude" in key]:
    del sys.modules[name]

import maptitude as reimported

assert reimported.SymOpError is held, "SymOpError was replaced by the re-import"
assert reimported.MaptitudeError is held_base, "MaptitudeError base was replaced"

try:
    reimported.parse_symop("not a symop")
except held as exc:
    assert type(exc) is held, "raised %r, expected %r" % (type(exc), held)
else:
    raise AssertionError("parse_symop did not raise")
"""


def test_exception_classes_survive_a_re_import() -> None:
    """A second run of the %init block must reuse the cached classes.

    SWIG emits %init into a Py_mod_exec slot, and the generated PyModuleDef
    sets ``m_size`` to 0, so CPython may run that slot more than once per
    process -- dropping the ``maptitude`` keys from ``sys.modules`` and
    importing again is enough to do it. While %init minted a fresh class on
    every run it left already-imported code holding superseded objects, so an
    ``except SymOpError`` bound before the re-import stopped catching what
    %exception raised. Identity is checked first, then a live raise, because
    the caught-or-not behavior is the part that actually broke.
    """
    result = subprocess.run(
        [sys.executable, "-c", _REIMPORT_SCRIPT],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, (
        f"re-import child exited {result.returncode}:\n{result.stderr}"
    )
