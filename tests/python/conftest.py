"""Shared fixtures and configuration for maptitude Python tests."""

import os
import sys

import pytest

# The map-loading helpers (MRC/CCP4 readers, symop extraction) are shared with
# the benchmark suite. Put the benchmarks directory on sys.path so tests import
# the single canonical implementation from ``helpers`` rather than duplicating it.
_BENCH_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "benchmarks")
)
if _BENCH_DIR not in sys.path:
    sys.path.insert(0, _BENCH_DIR)

pytest.importorskip("openeye.oechem", reason="OpenEye Toolkits not installed")


@pytest.fixture
def aspirin_mol():
    """Create an aspirin molecule (C9H8O4) for testing."""
    from openeye import oechem

    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, "CC(=O)OC1=CC=CC=C1C(=O)O")
    return mol


@pytest.fixture
def ethanol_mol():
    """Create an ethanol molecule (C2H6O) for testing."""
    from openeye import oechem

    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, "CCO")
    return mol
