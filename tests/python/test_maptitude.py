"""
Integration tests for maptitude Python bindings.

These tests verify that the SWIG bindings work correctly and that
OpenEye objects can be passed between Python and C++.
"""

import pytest


def _has_openeye():
    """Check if OpenEye toolkits are available."""
    try:
        from openeye import oechem
        return oechem.OEChemIsLicensed()
    except ImportError:
        return False


def test_import():
    """Verify the package can be imported."""
    import maptitude
    assert hasattr(maptitude, "__version__")
    assert maptitude.__version__ == "0.1.5"


def test_unit_cell_creation():
    """Test UnitCell construction and basic properties."""
    from maptitude import UnitCell

    cell = UnitCell(50.0, 60.0, 70.0, 90.0, 90.0, 90.0)
    assert cell.a == 50.0
    assert cell.b == 60.0
    assert cell.c == 70.0
    assert cell.alpha == 90.0
    assert cell.beta == 90.0
    assert cell.gamma == 90.0


def test_unit_cell_volume():
    """Test orthorhombic unit cell volume."""
    from maptitude import UnitCell

    cell = UnitCell(50.0, 60.0, 70.0, 90.0, 90.0, 90.0)
    assert abs(cell.Volume() - 210000.0) < 1e-6


def test_unit_cell_coordinate_roundtrip():
    """Test Cartesian <-> fractional coordinate conversion."""
    from maptitude import UnitCell

    cell = UnitCell(50.0, 60.0, 70.0, 90.0, 90.0, 90.0)
    frac = cell.CartesianToFractional(25.0, 30.0, 35.0)
    cart = cell.FractionalToCartesian(frac[0], frac[1], frac[2])

    assert abs(cart[0] - 25.0) < 1e-8
    assert abs(cart[1] - 30.0) < 1e-8
    assert abs(cart[2] - 35.0) < 1e-8


def test_symop_parse_identity():
    """Test parsing identity symmetry operator."""
    from maptitude import parse_symop

    op = parse_symop("x,y,z")
    result = op.Apply(0.25, 0.5, 0.75)
    assert abs(result[0] - 0.25) < 1e-10
    assert abs(result[1] - 0.5) < 1e-10
    assert abs(result[2] - 0.75) < 1e-10


def test_symop_parse_with_translation():
    """Test parsing symmetry operator with fractional translation."""
    from maptitude import parse_symop

    op = parse_symop("-x,y+1/2,-z")
    result = op.Apply(0.25, 0.3, 0.1)
    assert abs(result[0] - (-0.25)) < 1e-10
    assert abs(result[1] - 0.8) < 1e-10
    assert abs(result[2] - (-0.1)) < 1e-10


def test_parse_symops_multiple():
    """Test parsing multiple symmetry operators."""
    from maptitude import parse_symops

    ops = parse_symops("x,y,z\n-x,y+1/2,-z+1/2")
    assert len(ops) == 2


def test_residue_creation():
    """Test Residue construction."""
    from maptitude import Residue

    r = Residue("ALA", 123, "A", " ")
    assert r.name == "ALA"
    assert r.number == 123
    assert r.chain == "A"


def test_residue_hash():
    """Test that Residue can be used as a dict key."""
    from maptitude import Residue

    r1 = Residue("ALA", 1, "A")
    r2 = Residue("GLY", 2, "A")
    d = {r1: 0.95, r2: 0.88}
    assert len(d) == 2


def test_qscore_options_defaults():
    """Test QScoreOptions default values."""
    from maptitude import QScoreOptions

    opts = QScoreOptions()
    assert opts.sigma == 0.6
    assert opts.radial_step == 0.5
    assert opts.max_radius == 2.0
    assert opts.num_points == 8
    assert opts.normalize_map is True
    assert opts.isolate_points is True


def test_mapop_enum():
    """Test MapOp enum values."""
    from maptitude import MapOp

    assert MapOp.ADD != MapOp.SUBTRACT
    assert MapOp.MIN != MapOp.MAX


def test_density_score_result_tostring():
    """Test DensityScoreResult ToString method."""
    from maptitude import DensityScoreResult

    result = DensityScoreResult()
    result.overall = 0.95
    s = str(result)
    assert "0.950" in s
    assert "residues=" in s
    assert "atoms=" in s


def test_scattering_factor_carbon():
    """Test scattering factor lookup for carbon."""
    from maptitude import get_scattering_factors

    coeffs = get_scattering_factors(6, 0)
    assert coeffs is not None
    f0 = coeffs.Evaluate(0.0)
    assert abs(f0 - 6.0) < 0.3


def test_scattering_factor_table_size():
    """Test that the full scattering factor table is loaded."""
    from maptitude import get_scattering_factor_table

    table, count = get_scattering_factor_table()
    assert count == 209


def test_scattering_factor_fallback():
    """Test that unknown charge falls back to neutral."""
    from maptitude import get_scattering_factors

    neutral = get_scattering_factors(6, 0)
    fallback = get_scattering_factors(6, 5)
    assert neutral is not None
    assert fallback is not None
    assert abs(neutral.Evaluate(0.1) - fallback.Evaluate(0.1)) < 1e-10


@pytest.mark.skipif(
    not _has_openeye(),
    reason="OpenEye toolkits not available"
)
def test_density_scorer_rscc():
    """Test RSCC scoring with a simple molecule and grid."""
    from openeye import oechem, oegrid
    from maptitude import rscc

    # Create a simple molecule with coordinates
    mol = oechem.OEGraphMol()
    atom = mol.NewAtom(6)  # Carbon
    coords = oechem.OEFloatArray([0.0, 0.0, 0.0])
    mol.SetCoords(atom, coords)

    # Create observed and calculated grids
    obs_grid = oegrid.OEScalarGrid(10, 10, 10, -5.0, -5.0, -5.0, 1.0)
    calc_grid = oegrid.OEScalarGrid(10, 10, 10, -5.0, -5.0, -5.0, 1.0)

    # Fill both grids with a Gaussian centered at the atom
    import math
    for i in range(obs_grid.GetSize()):
        fx, fy, fz = obs_grid.ElementToSpatialCoord(i)
        r2 = fx**2 + fy**2 + fz**2
        val = math.exp(-r2 / 2.0)
        obs_grid[i] = val
        calc_grid[i] = val

    result = rscc(mol, obs_grid, 2.0, calc_grid=calc_grid)
    assert hasattr(result, "overall")


def _has_openeye():
    """Check if OpenEye toolkits are available."""
    try:
        from openeye import oechem, oegrid
        return True
    except ImportError:
        return False
