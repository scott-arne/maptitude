"""Behavior of the hand-written Python wrappers at the language boundary.

Covers what the wrappers own rather than delegate: they must not mutate the
options objects they are given, they must reject an argument the C++ side would
reject, and the errors they raise for a bad argument must be the documented
types rather than whatever a bare subscript or a SWIG overload happens to emit.
"""

from __future__ import annotations

import maptitude
import pytest

oechem = pytest.importorskip("openeye.oechem")
oegrid = pytest.importorskip("openeye.oegrid")

# Dims and midpoint of the grid the extents box [-4, -4, -4, 4, 4, 4] at a 0.5 A
# node interval used to produce. OESkewGrid has no extents-box constructor, so
# the geometry that box derived is set explicitly.
_DIM = 17
_SPACING = 0.5


def _make_grid():
    grid = oegrid.OESkewGrid()
    assert grid.SetDim(_DIM, _DIM, _DIM)
    edge = _DIM * _SPACING
    assert grid.SetUnitCell(edge, edge, edge, 90.0, 90.0, 90.0, _DIM, _DIM, _DIM)
    assert grid.SetMid(0.0, 0.0, 0.0)
    return grid


@pytest.fixture()
def scoring_inputs():
    mol = oechem.OEGraphMol()
    atom = mol.NewAtom(6)
    mol.SetCoords(atom, (0.0, 0.0, 0.0))
    obs = _make_grid()
    calc = _make_grid()
    ones = oechem.OEFloatArray(obs.GetSize())
    for i in range(obs.GetSize()):
        ones[i] = 1.0
    assert obs.SetValues(ones, obs.GetSize())
    assert calc.SetValues(ones, calc.GetSize())
    return mol, obs, calc


def test_rscc_does_not_mutate_caller_options(scoring_inputs) -> None:
    mol, obs, calc = scoring_inputs
    options = maptitude.RsccOptions()
    before = options.GetAtomRadiusMethod()
    maptitude.rscc(mol, obs, 2.0, calc_grid=calc, atom_radius="fixed", options=options)
    assert options.GetAtomRadiusMethod() == before


def test_rsr_does_not_mutate_caller_options(scoring_inputs) -> None:
    mol, obs, calc = scoring_inputs
    options = maptitude.RsrOptions()
    before = options.GetAtomRadiusMethod()
    maptitude.rsr(mol, obs, 2.0, calc_grid=calc, atom_radius="fixed", options=options)
    assert options.GetAtomRadiusMethod() == before


def test_coverage_does_not_mutate_caller_options(scoring_inputs) -> None:
    mol, obs, _ = scoring_inputs
    options = maptitude.CoverageOptions()
    options.SetSigma(2.0)
    maptitude.coverage(mol, obs, sigma=0.5, options=options)
    assert options.GetSigma() == pytest.approx(2.0)


def test_coverage_sigma_defaults_to_the_options_value(sigma_discriminating_inputs) -> None:
    """Omitting sigma must leave the options object's value in force.

    This does not test the sentinel regression (an explicit sigma=1.0 being
    swallowed). That is covered by test_coverage_accepts_an_explicit_sigma_of_one.
    """
    mol, obs = sigma_discriminating_inputs
    options = maptitude.CoverageOptions()
    options.SetSigma(3.0)
    with_options = maptitude.coverage(mol, obs, options=options)

    explicit = maptitude.coverage(mol, obs, sigma=3.0)
    assert with_options.overall == pytest.approx(explicit.overall)


@pytest.fixture()
def sigma_discriminating_inputs():
    """A grid on which coverage genuinely depends on sigma.

    `scoring_inputs` uses a uniform grid, where stddev is 0, the threshold is
    `mean` for every sigma, and coverage is 1.0 unconditionally. A sentinel test
    built on it passes with the bug still in place. This fixture instead fills a
    slab -1.0 <= x <= 1.0 with 1.0 and the rest with 0.0, which puts 1445 of the
    4913 voxels high: mean 0.294118, stddev 0.455645. An atom at the origin
    interpolates to exactly 1.0, so it is covered iff sigma < 1.549.
    """
    mol = oechem.OEGraphMol()
    atom = mol.NewAtom(6)
    mol.SetCoords(atom, (0.0, 0.0, 0.0))
    obs = _make_grid()
    slab = oechem.OEFloatArray(obs.GetSize())
    for i in range(obs.GetSize()):
        x, _, _ = obs.ElementToSpatialCoord(i)
        slab[i] = 1.0 if abs(x) <= 1.0 + 1e-9 else 0.0
    assert obs.SetValues(slab, obs.GetSize())
    return mol, obs


def test_the_sigma_fixture_actually_discriminates(sigma_discriminating_inputs) -> None:
    """Guard the guard: if this fails, the sentinel test below proves nothing."""
    mol, obs = sigma_discriminating_inputs
    at_one = maptitude.coverage(mol, obs, sigma=1.0)
    at_three = maptitude.coverage(mol, obs, sigma=3.0)
    assert at_one.overall == pytest.approx(1.0)
    assert at_three.overall == pytest.approx(0.0)


def test_coverage_accepts_an_explicit_sigma_of_one(sigma_discriminating_inputs) -> None:
    """sigma=1.0 was previously indistinguishable from omitting it."""
    mol, obs = sigma_discriminating_inputs
    options = maptitude.CoverageOptions()
    options.SetSigma(3.0)
    # An explicit 1.0 must override the options value, not be swallowed.
    # Old wrapper: `if sigma != 1.0` is False, so options' 3.0 wins -> 0.0.
    # Fixed wrapper: the explicit 1.0 reaches SetSigma -> 1.0.
    explicit_one = maptitude.coverage(mol, obs, sigma=1.0, options=options)
    assert explicit_one.overall == pytest.approx(1.0)


def test_rscc_rejects_wrong_options_class_on_override_path(scoring_inputs) -> None:
    """rscc must reject RsrOptions when atom_radius is set."""
    mol, obs, calc = scoring_inputs
    with pytest.raises(TypeError, match="options must be an RsccOptions"):
        maptitude.rscc(mol, obs, 2.0, calc_grid=calc, atom_radius="fixed", options=maptitude.RsrOptions())


def test_rsr_rejects_wrong_options_class_on_override_path(scoring_inputs) -> None:
    """rsr must reject RsccOptions when atom_radius is set."""
    mol, obs, calc = scoring_inputs
    with pytest.raises(TypeError, match="options must be an RsrOptions"):
        maptitude.rsr(mol, obs, 2.0, calc_grid=calc, atom_radius="fixed", options=maptitude.RsccOptions())


def test_coverage_rejects_wrong_options_class_on_override_path(scoring_inputs) -> None:
    """coverage must reject QScoreOptions when sigma is set."""
    mol, obs, _ = scoring_inputs
    with pytest.raises(TypeError, match="options must be a CoverageOptions"):
        maptitude.coverage(mol, obs, sigma=2.0, options=maptitude.QScoreOptions())


def test_non_override_path_still_rejects_wrong_options(scoring_inputs) -> None:
    """SWIG's own type check must remain when the wrapper does not copy."""
    mol, obs, calc = scoring_inputs
    # Non-override path: no atom_radius, so wrapper does not call _copy_rscc_options
    with pytest.raises(TypeError, match="Wrong number or type of arguments"):
        maptitude.rscc(mol, obs, 2.0, calc_grid=calc, options=maptitude.RsrOptions())


# ---------------------------------------------------------------------------
# Atom-radius method: the string and enum paths must agree, and both must say
# what went wrong.
# ---------------------------------------------------------------------------


def test_rscc_rejects_the_adaptive_radius_enum(scoring_inputs) -> None:
    """rscc has no adaptive radius model; it used to return the binned score.

    RuntimeError is how std::invalid_argument surfaces, which the README names
    as the channel for option-value errors.
    """
    mol, obs, calc = scoring_inputs
    options = maptitude.RsccOptions()
    options.SetAtomRadiusMethod(maptitude.AtomRadius.ADAPTIVE)
    with pytest.raises(RuntimeError, match="no adaptive atom-radius model"):
        maptitude.rscc(mol, obs, 2.0, calc_grid=calc, options=options)


def test_rscc_rejects_the_adaptive_radius_string(scoring_inputs) -> None:
    """The string path must reject exactly what the enum path rejects."""
    mol, obs, calc = scoring_inputs
    with pytest.raises(ValueError, match="rscc does not support atom_radius='adaptive'"):
        maptitude.rscc(mol, obs, 2.0, calc_grid=calc, atom_radius="adaptive")


@pytest.mark.parametrize("metric", ["rscc", "rsr"])
def test_an_unknown_radius_string_raises_value_error(metric, scoring_inputs) -> None:
    """A misspelling used to raise a bare KeyError('vdw').

    KeyError is outside the typed hierarchy the README presents as the complete
    set to handle, and it carried no list of what would have worked.
    """
    mol, obs, calc = scoring_inputs
    with pytest.raises(ValueError) as excinfo:
        getattr(maptitude, metric)(mol, obs, 2.0, calc_grid=calc, atom_radius="vdw")
    message = str(excinfo.value)
    assert "'vdw'" in message
    for accepted in ("'binned'", "'fixed'", "'scaled'"):
        assert accepted in message


@pytest.mark.parametrize("name", ["fixed", "scaled", "binned"])
def test_rscc_still_accepts_its_three_radius_strings(name, scoring_inputs) -> None:
    """The neutrality half: rejecting 'adaptive' must not narrow the rest."""
    mol, obs, calc = scoring_inputs
    maptitude.rscc(mol, obs, 2.0, calc_grid=calc, atom_radius=name)


def test_rsr_still_accepts_the_adaptive_radius_string(scoring_inputs) -> None:
    """rsr has a real adaptive branch and keeps it as its default."""
    mol, obs, calc = scoring_inputs
    maptitude.rsr(mol, obs, 2.0, calc_grid=calc, atom_radius="adaptive")


# ---------------------------------------------------------------------------
# fc_density symops normalization
# ---------------------------------------------------------------------------


@pytest.fixture()
def fc_inputs(scoring_inputs):
    mol, obs, _ = scoring_inputs
    return mol, obs, maptitude.UnitCell(20.0, 20.0, 20.0, 90.0, 90.0, 90.0)


def _grid_values(grid):
    # GetValues hands back a copy, so read it once rather than per node.
    values = grid.GetValues()
    return [values[i] for i in range(grid.GetSize())]


def test_fc_density_accepts_a_list_of_operator_strings(fc_inputs) -> None:
    """The documented shape. It used to reach SymOpVector and die there.

    The error was a raw SWIG overload dump listing four C++ vector constructors,
    which is neither a typed maptitude exception nor a usable description of the
    argument that was wrong.
    """
    mol, obs, cell = fc_inputs
    from_list = maptitude.fc_density(mol, obs, 2.0, cell, symops=["x,y,z", "-x,-y,z"])
    from_string = maptitude.fc_density(mol, obs, 2.0, cell, symops="x,y,z;-x,-y,z")
    assert _grid_values(from_list) == _grid_values(from_string)


def test_fc_density_still_accepts_symop_objects(fc_inputs) -> None:
    """The other documented shape, and the one the deleted branch guarded."""
    mol, obs, cell = fc_inputs
    parsed = list(maptitude.SymOp.ParseAll("x,y,z"))
    from_objects = maptitude.fc_density(mol, obs, 2.0, cell, symops=parsed)
    from_string = maptitude.fc_density(mol, obs, 2.0, cell, symops="x,y,z")
    assert _grid_values(from_objects) == _grid_values(from_string)


def test_fc_density_rejects_a_non_iterable_symops(fc_inputs) -> None:
    """An int used to be accepted: SymOpVector(42) builds 42 identity operators."""
    mol, obs, cell = fc_inputs
    with pytest.raises(TypeError, match="symops must be a string"):
        maptitude.fc_density(mol, obs, 2.0, cell, symops=42)


def test_fc_density_rejects_symops_entries_of_the_wrong_type(fc_inputs) -> None:
    mol, obs, cell = fc_inputs
    with pytest.raises(TypeError, match="symops entries must be SymOp or str, not int"):
        maptitude.fc_density(mol, obs, 2.0, cell, symops=[1, 2])


def test_wrap_and_pad_grid_never_returns_none(scoring_inputs) -> None:
    """The wrapper substitutes the original grid for the C++ nullptr.

    The README used to show `if padded is not None:`, a branch that can never
    run. Pinning the contract here is what makes that documentation checkable.
    """
    mol, obs, _ = scoring_inputs
    # The cell must be a whole number of the grid's node intervals, so derive each
    # edge from that axis's own interval rather than from one figure for the whole
    # grid:
    # wrap_and_pad_grid checks each edge against its axis independently. The atom
    # is at the origin and the grid is centred there, so no padding is needed and
    # the C++ function returns nullptr.
    gp = maptitude.get_grid_params(obs)
    padded = maptitude.wrap_and_pad_grid(
        mol=mol,
        grid=obs,
        cell_a=gp.x_dim * gp.x_spacing,
        cell_b=gp.y_dim * gp.y_spacing,
        cell_c=gp.z_dim * gp.z_spacing,
    )
    assert padded is obs
