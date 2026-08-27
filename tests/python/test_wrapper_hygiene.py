"""The Python wrappers must not mutate the options objects they are given."""

from __future__ import annotations

import maptitude
import pytest

oechem = pytest.importorskip("openeye.oechem")
oegrid = pytest.importorskip("openeye.oegrid")


@pytest.fixture()
def scoring_inputs():
    mol = oechem.OEGraphMol()
    atom = mol.NewAtom(6)
    mol.SetCoords(atom, (0.0, 0.0, 0.0))
    obs = oegrid.OEScalarGrid(oechem.OEDoubleArray([-4.0, -4.0, -4.0, 4.0, 4.0, 4.0]), 0.5)
    calc = oegrid.OEScalarGrid(oechem.OEDoubleArray([-4.0, -4.0, -4.0, 4.0, 4.0, 4.0]), 0.5)
    for i in range(obs.GetSize()):
        obs.SetValue(i, 1.0)
        calc.SetValue(i, 1.0)
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
    obs = oegrid.OEScalarGrid(oechem.OEDoubleArray([-4.0, -4.0, -4.0, 4.0, 4.0, 4.0]), 0.5)
    for i in range(obs.GetSize()):
        x, _, _ = obs.ElementToSpatialCoord(i)
        obs.SetValue(i, 1.0 if abs(x) <= 1.0 + 1e-9 else 0.0)
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
