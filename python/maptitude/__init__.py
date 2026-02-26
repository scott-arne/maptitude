"""
Maptitude - High-performance crystallographic electron density computation and scoring

This package provides C++ implementations of electron density calculation
and real-space scoring metrics with Python bindings via SWIG.

Example usage::

    from openeye import oechem, oegrid
    from maptitude import fc_density, rscc, rsr, qscore, UnitCell, parse_symops

    # Load structure and map
    mol = oechem.OEGraphMol()
    oechem.OEReadMolecule(oechem.oemolistream("model.pdb"), mol)
    obs_grid = oegrid.OEScalarGrid()
    oegrid.OEReadGrid(obs_grid, "2fofc.map")

    # Define crystal form
    cell = UnitCell(50.0, 60.0, 70.0, 90.0, 90.0, 90.0)
    symops = parse_symops("x,y,z\\n-x,y+1/2,-z+1/2")

    # Compute model density
    calc_grid = fc_density(mol, obs_grid, 2.0, cell, symops=symops)

    # Score fit
    result = rscc(mol, obs_grid, 2.0, calc_grid=calc_grid)
    print(f"Overall RSCC: {result.overall}")
"""

import warnings

# Version info
__version__ = "0.1.0"
__version_info__ = (0, 1, 0)


def _check_openeye_version():
    """Check that the OpenEye version matches what was used at build time."""
    try:
        from . import _build_info
    except ImportError:
        # Build info not available (development install)
        return

    # Only check if using shared/dynamic linking
    if getattr(_build_info, 'OPENEYE_LIBRARY_TYPE', 'STATIC') != 'SHARED':
        return

    build_version = getattr(_build_info, 'OPENEYE_BUILD_VERSION', None)
    if not build_version:
        return

    try:
        from openeye import oechem
        runtime_version = oechem.OEToolkitsGetRelease()
        if runtime_version and build_version:
            build_parts = build_version.split('.')[:3]
            runtime_parts = runtime_version.split('.')[:3]
            if build_parts != runtime_parts:
                warnings.warn(
                    f"OpenEye version mismatch: maptitude was built with OpenEye Toolkits "
                    f"{build_version} but runtime has OpenEye Toolkits {runtime_version}. "
                    f"This may cause compatibility issues.",
                    RuntimeWarning
                )
    except ImportError:
        warnings.warn(
            "openeye-toolkits package not found. "
            "This wheel requires openeye-toolkits to be installed. "
            "Install with: pip install openeye-toolkits",
            ImportWarning
        )


# Check OpenEye version on import
_check_openeye_version()

from .maptitude import (
    # Core types
    Residue,
    UnitCell,
    SymOp,
    QScoreOptions,
    RsccOptions,
    RsrOptions,
    CoverageOptions,
    DensityScoreResult,
    GridParams,
    MapOp_ADD,
    MapOp_SUBTRACT,
    MapOp_MIN,
    MapOp_MAX,
    AtomRadius_FIXED,
    AtomRadius_SCALED,
    AtomRadius_BINNED,
    AtomRadius_ADAPTIVE,
    RadialSampling_FIXED,
    RadialSampling_ADAPTIVE,

    # Classes
    DensityCalculator,

    # Density scoring functions
    rscc,
    rsr,
    qscore,
    ediam,
    coverage,

    # Grid utility functions
    get_grid_params,
    interpolate_density,
    interpolate_density_batch,
    interpolate_density_periodic,
    interpolate_density_periodic_batch,
    get_atom_grid_points,
    scale_map,
    combine_maps,
    diff_to_calc,
    wrap_and_pad_grid,

    # Scattering factor types and functions
    CromerMannCoeffs,
    ScatteringFactorEntry,
    get_scattering_factors,
    get_scattering_factor_table,

    # Convenience functions
    fc_density,
    parse_symop,
    parse_symops,

    # Container types
    DoubleVector,
    UnsignedIntVector,
    SymOpVector,
)


class MapOp:
    """Grid combination operations."""
    ADD = MapOp_ADD
    SUBTRACT = MapOp_SUBTRACT
    MIN = MapOp_MIN
    MAX = MapOp_MAX


class AtomRadius:
    """Atom radius methods for density scoring."""
    FIXED = AtomRadius_FIXED
    SCALED = AtomRadius_SCALED
    BINNED = AtomRadius_BINNED
    ADAPTIVE = AtomRadius_ADAPTIVE


class RadialSampling:
    """Radial sampling strategies for Q-score computation."""
    FIXED = RadialSampling_FIXED
    ADAPTIVE = RadialSampling_ADAPTIVE


__all__ = [
    "__version__",
    "__version_info__",
    # Core types
    "Residue",
    "UnitCell",
    "SymOp",
    "QScoreOptions",
    "RsccOptions",
    "RsrOptions",
    "CoverageOptions",
    "AtomRadius",
    "RadialSampling",
    "DensityScoreResult",
    "GridParams",
    "MapOp",
    # Classes
    "DensityCalculator",
    # Density scoring functions
    "rscc",
    "rsr",
    "qscore",
    "ediam",
    "coverage",
    # Grid utilities
    "get_grid_params",
    "interpolate_density",
    "interpolate_density_batch",
    "interpolate_density_periodic",
    "interpolate_density_periodic_batch",
    "get_atom_grid_points",
    "scale_map",
    "combine_maps",
    "diff_to_calc",
    "wrap_and_pad_grid",
    # Scattering factors
    "get_scattering_factors",
    "get_scattering_factor_table",
    # Convenience functions
    "fc_density",
    "parse_symop",
    "parse_symops",
]
