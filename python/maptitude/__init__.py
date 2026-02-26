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
    RSCCOptions,
    DensityScoreResult,
    GridParams,
    MapOp_Add,
    MapOp_Subtract,
    MapOp_Min,
    MapOp_Max,
    AtomRadius_Fixed,
    AtomRadius_Scaled,
    AtomRadius_Binned,
    RadialSampling_Fixed,
    RadialSampling_Adaptive,

    # Classes
    DensityCalculator,

    # Density scoring functions
    RSCC,
    RSR,
    QScore,
    EDIAm,
    Coverage,

    # Grid utility functions
    GetGridParams,
    InterpolateDensity,
    InterpolateDensityBatch,
    InterpolateDensityPeriodic,
    InterpolateDensityPeriodicBatch,
    GetAtomGridPoints,
    ScaleMap,
    CombineMaps,
    DiffToCalc,
    WrapAndPadGrid,

    # Convenience functions
    fc_density,
    rscc,
    rsr,
    qscore,
    ediam,
    coverage,
    parse_symop,
    parse_symops,
    scale_map,
    combine_maps,
    diff_to_calc,
    wrap_and_pad_grid,

    # Container types
    DoubleVector,
    UnsignedIntVector,
    SymOpVector,
)


class MapOp:
    """Grid combination operations."""
    Add = MapOp_Add
    Subtract = MapOp_Subtract
    Min = MapOp_Min
    Max = MapOp_Max


class AtomRadius:
    """Atom radius methods for RSCC scoring."""
    Fixed = AtomRadius_Fixed
    Scaled = AtomRadius_Scaled
    Binned = AtomRadius_Binned


class RadialSampling:
    """Radial sampling strategies for Q-score computation."""
    Fixed = RadialSampling_Fixed
    Adaptive = RadialSampling_Adaptive


__all__ = [
    "__version__",
    "__version_info__",
    # Core types
    "Residue",
    "UnitCell",
    "SymOp",
    "QScoreOptions",
    "RSCCOptions",
    "AtomRadius",
    "RadialSampling",
    "DensityScoreResult",
    "GridParams",
    "MapOp",
    # Classes
    "DensityCalculator",
    # Density scoring functions
    "RSCC",
    "RSR",
    "QScore",
    "EDIAm",
    "Coverage",
    # Grid utilities
    "GetGridParams",
    "InterpolateDensity",
    "InterpolateDensityBatch",
    "InterpolateDensityPeriodic",
    "InterpolateDensityPeriodicBatch",
    "GetAtomGridPoints",
    "ScaleMap",
    "CombineMaps",
    "DiffToCalc",
    "WrapAndPadGrid",
    # Convenience functions
    "fc_density",
    "rscc",
    "rsr",
    "qscore",
    "ediam",
    "coverage",
    "parse_symop",
    "parse_symops",
    "scale_map",
    "combine_maps",
    "diff_to_calc",
    "wrap_and_pad_grid",
]
