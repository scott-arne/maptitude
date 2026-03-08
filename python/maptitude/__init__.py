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

import os
import re
import warnings

# Version info
__version__ = "0.1.1"
__version_info__ = (0, 1, 1)


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


def _ensure_library_compat():
    """Create compatibility symlinks when OpenEye library versions differ from build time.

    When the package is built with shared OpenEye libraries, the compiled extension
    records exact versioned library filenames (e.g., liboechem-4.3.0.1.dylib).
    If the user upgrades openeye-toolkits, these filenames change (e.g., to
    liboechem-4.3.0.2.dylib) and the dynamic linker fails to load the extension.

    This function detects version mismatches and creates symlinks from the expected
    (build-time) library names to the actual (runtime) library files in the package
    directory, which is on the extension's rpath (@loader_path / $ORIGIN).
    """
    try:
        from . import _build_info
    except ImportError:
        return False

    if getattr(_build_info, 'OPENEYE_LIBRARY_TYPE', 'STATIC') != 'SHARED':
        return False

    expected_libs = getattr(_build_info, 'OPENEYE_EXPECTED_LIBS', [])
    if not expected_libs:
        return False

    try:
        from openeye import libs
        oe_lib_dir = libs.FindOpenEyeDLLSDirectory()
    except (ImportError, Exception):
        return False

    if not os.path.isdir(oe_lib_dir):
        return False

    pkg_dir = os.path.dirname(__file__)
    created_any = False

    for expected_name in expected_libs:
        # Check if the expected library exists in the OpenEye lib directory
        if os.path.exists(os.path.join(oe_lib_dir, expected_name)):
            continue

        # Check for an existing symlink in our package directory
        symlink_path = os.path.join(pkg_dir, expected_name)
        if os.path.islink(symlink_path):
            if os.path.exists(symlink_path):
                continue  # Valid symlink already exists
            # Stale symlink (target was removed/upgraded) - remove it
            try:
                os.unlink(symlink_path)
            except OSError:
                continue
        elif os.path.exists(symlink_path):
            continue  # Regular file exists

        # Extract the base library name (e.g., "liboechem" from "liboechem-4.3.0.1.dylib")
        # Handles both dash-versioned (liboechem-4.3.0.1.dylib) and
        # dot-versioned (libzstd.1.dylib) naming conventions
        match = re.match(r'(lib\w+?)(-[\d.]+)?(\.[\d.]*\w+)', expected_name)
        if not match:
            continue
        base_name = match.group(1)

        # Find the actual library with a potentially different version
        actual_path = None
        for f in os.listdir(oe_lib_dir):
            if f.startswith(base_name + '-') or f.startswith(base_name + '.'):
                actual_path = os.path.join(oe_lib_dir, f)
                break

        if actual_path:
            symlink_path = os.path.join(pkg_dir, expected_name)
            try:
                os.symlink(actual_path, symlink_path)
                created_any = True
            except OSError:
                pass

    return created_any


# Create compatibility symlinks before loading the C extension
_ensure_library_compat()

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
