# Maptitude

High-performance crystallographic electron density computation and scoring for Python.

Maptitude provides C++ implementations of model density calculation and five real-space scoring metrics, accessible
through a Python API that works directly with OpenEye Toolkits molecule and grid objects.

## Features

This project provides the following for both C++ and Python:

- **Five scoring metrics** for evaluating how well a molecular model fits an experimental electron density map:
  - **RSCC**: Real-Space Correlation Coefficient (Pearson correlation between observed and calculated density)
  - **RSR**: Real-Space R-Factor (residual between observed and calculated density)
  - **Q-score**: Radial density profile analysis (Pintilie et al., 2020)
  - **EDIAm**: Electron Density Index Averaged, Modified (density at atom centers and bond midpoints)
  - **Coverage**: Fraction of atoms observed above a density threshold
- **Model density calculation** via Fourier synthesis with Cromer-Mann scattering factors, symmetry expansion, bulk
  solvent correction, and per-shell amplitude scaling
- **Per-residue and per-atom resolution** for all metrics
- **Grid utilities** including trilinear interpolation, periodic-aware sampling, map scaling, and map combination
- **OpenMP parallelism** for structure factor accumulation

## Requirements

- Python 3.10 or later
- OpenEye Toolkits 2025.2 or later (with a valid license)

## Installation

Install from a wheel:

```bash
pip install maptitude
```

Building from source requires CMake 3.21+, SWIG 4.0+, FFTW3, and the OpenEye C++ Toolkits.

## Quick Start

```python
from openeye import oechem, oegrid
from maptitude import fc_density, rscc, UnitCell, parse_symops

# Load structure and observed density map
mol = oechem.OEGraphMol()
oechem.OEReadMolecule(oechem.oemolistream("model.pdb"), mol)

obs_grid = oegrid.OEScalarGrid()
oegrid.OEReadGrid(obs_grid, "2fofc.map")

# Score the fit
result = rscc(mol, obs_grid, 2.0)
print(f"Overall RSCC: {result.overall:.3f}")
```

## Usage

### Scoring Metrics

All scoring functions return a `DensityScoreResult` with an overall score and per-residue and per-atom breakdowns.

#### RSCC (Real-Space Correlation Coefficient)

```python
from maptitude import rscc, RsccOptions, AtomRadius

result = rscc(mol, obs_grid, resolution, calc_grid=calc_grid)

# With custom options
opts = RsccOptions()
opts.atom_radius_method = AtomRadius.ADAPTIVE
result = rscc(mol, obs_grid, resolution, calc_grid=calc_grid, options=opts)
```

#### RSR (Real-Space R-Factor)

```python
from maptitude import rsr

result = rsr(mol, obs_grid, resolution, calc_grid=calc_grid)
```

#### Q-score

```python
from maptitude import qscore, QScoreOptions

# Uses only the observed map (no calc_grid needed)
result = qscore(mol, obs_grid, resolution)

# Customize sampling parameters
opts = QScoreOptions()
opts.sigma = 0.6  # Gaussian width (Angstroms)
opts.num_points = 8  # Points per radial shell
opts.max_radius = 2.0  # Maximum sampling radius (Angstroms)
result = qscore(mol, obs_grid, resolution, options=opts)
```

#### EDIAm

```python
from maptitude import ediam

result = ediam(mol, obs_grid, resolution)
```

#### Coverage

```python
from maptitude import coverage, CoverageOptions

result = coverage(mol, obs_grid)

# Adjust the density threshold
opts = CoverageOptions()
opts.sigma = 1.5  # Number of standard deviations above mean for density
result = coverage(mol, obs_grid, options=opts)
```

### Working with Results

```python
result = rscc(mol, obs_grid, resolution, calc_grid=calc_grid)

# Overall score
print(f"RSCC: {result.overall:.3f}")

# Per-residue scores
for residue, score in result.by_residue.items():
    print(f"  {residue.chain}/{residue.name} {residue.number}: {score:.3f}")

# Per-atom scores
for atom_idx, score in result.by_atom.items():
    print(f"  Atom {atom_idx}: {score:.3f}")

# Formatted summary
print(result)
```

### Atom Filtering

All scoring functions accept an optional `mask` parameter to restrict the overall score to a subset of atoms:

```python
from openeye import oechem

# Score only chain A
mask = oechem.OEHasChainID("A")
result = rscc(mol, obs_grid, resolution, mask=mask, calc_grid=calc_grid)
```

### Grid Operations

```python
from maptitude import (
    scale_map, combine_maps, diff_to_calc, wrap_and_pad_grid,
    interpolate_density, get_grid_params, MapOp,
)

# Scale a map in place
scale_map(grid, 2.0)

# Combine two maps
summed = combine_maps(grid_a, grid_b, MapOp.ADD)
diff = combine_maps(grid_a, grid_b, MapOp.SUBTRACT)

# Derive calculated density from observed and difference maps
# calc = obs - 2 * diff
calc_grid = diff_to_calc(obs_grid, diff_grid)

# Handle CCP4 unit-cell maps where coordinates extend beyond the cell
padded = wrap_and_pad_grid(grid, mol, cell_a, cell_b, cell_c, padding=3.0)
if padded is not None:
    grid = padded  # Use the padded grid

# Sample density at a point
value = interpolate_density(grid, x, y, z)

# Inspect grid geometry
params = get_grid_params(grid)
print(f"Dimensions: {params.x_dim} x {params.y_dim} x {params.z_dim}")
print(f"Spacing: {params.spacing} A")
```

### Crystallographic Types

```python
from maptitude import UnitCell, parse_symop, parse_symops

# Unit cell
cell = UnitCell(50.0, 60.0, 70.0, 90.0, 90.0, 90.0)
print(f"Volume: {cell.Volume():.1f} A^3")

# Coordinate conversion
frac = cell.CartesianToFractional(25.0, 30.0, 35.0)
cart = cell.FractionalToCartesian(frac[0], frac[1], frac[2])

# Symmetry operators
op = parse_symop("x,y,z")
result = op.Apply(0.25, 0.5, 0.75)

ops = parse_symops("x,y,z\n-x,y+1/2,-z+1/2")
```

### Scattering Factors

```python
from maptitude import get_scattering_factors

# Look up Cromer-Mann coefficients for carbon (Z=6)
coeffs = get_scattering_factors(6, charge=0)
f0 = coeffs.Evaluate(0.0)  # Scattering factor at sin(theta)/lambda = 0
```

## Configuration Reference

### RsccOptions / RsrOptions

| Property              | Type         | Default                            | Description                                        |
|-----------------------|--------------|------------------------------------|----------------------------------------------------|
| `atom_radius_method`  | `AtomRadius` | `BINNED` (RSCC) / `ADAPTIVE` (RSR) | How atom scoring radii are determined              |
| `fixed_atom_radius`   | `float`      | `1.5`                              | Radius in Angstroms (when method is `FIXED`)       |
| `atom_radius_scaling` | `float`      | `1.0`                              | Multiplier for vdW radii (when method is `SCALED`) |

**AtomRadius methods:**

| Value                 | Description                                      |
|-----------------------|--------------------------------------------------|
| `AtomRadius.FIXED`    | Same radius for all atoms                        |
| `AtomRadius.SCALED`   | Atom vdW radius multiplied by a scaling factor   |
| `AtomRadius.BINNED`   | Resolution-dependent radius bins                 |
| `AtomRadius.ADAPTIVE` | B-factor and resolution dependent (Tickle, 2012) |

### QScoreOptions

| Property          | Type             | Default | Description                                      |
|-------------------|------------------|---------|--------------------------------------------------|
| `sigma`           | `float`          | `0.6`   | Gaussian reference width (Angstroms)             |
| `radial_step`     | `float`          | `0.5`   | Step between radial shells (Angstroms)           |
| `max_radius`      | `float`          | `2.0`   | Maximum sampling radius (Angstroms)              |
| `num_points`      | `int`            | `8`     | Sample points per radial shell                   |
| `normalize_map`   | `bool`           | `True`  | Normalize the map before scoring                 |
| `isolate_points`  | `bool`           | `True`  | Exclude shell points near neighboring atoms      |
| `radial_sampling` | `RadialSampling` | `FIXED` | Radial sampling strategy (`FIXED` or `ADAPTIVE`) |

### CoverageOptions

| Property | Type    | Default | Description                                              |
|----------|---------|---------|----------------------------------------------------------|
| `sigma`  | `float` | `1.0`   | Standard deviations above mean for the density threshold |

## References

- Pintilie, G. et al. (2020). "Measurement of atom resolvability in cryo-EM maps with Q-scores." *Nature Methods*, 17,
  328--334.
- Tickle, I. J. (2012). "Statistical quality indicators for electron-density maps." *Acta Crystallographica Section D*,
  68, 454--467.

## License

MIT License. See [LICENSE](LICENSE) for details.
