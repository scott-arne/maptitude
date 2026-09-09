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
- **CCP4/MRC map I/O** that preserves the MRC2000 `ORIGIN` record and the symmetry block, both of which the
  OpenEye grid readers drop
- **OpenMP parallelism** for structure factor accumulation

## Requirements

- Python 3.11 or later
- OpenEye Toolkits 2026.1 or later (with a valid license)

## Installation

Install from a wheel:

```bash
pip install maptitude
```

Building from source requires CMake 3.21+, SWIG 4.2+, and the OpenEye C++ Toolkits.
The FFT is PocketFFT, vendored in `third_party/pocketfft`; nothing else needs installing.
SWIG 4.2 is what the default stable-ABI build needs; 4.0 is enough with
`-DMAPTITUDE_USE_STABLE_ABI=OFF`.

## Quick Start

```python
from openeye import oechem
from maptitude import UnitCell, fc_density, get_unit_cell, read_map, rscc

# Load structure and observed density map. read_map places the grid on the
# origin the file declares; oegrid.OEReadGrid ignores that record.
mol = oechem.OEGraphMol()
oechem.OEReadMolecule(oechem.oemolistream("model.pdb"), mol)

obs_grid, symops = read_map("2fofc.ccp4")

# Calculate the model density to score against. get_unit_cell reports the
# grid's edges; the angles are yours to supply.
edges = get_unit_cell(obs_grid)
cell = UnitCell(edges.a, edges.b, edges.c, 90.0, 90.0, 90.0)
calc_grid = fc_density(mol, obs_grid, 2.0, cell)

# Score the fit
result = rscc(mol, obs_grid, 2.0, calc_grid=calc_grid)
print(f"Overall RSCC: {result.overall:.3f}")
```

> **Orthorhombic cells only.** `DensityCalculator` and `fc_density` currently
> support orthorhombic unit cells (all angles 90 degrees). Monoclinic and
> triclinic cells raise `CellError`. General lattice support is planned.

## Usage

### Scoring Metrics

All scoring functions return a `DensityScoreResult` with an overall score and per-residue and per-atom breakdowns.

#### RSCC (Real-Space Correlation Coefficient)

```python
from maptitude import rscc, RsccOptions, AtomRadius

result = rscc(mol, obs_grid, resolution, calc_grid=calc_grid)

# With custom options
opts = RsccOptions()
opts.atom_radius_method = AtomRadius.SCALED
opts.atom_radius_scaling = 1.5
result = rscc(mol, obs_grid, resolution, calc_grid=calc_grid, options=opts)
```

`rscc` supports `FIXED`, `SCALED`, and `BINNED`. `AtomRadius.ADAPTIVE` raises
`RuntimeError`; it is available on `rsr` only.

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

### Map Files

`read_map` and `write_map` handle CCP4 and MRC files directly. They exist because
`OEReadGrid` returns the payload, the cell and the space group but never consults
the MRC2000 `ORIGIN` record and does not expose the symmetry block, so a map read
through it can land in the wrong place and loses its symmetry.

```python
from maptitude import read_map, write_map, parse_symops, OriginSource

# Returns a MapFile, which unpacks as (grid, symops).
grid, symops = read_map("2fofc.ccp4")

# The symmetry text is in the newline-delimited form parse_symops accepts.
ops = parse_symops(symops)

# Write it back, restoring the ORIGIN and symmetry records.
write_map("out.ccp4", grid, symops)
```

`write_map` writes and verifies a temporary sibling of the destination and renames
it only once the file reads back as the grid it came from, so a grid that cannot be
written faithfully raises with the destination untouched, rather than leaving a
wrong map on disk or destroying a good one. A relative path is resolved against the
working directory once, on entry. The default space group is P1.

When a file sets both `ORIGIN` and `NxSTART` to nonzero values they can disagree
about where the map belongs. `OriginSource` picks the winner:

| Value                          | Places the map using          |
|--------------------------------|-------------------------------|
| `OriginSource.ORIGIN_RECORD`   | The MRC2000 `ORIGIN` record (default) |
| `OriginSource.NXSTART`         | The `NxSTART` node offsets    |

```python
grid, symops = read_map("map.ccp4", OriginSource.NXSTART)
```

Two limits are worth knowing before you rely on a round trip. The file's text
labels are not preserved: `write_map` splices back the symmetry block, `NSYMBT`
and `ORIGIN`, and the label block is not among them. And a box that
`wrap_and_pad_grid` actually padded is refused, because its declared cell equals
its full sampled extent; the refusal is a `GridError` naming the dimensions, and
the destination is left untouched.

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

# Handle CCP4 unit-cell maps where coordinates extend beyond the cell.
# Always returns a grid: the original one when no padding was needed.
grid = wrap_and_pad_grid(grid, mol, cell_a, cell_b, cell_c, padding=3.0)

# Sample density at a point
value = interpolate_density(grid, x, y, z)

# Inspect grid geometry
params = get_grid_params(grid)
print(f"Dimensions: {params.x_dim} x {params.y_dim} x {params.z_dim}")
print(f"Spacing: {params.x_spacing} x {params.y_spacing} x {params.z_spacing} A")
print(f"First node: {params.x_origin}, {params.y_origin}, {params.z_origin}")
```

#### Sampling

`interpolate_density` samples the node span, `[origin, origin + (n - 1) * spacing]`
on each axis, and returns `default_value` outside it. The periodic variants treat
the grid as one period of a map that tiles space, so a point in an axis's final
interval blends across the period boundary and no point is outside; there
`default_value` only covers a non-finite coordinate. They need the cell to be
commensurate with the spacing and raise `CellError` when it is not:

```python
from maptitude import interpolate_density_periodic, get_unit_cell

cell = get_unit_cell(grid)
value = interpolate_density_periodic(grid, x, y, z, cell.a, cell.b, cell.c)
```

Both have batch forms that take a flat `[x0, y0, z0, x1, y1, z1, ...]` sequence
and the point count, and return a tuple of densities:

```python
from maptitude import interpolate_density_batch, interpolate_density_periodic_batch

points = [x0, y0, z0, x1, y1, z1]
values = interpolate_density_batch(grid, points, 2)
values = interpolate_density_periodic_batch(grid, points, 2, cell.a, cell.b, cell.c)
```

`get_atom_grid_points` returns the grid nodes within a radius of a point, which is
what the scoring metrics sample over:

```python
from maptitude import get_atom_grid_points

nodes = get_atom_grid_points(grid, x, y, z, radius=1.5)
```

#### Geometry checks

```python
from maptitude import grid_contains, same_grid_geometry

# Is a Cartesian point inside the sampled box?
grid_contains(get_grid_params(grid), x, y, z)

# Do two grids sample the same box the same way? Required before combine_maps.
same_grid_geometry(grid_a, grid_b)          # default tolerance 1e-6
same_grid_geometry(grid_a, grid_b, tol=1e-4)
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
coeffs = get_scattering_factors(6, formal_charge=0)
f0 = coeffs.Evaluate(0.0)  # Scattering factor at sin(theta)/lambda = 0
```

`get_scattering_factor_table` returns the whole table as `(entries, count)`, where
each entry carries `atomic_number`, `formal_charge`, and its `coeffs`:

```python
from maptitude import get_scattering_factor_table

entries, count = get_scattering_factor_table()
for entry in entries:
    print(entry.atomic_number, entry.formal_charge, entry.coeffs.Evaluate(0.0))
```

## Configuration Reference

### RsccOptions / RsrOptions

| Property              | Type         | Default                            | Description                                        |
|-----------------------|--------------|------------------------------------|----------------------------------------------------|
| `atom_radius_method`  | `AtomRadius` | `BINNED` (RSCC) / `ADAPTIVE` (RSR) | How atom scoring radii are determined              |
| `fixed_atom_radius`   | `float`      | `1.5`                              | Radius in Angstroms (when method is `FIXED`)       |
| `atom_radius_scaling` | `float`      | `1.0`                              | Multiplier for vdW radii (when method is `SCALED`) |

**AtomRadius methods:**

| Value                 | Description                                      | Supported by |
|-----------------------|--------------------------------------------------|--------------|
| `AtomRadius.FIXED`    | Same radius for all atoms                        | `rscc`, `rsr` |
| `AtomRadius.SCALED`   | Atom vdW radius multiplied by a scaling factor   | `rscc`, `rsr` |
| `AtomRadius.BINNED`   | Resolution-dependent radius bins                 | `rscc`, `rsr` |
| `AtomRadius.ADAPTIVE` | B-factor and resolution dependent (Tickle, 2012) | `rsr` only   |

`ADAPTIVE` is the only value that is not shared. `rscc` has no adaptive radius
model, so passing it — as the enumerator or as the string `"adaptive"` — raises
rather than silently falling back to `BINNED`.

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

## Exceptions

All *domain* exceptions raised by maptitude derive from `MaptitudeError`, so a single
`except maptitude.MaptitudeError` catches every failure that describes your input. Four builtin
types are also raised, for failures in how the call was made or in the memory behind it rather than
in the structure or map it was given -- they are listed in the second table, and
`except maptitude.MaptitudeError` does not catch them.

| Exception        | Raised when                                                          |
|------------------|----------------------------------------------------------------------|
| `MaptitudeError` | Base class. Never raised directly.                                   |
| `StructureError` | The molecule is unsuitable: no heavy atoms, missing coordinates.     |
| `GridError`      | Grid geometry is wrong, a required grid is missing, a numeric argument is outside its usable range, or the copy that returns a grid to Python did not preserve the source geometry. |
| `SymOpError`     | A symmetry-operator string cannot be parsed.                         |
| `CellError`      | Unit-cell parameters are invalid or describe an unsupported lattice, including a map file whose sampling is not axis-aligned. |

| Exception      | Raised when                                                                       |
|----------------|------------------------------------------------------------------------------------|
| `RuntimeError` | An option value is out of range -- every option setter validates in C++ and its `std::invalid_argument` surfaces here. Also when the path that returns a grid to Python fails for a reason that is not the source geometry: an OpenEye wrapper that does not carry the pointer the copy needs, or a C++ exception other than `std::bad_alloc` raised while copying. No argument value selects those paths, but a bare `RuntimeError` is not a `MaptitudeError`, so `except maptitude.MaptitudeError` does not catch them. |
| `ValueError`   | A string argument the Python wrappers resolve against a table names nothing: `atom_radius="vdw"`, or `atom_radius="adaptive"` to `rscc`. Also an `int` outside a scoped enum's declared range that fits in a C `long`, which the `MapOp` and `OriginSource` typemaps check: `combine_maps(lhs, rhs, 47)` and `read_map(path, 47)`. One too large to fit raises `OverflowError` instead. |
| `TypeError`    | An argument has the wrong type: a `mask` that is not an OpenEye atom predicate, an `options` object of the wrong class, a `symops` that is not a string or a sequence of `SymOp`. |
| `MemoryError`  | The copy that returns a grid to Python could not be allocated. |

Those four are what maptitude raises itself. SWIG's argument conversion runs in front of the
library and raises builtins of its own, so an integer too large for the C++ parameter it feeds
arrives as an `OverflowError` from the conversion rather than a `RuntimeError` from the setter:
`num_points = 2**40` is one. One rejection appears in two of the table's rows: `rscc` refuses the
adaptive radius model as a `ValueError` when it is spelled `atom_radius="adaptive"` and as a
`RuntimeError` when it is spelled `atom_radius=AtomRadius.ADAPTIVE` or carried on an `RsccOptions`.
The split is deliberate -- the string is resolved against a three-entry table in Python, where an
unmatched name is a `ValueError` like any other, while the enum value names a real method that this
metric does not implement and is refused in C++ along with every other option-value rejection --
but it means catching both types is the only way to catch the rejection whatever spelling reaches
you.

```python
import maptitude

try:
    result = maptitude.rscc(mol, obs_grid, resolution=2.0, calc_grid=calc_grid)
except maptitude.GridError as exc:
    print(f"grid problem: {exc}")
except maptitude.MaptitudeError as exc:
    print(f"maptitude failed: {exc}")
```

### Rejected input

These are the input classes the library refuses rather than scoring. Where the same input was
accepted before this release, it variously returned a plausible wrong value or a NaN, hung,
exhausted memory, divided by zero, or read uninitialized memory; that list of mechanisms is
illustrative, not exhaustive. Some rows widen a rejection that already existed rather than adding
one: at `v0.2.4` a non-positive `resolution` and a metric on a molecule with no heavy atoms both
raised already.

| Input                                                                | Result                        |
|----------------------------------------------------------------------|-------------------------------|
| A non-finite or non-positive `resolution`, at any entry point         | `GridError`                   |
| A `resolution` and unit cell needing a Miller-index box over 2e8 points | `GridError`                 |
| A node interval at or above twice its own axis's cell edge            | `GridError`                   |
| `n_scale_shells` outside `[1, 1000]`                                  | `GridError`                   |
| A non-orthorhombic unit cell                                          | `CellError`                   |
| A zero, negative, or non-finite cell edge to `wrap_and_pad_grid`      | `CellError`                   |
| A cell edge that is neither `n` nor `n - 1` node intervals, or is further from that count than the float-coordinate allowance, to a periodic entry point or `wrap_and_pad_grid` | `CellError` |
| A molecule with no heavy atoms                                        | `StructureError`              |
| A Q-score radial sweep that cannot terminate or produce a shell       | `GridError`                   |
| An `AtomRadius` or `RadialSampling` value the enum does not declare   | `RuntimeError`                |
| `AtomRadius.ADAPTIVE` to `rscc`                                       | `RuntimeError`                |
| `atom_radius="adaptive"` to `rscc`, or any unrecognised radius name   | `ValueError`                  |
| A non-finite `CoverageOptions.sigma`                                  | `RuntimeError`                |
| A non-finite, zero, or negative `QScoreOptions.sigma`                 | `RuntimeError`                |

`CoverageOptions.sigma` and `QScoreOptions.sigma` are validated differently on purpose. Q-score's
sigma is a Gaussian width and has to be positive. Coverage's is a multiplier in the threshold
`mean + sigma * stddev`, where zero means "threshold at the mean" and a negative value means
"threshold below the mean" -- both meaningful requests. Only NaN and the infinities are refused
there, because the sigma enters that threshold directly and none of the three leaves one worth
comparing a density against. NaN makes the threshold NaN, so every `rho >= threshold` is false and
coverage returns a plausible `0.0`. On a map with nonzero spread, `+inf` makes the threshold `+inf`
and returns `0.0` as well, while `-inf` makes it `-inf` and returns a perfect `1.0`; on a flat map
`sigma * stddev` is NaN for either infinity, so both return `0.0`.

Under `RadialSampling.ADAPTIVE` the maximum sampling radius comes from the atom rather than from
the options, so a sweep that fails on one atom's radius scores that atom `NaN` and leaves the rest
of the molecule scored. A sweep no atom in the molecule can satisfy is a property of the resolution
and the grid spacing instead, and raises `GridError`.

## References

- Pintilie, G. et al. (2020). "Measurement of atom resolvability in cryo-EM maps with Q-scores." *Nature Methods*, 17,
  328--334.
- Tickle, I. J. (2012). "Statistical quality indicators for electron-density maps." *Acta Crystallographica Section D*,
  68, 454--467.

## License

MIT License. See [LICENSE](LICENSE) for details.

The wheels bundle one third-party component: PocketFFT, vendored unmodified in
[`third_party/pocketfft`](third_party/pocketfft) under the BSD-3-Clause license.
Its license text ships alongside maptitude's own in the wheel metadata.
