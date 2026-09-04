# Changelog

All notable changes to this project are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project is pre-1.0: breaking changes may land in a minor release.

## [0.5.0]

The map I/O release. Maptitude now reads and writes CCP4/MRC files through its
own `read_map` and `write_map`, which preserve the two header records the
OpenEye grid I/O discards: the MRC2000 `ORIGIN` placement and the symmetry
block. The Python loaders that stood in for this in the benchmark suite are
gone.

### Added

- `read_map(path, tiebreak=OriginSource.ORIGIN_RECORD)`, returning a `MapFile`
  of `(grid, symops)`. The grid arrives positioned at the origin the file
  declares, and the symmetry text arrives in the newline-delimited form
  `parse_symops` accepts. Available from C++ and Python.
- `write_map(path, grid, symops="")`. It writes and verifies a temporary
  sibling of the destination and renames it only once the file reads back as
  the grid it came from, so a grid that cannot be written faithfully raises
  with the destination untouched. Available from C++ and Python.
- `OriginSource`, the tiebreak selecting which record places a map when both
  `ORIGIN` and `NxSTART` are nonzero. `ORIGIN_RECORD` is the default.
- `compare_written_map` and `MAP_VERIFY_TOL` in `MapIO.h`, the predicate behind
  `write_map`'s verification and its relative tolerance; `compare_placement_records`
  and `MAP_PLACEMENT_SLACK`, the predicate behind `write_map`'s placement
  verification and its fractional slack. Public so a test can drive the
  comparison directly rather than having to manufacture a file that fails on
  each quantity. C++ only.
- Both `read_map` and `write_map` reject a path containing an embedded NUL, so a
  caller cannot write through a truncated path. Pinned at
  `RaisesOnAPathWithAnEmbeddedNul` and `RefusesADestinationWithAnEmbeddedNul` in
  the C++ suite, `test_write_map_refuses_a_path_with_an_embedded_nul` in Python.
- `MapFile` is importable from the `maptitude` package.

### Changed

- **EM density moves.** `read_map` places a map where its `ORIGIN` record says
  it belongs. On the one EM asset in the test suite that is 220 A from where the
  bare `OEReadGrid` puts it — the full `ORIGIN` vector, since that file's
  `NxSTART` is zero — so every score reads density from a different place than
  it did before.
- **A symmetry record that does not parse now raises `SymOpError`** where
  `OEReadGrid` did not read the block at all. A block whose framing is broken —
  `NSYMBT` negative or not a multiple of 80, or a file that ends before the
  block does — raises `GridError`.
- **`write_map` defaults an absent space group to P1,** so an EM map written
  from a grid with no space group lands with `ISPG` 1 rather than 0. An unset
  space group is what makes `OEWriteGrid` double the cell and regrid the
  payload, so the alternative is a corrupted file.
- **The scoped enums are checked at the Python boundary.** `combine_maps`'s `op` and
  `read_map`'s `tiebreak` previously accepted any integer and cast it, so an
  out-of-range op reached the C++ `switch`, matched no case, and returned an all-zero
  grid with no error. Now an `int` outside the declared range that fits in a C `long`
  raises `ValueError`; one too large to fit raises `OverflowError`; a non-`int` raises
  `TypeError`, and so does a `bool`, which subclasses `int` but names no enumerator —
  `combine_maps(lhs, rhs, True)` used to perform `SUBTRACT` silently.
- **The grid return path raises `GridError` or `CellError` where it raised
  `RuntimeError`.** The typemap that hands a C++ grid back to Python covers four
  pre-existing functions, not only the new `read_map`:
  `DensityCalculator.Calculate`, `combine_maps`, `diff_to_calc` and
  `wrap_and_pad_grid`. `GridError` and `CellError` derive from `MaptitudeError`,
  which derives from `Exception` and **not** from `RuntimeError`, so a handler
  catching `RuntimeError` around any of those four no longer catches this
  failure.

### Removed

- `benchmarks/helpers.py`'s `load_mrc_grid`, `load_ccp4_grid` and
  `extract_ccp4_symops`, replaced by `read_map`. `load_ccp4_grid` returned the
  three cell edges as a tuple; the skew carrier holds them, so call
  `get_unit_cell(grid)` instead.
- The `sys.path` inversion in `tests/python/conftest.py` that put `benchmarks/`
  on the test path. The only test that used it imported exactly the two loaders
  now retired.

### Known limitations

- `write_map` refuses a box that `wrap_and_pad_grid` actually padded — when no
  padding is needed it hands back the caller's own grid, which writes normally.
  The padded box declares a cell equal to its full sampled extent, which forces
  the written `NX` to equal `NC` and inflates the node count by one per axis on
  the way back in. The refusal is a `GridError` naming the dimensions it wrote
  against the ones it expected; the destination is left untouched. Lifting it
  means rewriting the declared cell, which would stop the written box comparing
  same-geometry with the box in memory.

## [0.4.0]

The anisotropic-carrier release. Every density grid maptitude accepts, returns or
stores is now an `OESkewGrid`, which carries a unit cell and so a separate node
interval per axis, in place of `OEScalarGrid` and its one scalar spacing. On an
isotropic map the two describe the same lattice; on an anisotropic one the scalar
path applied a single spacing to all three axes, resampling values onto a lattice
the experiment never produced. This is a breaking release across C++ signatures,
the SWIG layer and the Python API, with no shim, no dual code path and no
deprecation window. Downstream code that subscripts a grid, constructs one from
an extents box, or reads a spacing off one is affected; **The Python break**
below is written for that reader.

### Added

- `get_unit_cell`, `grid_contains` and `same_grid_geometry` — the unit-cell
  accessor, the in-grid predicate and the geometry comparison the per-axis
  carrier needs. All three are reachable from Python.
- `grid_node_origin`, `grid_fractional_index` and `interpolate_density_at`, C++
  only. The first two write their results through out-parameters;
  `interpolate_density_at` is deliberately not exposed through SWIG.
- `require_commensurate_cell` in `Grid.h`, the check the periodic entry points
  make on the cell edges they are given, as a public function. It accepts a cell
  edge spanning a whole number of node intervals — `n` of them or `n - 1`,
  according to whether the grid's closing plane repeats its first node — rejects
  anything else with `CellError`, and returns the accepted count per axis as a
  `CellPeriods`. A grid sampling more than one cell is rejected too: the counts
  are the two a grid whose nodes span a single period can have, not every count
  that could divide the node span. A caller about to sample a grid periodically
  can reject an unusable cell before doing any other work, which is how
  `wrap_and_pad_grid` uses it ahead of the centroid shift;
  `interpolate_density_periodic_batch` keeps the returned periods so the check
  runs once per batch rather than once per point. C++ only; a Python caller
  gets the behaviour through the functions that call it.
- `PAD_INTERVAL_COUNT_TOL` in `GridOps.h`, the relative tolerance
  `wrap_and_pad_grid` uses to recognise a whole number of node intervals in the
  extent it has to cover. It was a constant inside the function; it is public
  because it bounds how far short of the requested extent the padded grid may
  fall, and that bound is part of the function's contract.
- **A configure-time warning when the installed SWIG differs from the one that
  generated the checked-in proxy.** `python/maptitude/maptitude.py` is tracked,
  and the build's post-build copy rewrites it with the output of whatever SWIG
  is installed, so a contributor on another SWIG restaged a wholly regenerated
  proxy alongside their real change with nothing to say so — observed with
  4.4.1 overwriting a 4.5.0 proxy. The check reads the version from the tracked
  file's own banner, so there is no second constant to keep in sync, and it
  warns rather than failing, so a mismatch reports itself without blocking the
  build.

### Changed

- **`OESkewGrid` replaces `OEScalarGrid` as the density carrier** everywhere it
  appeared: C++ signatures, the SWIG typemaps, and the Python API. There is no
  overload taking the old type and no conversion layer, so every caller moves at
  once. The `OEScalarGrid` SWIG machinery is deleted.
- **`GridParams` carries `x_spacing`, `y_spacing` and `z_spacing`; the scalar
  `spacing` field is gone.** Code that read `gp.spacing` reads three fields now,
  and code that assumed one spacing has to decide which axis it meant.
- **`vector_to_grid` raises `GridError` when the vector's length is not the
  grid's element count**, where it previously copied the shorter of the two and
  reported nothing. Too short, too long and empty all raise
  (`tests/cpp/test_grid_ops.cpp:351`).
- **`MAX_FFT_GRID_POINTS` (`2e8`) bounds the calculated-density FFT sampling
  grid.** A request whose per-axis counts multiply out past that now raises
  instead of attempting the allocation.
- **A spacing that spans the cell edge in one node interval is still accepted,
  but the accepting case is now stated as such** — the fixture that pinned it is
  `DensityCalculatorValidationTest.StillAcceptsOneNodeIntervalSpanningTheCellEdge`
  (`tests/cpp/test_input_validation.cpp:725`). This is a narrowing of what the
  suite claims, not of what the library accepts.
- **`combine_maps` rejects grids whose spacing differs by less than the old
  tolerance.** The geometry gate is `same_grid_geometry`, which compares the
  derived per-axis geometry rather than a single spacing within a loose
  tolerance, so a pair that used to combine can now raise
  (`GridOpsTest.CombineRejectsGridsWhoseSpacingDiffersBelowTheOldTolerance`).
- **The test suite gained an independent third-party reader as a control.**
  `tests/python/test_gemmi_control.py` decodes the same CCP4/MRC assets with
  gemmi and compares them to what maptitude reads through OpenEye, which is the
  only check in the suite that can see a format misreading maptitude and OpenEye
  share. gemmi is a **dev/test-only** dependency: it is declared in the `dev`
  extra and the wheel's runtime dependencies are unchanged. CI installs it in the
  three `Test wheel` steps of `.github/workflows/build-wheels.yml`, which
  previously installed only `pytest`, `numpy` and the built wheel and would have
  failed collection on a clean runner. That workflow runs on `v*` tags, so the
  gap was latent until a release tag was pushed. On the four assets it covers the
  control rules out a shared format misreading insofar as it shows in the dims,
  the unit cell and the node values it walks; an off-by-one at the closing plane;
  a node origin misplaced by more than 1e-3 A on `390_emd_30342_A_z4.mrc`, the
  only asset its origin test reads; and an x/y transposition of the grid *dims*
  on that same asset, the only one whose x and y dims differ. It does **not**
  rule out an x/y transposition of the grid *spacings*: on three of the four
  assets the x and y spacings are exactly equal, and on the fourth they differ by
  1.0e-08 relative, a hundredfold under the `rel=1e-6` the spacing test asserts
  at, so such a swap passes the module. Nor does it catch a misplacement both
  readers share; the MRC `ORIGIN` record is that case and is recorded as such.
- **Periodic interpolation is now genuinely periodic across the cell boundary.**
  `interpolate_density_periodic`, `interpolate_density_periodic_batch` and the
  padded grid built by `wrap_and_pad_grid` wrap on the period the caller's cell
  names rather than on the node count: for a cell edge spanning `p` node
  intervals, node `p - 1`'s upper neighbour is node 0, so a point in that last
  interval blends the two. Which period applies is a real choice, not a
  formality. `1d26_2fofc.ccp4`'s own cell is `p = n - 1 = 48`, because the
  reader appends a closing plane. The two periods part company over `x` between
  48 and 49 node intervals above the first node — past the period on the
  `p = 48` reading, the period's own last interval on the `p = 49` one. Holding
  `y` and `z` 3.5 spacings above their first nodes and taking `x` a quarter, a
  half, three quarters and 0.99 of the way through that interval, `p = 48` gives
  `-0.0017`, `0.0213`, `0.0442` and `0.0662` while `p = 49` holds flat at
  `-0.0246`, the plateau produced by blending the duplicated closing plane
  against its source. The interval below it separates nothing: node 48 repeats
  node 0 to the bit at all 49x25 nodes of that face, so `p = 48`, `p = 49` and
  the non-periodic entry point all return `0.0423`, `0.0200`, `-0.0023` and
  `-0.0237` there.
  Previously the wrap produced coordinates one node interval wider than the
  domain the interpolator accepts, and every point in that final interval came
  back as `default_value`: on a 10-node, 1.0 A, cell-10 grid, `x` anywhere in
  `[9, 10)` returned the default. `wrap_and_pad_grid` baked that gap into the
  padded map as zero density — 440 of 1331 voxels on a uniformly filled test
  grid. There is no longer an outside to fall into: `default_value` is returned
  only when a point's fractional index is not finite, which covers a non-finite
  coordinate and also a finite one large enough that dividing it by the spacing
  overflows.
- **The periodic entry points and `wrap_and_pad_grid` raise `CellError` for a
  cell edge that is not a whole number of the grid's node intervals**, and that
  number has to be the node count or one less — `n * spacing` or
  `(n - 1) * spacing` per axis, according to whether the closing plane repeats
  the first node. Making the last node adjacent to the first only reproduces
  the crystal when the cell closes on nodes the grid holds; wrapping an
  incommensurate cell returned density from the wrong place with nothing to
  mark it as wrong. Callers passing a cell edge that is neither must now
  correct it. `wrap_and_pad_grid` makes the check before it shifts the
  molecule, so a rejected cell leaves the caller's coordinates where
  they were; an incommensurate cell that happened to need no padding previously
  translated the molecule by a vector that was not a lattice vector of the map
  and reported nothing. The `CellError` text names which of the two tests
  refused the edge, because they call for opposite responses: a wrong interval
  count means the cell does not belong to this grid, while a count that lands
  and a residual that does not means the grid's own coordinates no longer pin
  its spacing that finely. Callers matching on the wording must update.
- **A point on a grid's own corner interpolates instead of falling out of the
  grid.** The node span is derived from `ElementToSpatialCoord`, and its error
  scales with the magnitude of the floats OpenEye holds the geometry in: a
  five-node grid nominally starting at 12.3 A reports that node 1.9e-7 A away,
  where the same grid at the Cartesian origin is off by 1.5e-15 A. A caller
  querying the coordinate they built the grid from was told the point lay
  outside. The domain test and the commensurability test now share one tolerance
  scheme, a counted number of half-ulps of float times the largest magnitude the
  axis's geometry takes — the further endpoint, or the sampled extent
  `n * spacing`, whichever is greater. Discrimination is preserved, but it is no
  longer a constant: half a node interval is
  `spacing / (12 * 2^-24 * magnitude)` times the boundary tolerance, so the
  ratio falls as a grid's coordinates grow. It is about 5500x for a 256-node
  0.9 A map at the Cartesian origin and about 420x for a 20-node map of that
  spacing 3000 A out. The boundary tolerance would reach half a node interval
  only at a magnitude of about 1.4 million node intervals, which for a 0.9 A map
  is over 100 micrometres.
- **The commensurability allowance is capped at a quarter of a node interval.**
  Scaling it by the axis magnitude let it grow without bound, while the count
  test it feeds rounds the edge-to-spacing ratio and so already leaves the edge
  within half an interval of the count it picks. An allowance at half an interval
  was therefore satisfied by every edge that got that far, and the
  commensurability check stopped adding anything to the count test: measured, a
  12-node axis at 2.0 A spacing with origin 3e6 A accepted edges of 23.02 A and
  24.98 A alongside its true 24.0 A extent, and neither is a whole number of its
  node intervals. The allowance is now the smaller of the counted rounding budget and a
  quarter of the node interval, so it stops growing at an axis magnitude of
  `2^19` node intervals — about 47 micrometres for a 0.9 A map. A cell edge more
  than a quarter of an interval from the count it rounds to now raises
  `CellError` where such an edge on a distant grid was previously accepted.
  A quarter of an interval is the allowance's ceiling, and is its width only from
  that magnitude up; below it the counted budget governs and is narrower — near
  the Cartesian origin, far narrower, so a 10-node 1.0 A grid at origin 0 gets
  4.77e-06 A and refuses a 10.0001 A edge.
- **Far enough from the Cartesian origin, a grid can be refused the cell
  `GetUnitCell` reports for it.** Float node coordinates quantize, which
  perturbs the spacing the geometry derivation recovers, and the extent test
  scales that perturbation by `n / (n - 1)` — so a low node count amplifies it,
  and past a quarter interval the cap stops forgiving the difference. Acceptance
  is not monotone in the translation and the onset depends on the axis's shape,
  so there is no single threshold worth quoting; the earliest case found is a
  two-node axis at `1.34 * 2^20` node intervals, whose requested 1.3 A spacing
  quantized to 1.5 A. This is a consequence of the cap rather than a regression
  it introduced, and a deliberate trade: uncapped, the allowance reaches the
  half-interval the count test already guarantees at `2^20` node intervals, and
  the extent test discriminates nothing at all from there up. For a 0.9 A map
  `2^20` intervals is some 94 micrometres from the origin; the furthest node in
  any of this repository's five test maps sits 89 node intervals out.
  `GridOpsTest.InterpolateDensityPeriodicRefusesADistantGridsOwnCell` pins the
  refusal.
- **The grid `wrap_and_pad_grid` builds covers the extent it was sized for.**
  The node count truncated the interval count rather than rounding it up, so an
  atom extent that was not a whole number of node intervals produced a grid
  short by up to one interval per axis, leaving the outermost atoms outside the
  padded grid's own domain. The interval count is snapped to a whole number it
  is within `PAD_INTERVAL_COUNT_TOL` of before rounding up, so an extent that is
  an exact multiple does not buy a spurious node; the padded span may therefore
  fall short by up to that fraction of the requested extent, which is deliberate
  and is part of the function's documented contract.
- **`wrap_and_pad_grid` requires a finite, non-negative `padding`, and rejects an
  extent it cannot size a grid for.** The padded grid's node count is derived
  from the atom extent widened by the padding, and the interval count it
  produces becomes a grid dimension through a conversion that is undefined for a
  value outside the range that type represents. Neither term was bounded. The
  cell edge the padded grid is built with is that dimension times the source
  grid's node interval, narrowed to a `float` under the same rule; bounding the
  dimension does not bound that product. A NaN padding compared false against
  every grid face and returned `nullptr`, reporting that the atoms already fit;
  it now raises `GridError`. So does a negative padding, which shrank the box
  the atoms had to fit inside rather than widening it, and so does any atom
  extent and padding that between them need more node intervals than a grid
  dimension can hold, or that size a cell edge past what a `float` represents.
  Zero padding remains admissible. Like the cell edges, the padding's own
  finiteness and sign are checked before the centroid shift, so those two
  rejections leave the caller's molecule where it was. The interval-count and
  cell-edge limits are only known once the padded geometry has been derived,
  which is after the shift, so a molecule that needed a shift and is then
  rejected there is left translated by whole cell vectors -- to within one float
  store's rounding at the coarser of its old and new coordinates -- a position
  crystallographically equivalent to the one its own float coordinates fixed,
  not a corrupted one.
- **`wrap_and_pad_grid` rejects a centroid shift it cannot store, and rejects it
  before moving any atom.** The shift is computed in double and each shifted
  coordinate was narrowed to a `float` as that atom was written back, with
  nothing bounding the narrowing; a floating-point conversion is undefined for a
  value outside the destination's range, and on this arm64 host it saturates to
  an infinity. How far the molecule starts from the grid centre is not what did
  it -- the shift brings the centroid back to within half a cell edge of that
  centre from any starting distance at which a double still resolves the cell
  edge, which the reproduction below is well inside, so an atom ends up near the
  grid's centre, offset by its own displacement from the centroid. The infinity
  was written when those together passed the float maximum. In the case this was
  reproduced from, the grid's centre was itself a quarter of a cell edge short
  of that maximum, and what the caller then saw was the sizing guard below
  rejecting the infinite extent the write produced -- an error raised over
  coordinates the call had already replaced; and because that extent is measured
  over the heavy atoms alone, an overflow confined to the rest of the molecule
  had nothing later in the function looking at it. Every atom's shifted coordinate is now computed
  before any of them is written, and a component that is not a finite value a
  `float` can hold raises `GridError` with no atom modified, so a shift that
  cannot be stored no longer leaves the molecule part-way through one. A shift
  every atom's float coordinates can hold is applied exactly as before.
- **The advertised SWIG floor is 4.2, correcting a 4.0 the build already refused
  to honour.** `openeye_add_swig_module` runs `find_package(SWIG 4.2 REQUIRED)`
  when stable-ABI output is requested, and `MAPTITUDE_USE_STABLE_ABI` defaults
  to `ON` and is pinned `ON` for wheels, so the default build and every wheel
  build stop at configure time below 4.2 while the README and the
  `build-system` requires both advertised 4.0. Measured with a SWIG reporting
  4.1.0: `CMake Error ... required is at least "4.2"`, naming the module line
  that asked. Both declarations now say 4.2, and the README records that 4.0
  is still enough under `-DMAPTITUDE_USE_STABLE_ABI=OFF`, which is the version
  that configuration's module path asks for.

  The top-level `find_package(SWIG 4.0)` probe deliberately keeps its 4.0.
  Raising it to match was tried and measured worse: because that probe is not
  `REQUIRED` and only sets the `SWIG_FOUND` guard, asking it for 4.2 turned the
  error above into a configure that succeeds and silently omits the Python
  bindings. A permissive probe followed by the module's own `REQUIRED` check is
  what makes the failure loud.

  This corrects the advertised number to the one the build enforces; it is not a
  report that 4.2 was tested. The only SWIG this project is built and tested
  against is 4.5.0, which is also the current release.

### Removed

- **`OEGridSameGeometry` as the geometry comparison.** `same_grid_geometry`
  replaces it and is not a drop-in: where the old check returned `false` for a
  grid whose geometry cannot be derived, the new one raises. Measured, a grid
  built with a zero unit cell — whose node coordinates come back `nan` — gives
  `GridError: Grid element 0 (walking axis x) has a non-finite x coordinate
  (nan); the geometry cannot be derived`. `GridError` derives from
  `MaptitudeError`, which derives from `Exception` and **not** from
  `RuntimeError`, so a handler catching `RuntimeError` will not catch it.
- **The implicit `OEScalarGrid` argument conversion.** Passing an `OEScalarGrid`
  to a maptitude function raises at the boundary instead of being converted:
  measured, `TypeError: Expected OESkewGrid-derived object.`
- `grid_bounds`, a bounding-box helper added earlier in this unreleased cycle and
  withdrawn before release. It reported the scalar carrier's box — half a node
  interval outside the first and last nodes on each face — which is not the
  domain any maptitude function interpolates over. Use `grid_contains`, or
  `get_grid_params` and the node span.

The Python-surface removals that follow from the carrier switch are grouped in
**The Python break** below rather than repeated here.

### Fixed

- **`wrap_and_pad_grid` aimed the centroid shift half a node interval off centre
  on every axis.** The padded grid is now centred where the arithmetic intended.
  The value this moves is the shift asserted by
  `GridOpsTest.WrapAndPadGridShiftsCoordinates`, whose tolerance was exactly the
  size of the error.

- **`grid_contains` answered yes almost everywhere for a `GridParams` with a
  zero-node axis.** `GridParams` is a public struct callers fill in directly, so
  `get_grid_params`'s "at least two nodes on every axis" check does not stand on
  that path. The node count was subtracted as `unsigned`, so a `dim` of 0 wrapped
  `0u - 1u` to `2^32 - 1`: that put the axis's far endpoint there and widened the
  span tolerance to some 1536 A, giving a measured accepted range of
  `[-1536 A, 4.29e9 A]`. The subtraction is now done in `double`, and a zero-node
  axis at the Cartesian origin admits nothing.
  `GridContainsDegenerateParams.AZeroNodeAxisNearTheOriginAdmitsNothing` pins it.
  Grids from `get_grid_params` are unaffected — for any `dim` of 1 or more the
  two forms agree exactly. The uncapped boundary tolerance every axis carries is
  unchanged, so a zero-node axis far from the origin still admits a band.

### The Python break

The carrier switch removes five classes of method from the grid objects a Python
caller handles, because `OESkewGrid` does not have them. Each item below is the
call that stopped existing and the path that replaces it. Every absence in this
section was confirmed by attribute lookup against `OESkewGrid` on OpenEye
2026.1.0; every replacement was confirmed by running it.

**Value access.** `__getitem__`, `__setitem__`, `__iter__`, `GetValue`,
`SetValue`, `SetAll` and `CheckValues` are gone. Subscripting is the most common
downstream idiom this release breaks — `grid[i]` now raises `TypeError:
'OESkewGrid' object is not subscriptable`. Read through `GetValues()` and write
through `SetValues(OEFloatArray, int)`, against the x-fastest linearization
`el = iz * n_x * n_y + iy * n_x + ix`:

```python
# read: hoist GetValues() out of the loop -- it returns a list copy of the
# whole grid, so calling it per element copies the grid per element, and
# writing into the copy does not touch the grid.
values = grid.GetValues()
for iz in range(n_z):
    for iy in range(n_y):
        for ix in range(n_x):
            v = values[iz * n_x * n_y + iy * n_x + ix]

# write: fill one OEFloatArray, then one SetValues call.
arr = oechem.OEFloatArray(grid.GetSize())
for i in range(grid.GetSize()):
    arr[i] = compute(i)
assert grid.SetValues(arr, grid.GetSize())
```

**Containment.** `IsInGrid` is gone. `grid_contains` replaces it and is a
narrower predicate — see the in-grid domain item under **Why recorded values
moved** — so it is not a rename.

**Index and coordinate maps.** `GridIdxToElement`, `GridIdxToSpatialCoord`,
`SpatialCoordToElement`, `SpatialCoordToGridIdx`, `GetXIdx`, `GetXInc` and their
`GetY`/`GetZ` counterparts are gone, and so are the bare coordinate accessors
`GetX`, `GetY` and `GetZ`. Those three are named separately because they are the
inverse call, not a counterpart: `GetX(ix)` takes a node index and returns its
coordinate, where `GetXIdx(x)` takes a coordinate and returns a node index.
Derive the whole mapping from `get_grid_params(grid)`, whose per-axis origin, dim
and spacing give node positions directly — node `ix` sits at
`x_origin + ix * x_spacing` — and use the linearization above.

The coordinate-to-index direction needs its own recipe, and the binning is the
part that catches people out:

```python
import math

gp = maptitude.get_grid_params(grid)

def x_idx(x):
    # Half-spacing bins centred on nodes, matching the removed GetXIdx.
    raw = math.floor((x - gp.x_origin) / gp.x_spacing + 0.5)
    return min(max(raw, 0), gp.x_dim - 1)
```

`round()` is the wrong reach. On the geometry this was measured over — an
`OEScalarGrid` built from the extents box `[-1.5, 1.5]` at spacing 1.0, which
OpenEye expanded to four nodes at -1.5, -0.5, 0.5 and 1.5 — the clamped `floor`
above reproduced `GetXIdx` on all 20 coordinates probed, and on the 14 of those
`SpatialCoordToGridIdx` answered for it returned that same binning applied per
axis, `SpatialCoordToElement` its linearization. Substituting
`round((x - gp.x_origin) / gp.x_spacing)` for the `floor` disagreed on two of
the 20, at exactly `x = -1.0` and `x = 1.0`, because Python rounds half to even
where OpenEye takes the upper bin. It is a silent one-voxel shift, and it does
not announce itself by striking every half-way coordinate: on a 10-node 1.0 A
grid probed at all nine coordinates half way between adjacent nodes, `round`
missed the five whose lower node index is even and matched the four whose lower
index is odd.

The two calls differ out of range, and only one of the two behaviours has a
drop-in. `GetXIdx` clamped — on that grid it returned `0` for `x = -3.0` and `3`
for `x = 3.0` — which is what the `min`/`max` above reproduces.
`SpatialCoordToGridIdx` instead raised `IndexError: spatial coordinate out of
range` outside the old half-spacing box, accepting `-2.0` and rejecting `2.0`.
No maptitude predicate reproduces that domain: `grid_contains(gp, x, y, z)`
tests the node span rather than the box and is therefore narrower — on that grid
it is `True` at `-1.5` and `1.5` and `False` at `-1.75`, where
`SpatialCoordToGridIdx` still answered `(0, 2, 2)`. Code that needs the old
domain has to widen the node span by half a spacing on each face itself, and
close only the lower one. On three grids probed on all three axes, one of them
with its origin away from zero, the accepted interval per axis was
`[origin_i - s_i/2, origin_i + (n_i - 1) * s_i + s_i/2)`: the low face answered
and the high face raised, the largest accepted coordinate sitting one float step
below it.
`grid_fractional_index`, which computes the unclamped fractional index, is C++
only and is not on the Python surface.

**Box geometry.** The `GetXMin`/`GetXMax` family, `SetXDim`, `SetXMid` and
`IsXMidSet`, with their `Y` and `Z` counterparts, are gone. For reading, the node
span `[origin_i, origin_i + (n_i - 1) * spacing_i]` from `get_grid_params`
replaces the box, and it is not the same region: the old box lay half a node
interval outside the first and last nodes on each face.

That span does not migrate the setters. `OESkewGrid` has no per-axis mutator, so
code that changed one axis reads the other two back and passes all three:

```python
gp = maptitude.get_grid_params(g)

# SetXMid(mx) becomes: recompute the two midpoints you are keeping.
mid = [gp.x_origin + (gp.x_dim - 1) * gp.x_spacing / 2.0,
       gp.y_origin + (gp.y_dim - 1) * gp.y_spacing / 2.0,
       gp.z_origin + (gp.z_dim - 1) * gp.z_spacing / 2.0]
assert g.SetMid(mx, mid[1], mid[2])

# SetXDim(nx) becomes: pass the other two dims, then rebuild the cell.
assert g.SetDim(nx, gp.y_dim, gp.z_dim)
assert g.SetUnitCell(nx * gp.x_spacing, gp.y_dim * gp.y_spacing,
                     gp.z_dim * gp.z_spacing, 90.0, 90.0, 90.0,
                     nx, gp.y_dim, gp.z_dim)
```

That second `SetUnitCell` is not optional and nothing will remind you. `SetDim`
leaves the unit cell at its old edges while `get_grid_params` reports the new
geometry and stays self-consistent, so the grid still reads plausibly and fails
only where the cell is used. Measured on a 4x6x8 grid spaced
`(1.0, 0.5, 0.25)`, `SetDim(7, 6, 8)` alone left `get_unit_cell` reporting
`a = 4.0` against a sampled extent of 7.0, and `interpolate_density_periodic`
raised `CellError: ... got 4 A against 7 A (7 nodes at 1 A, every node
distinct) or 6 A (6 intervals, ...) ...`; reapplying `SetUnitCell` as above
cleared it and the call returned. `SetMid` on its own does
not need that second call: on the same grid it left all three dims and spacings
and the cell itself untouched, moving only the origin.

**Construction.** Both `OEScalarGrid` constructors are gone — the extents-box
form `OEScalarGrid(OEDoubleArray, spacing)` and the seven-argument
`OEScalarGrid(nx, ny, nz, mx, my, mz, spacing)`. One recipe replaces both:

```python
g = oegrid.OESkewGrid()
assert g.SetDim(nx, ny, nz)
assert g.SetUnitCell(nx * sx, ny * sy, nz * sz, 90.0, 90.0, 90.0, nx, ny, nz)
assert g.SetMid(mx, my, mz)
```

A caller who was passing an extents box `minmax` derives the recipe's arguments
from it. The dim needs care: `int(span / sp) + 1` truncates a quotient that binary
division has left just below an integer, and on the boxes measured below that cost
a whole boundary plane on 471 of 31 680 axes, with nothing raised to signal it.
Snap to the integer instead when the spacing divides the span:

```python
import math

def dim_from_span(span, sp):
    # The old constructor stepped exactly at span == k * sp, so a quotient that
    # lands a few ulps low must not truncate to k - 1.
    q = span / sp
    n = round(q)
    if n >= 0 and abs(n * sp - span) <= 1e-9 * max(1.0, abs(span)):
        return n + 1
    return math.floor(q) + 1

dim = [dim_from_span(minmax[i + 3] - minmax[i], sp) for i in range(3)]
mid = [(minmax[i] + minmax[i + 3]) / 2.0 for i in range(3)]
```

This is not exact everywhere. The old constructor's dim follows the floating-point
value of `minmax[i + 3] - minmax[i]` rather than the span written in the source, so
where the box sits reaches the dim through the rounding of that subtraction: at a
spacing of 0.05 the box `[0.0, 0.05]` gave dim 2 while `[2.5, 2.55]` gave dim 1,
because `2.55 - 2.5` is `0.04999999999999982`. Across 2 430 boxes at nine
placements, two that produced the same float span at the same spacing never
disagreed, so placement reaches the dim only through that rounding.

Over 10 560 boxes at twelve spacings and ten placements — 31 680 axes — the recipe
was exact on all 960 boxes whose three lower corners were `0.0`, and differed from
the old constructor on 888 axes. **Every one of those 888 was one node too many,
never one too few**: across that set the recipe never dropped a boundary plane, and
where it differed it left a redundant one — the opposite of the `int(span / sp) + 1`
failure above. A difference needs the snap's own condition to fire, and that
condition is a tolerance band rather than exact divisibility: it fires when the
span lands within `1e-9 * max(1, abs(span))` of the nearest integer multiple of
the spacing. Landing inside the band is not the same as sitting a whole number of
spacings from the opposite edge, and on a separate 9 600-box sweep the band is
what did the work — the predicate fired on all 251 differences there, while only
9 of those spans were an exact multiple. A span `1e-10` short of one spacing is
already inside the band, because the band is never narrower than `1e-9`: at
`sp = 1.0` the box `[0.0, 0.9999999999]` gives dim 1 from the old constructor and
2 from the recipe. Where a positive span misses the band the recipe reduces to
`int(span / sp) + 1` — over 3 000 boxes with arbitrary unrounded corners the
predicate never fired once, and both recipes agreed with the old constructor on
every axis. A caller who needs the old dims exactly should not recompute
them from the box at all: the old grid is still in hand during migration, and
`GetXDim()`, `GetYDim()` and `GetZDim()` read the true dims straight off it.

Rebuilding through the recipe, `get_grid_params` reported back exactly the dims
passed to `SetDim` on all 600 rebuilds measured, and the two carriers' midpoints
matched exactly on all 600 — `GetXMid()`, `GetYMid()` and `GetZMid()` returned
identical values on every one.

What costs precision is recomputing that midpoint instead of reading it.
Deriving it as `x_origin + (x_dim - 1) * x_spacing / 2` from `get_grid_params`
reproduced the old grid's value on only 271 of those same 600, drifting up to
7.7e-06. The spacing `get_grid_params` reports is re-derived from the cell
rather than the value handed to the recipe — on one 36-node axis at a nominal
0.35, `z_spacing` came back as `0.3500004359654018` — and `(dim - 1) / 2`
multiplies that gap straight into a derived midpoint; on that axis alone
17.5 × 4.4e-07 is 7.6e-06. Read the midpoint off the grid.

Node *coordinates* differ even when they are read rather than derived. Over the
428 rebuilds whose dims matched the old constructor's — grids whose midpoints
and dims were therefore identical — `ElementToSpatialCoord(0)` returned a
different origin on 178 of them, by up to 2.1e-06. Both carriers round
coordinates to 32-bit float. An assertion comparing rebuilt node coordinates
against old ones needs a tolerance, not `==`.

**Writing a constructed grid.** A grid built by that recipe needs two more lines
before `OEWriteGrid`, or the file it writes is wrong — and `OEWriteGrid` returns
`True` either way. Measured on seven recipe-built geometries, reading the origin
out of the raw CCP4 header rather than through a reader, every unfixed write was
wrong, in one of three ways: five of the seven re-expand the map by one node per
axis; one keeps the right node count but misplaces the origin; and the one
anisotropic grid is resampled onto a single uniform interval, losing its per-axis
geometry. That last is the worst of the three for this release's audience — a
21x20x11 grid spaced `(0.5, 0.5, 0.25)` was written as 41x39x11 nodes at a uniform
0.25 A, the finest of its three spacings, growing the stored map from 4 620 values
to 17 589.

```python
assert g.SetUnitCell(2 * nx * sx, 2 * ny * sy, 2 * nz * sz, 90.0, 90.0, 90.0,
                     2 * nx, 2 * ny, 2 * nz)   # doubled cell, unchanged spacing
assert g.SetSpaceGroup(1)                      # P1
```

The two lines restore the **dimensions** on all seven geometries, and the per-axis
spacings on the anisotropic one; on all seven they leave the grid in memory
unchanged — same dims, per-axis spacings, node origin and midpoint. They restore
the **origin** exactly on 17 of the 21 axes measured, and those 17 are exactly the
axes where `origin_i / spacing_i` is a whole number. CCP4 stores the origin as the
integer node index `NxSTART`, and OpenEye left the float `ORIGIN` field at zero in
all fourteen writes, so an axis whose node origin falls half way between two
multiples of its spacing cannot be written exactly: 20³ at 0.5 A centred on 0 has
an in-memory node origin of `-4.75` and reads back at `(-5.0, -5.0, -4.5)` with
both lines applied — half a node interval off, and not with the same sign on every
axis.

Where the origin lands *without* the two lines depends on the geometry, so the
figure below is given for the one this repository actually writes.
`tests/data/make_test_map.py` builds 21³ at 0.5 A centred on 0 — odd dimension,
origin-centred, the phase's only real write-path exemplar — and there the unfixed
file is exactly one node interval high on every axis: the node origin is `-5.0`,
the exact node index is `-10`, and the unfixed header carries `-9` where the fixed
header carries `-10`. That `+ 1` is specific to this exemplar and is not the
general rule; on four of the seven geometries the unfixed writer already emitted
the exact index and the origin came out right, and it was the node count that was
wrong instead.

This is an OpenEye write-path behaviour that the carrier switch exposes, not a
maptitude API change. A grid that came from `OEReadGrid` does not need either
line, so code that reads, scores and discards is unaffected; code that builds a
grid and writes it is.

**The breaks that are not vanished methods.** `GetSpacing()` and `SetSpacing()`
both survive the carrier switch, and neither is safe on an anisotropic map: one
reports a spacing that is not the map's, the other sets one that is not the one
asked for. `SetSpacing()` damages an isotropic map too, by a second route — it
sets the nodes correctly there and leaves the unit cell behind.

`GetSpacing()` still returns a single scalar, where there are now three spacings.
Measured on three constructed grids whose smallest node interval was on x, on y
and on z in turn, it returned the smallest of the three in every case — so it is
right on an isotropic map and **wrong and silent** on an anisotropic one.

`SetSpacing()` is the same error in the write direction, and it returns `True`.
It sets the *smallest* of the three intervals to the value asked for and scales
the other two by that same factor: over seven calls spanning five geometries the
result was the old spacings multiplied by `requested / min(old)` on every axis,
to float precision. A 4x6x8 grid spaced `(1.0, 0.5, 0.25)` given
`SetSpacing(0.5)` came back spaced `(2.0, 1.0, 0.5)` — dims and midpoint
unchanged, every interval doubled, the node extent doubled from
`(3.0, 2.5, 1.75)` to `(6.0, 5.0, 3.5)`.

The cell is the second route. Its edges were unchanged after the call in all
seven, so the five whose scale factor was not 1 ended up with a cell their nodes
no longer tile: the `(1.0, 0.5, 0.25)` grid above kept a `(4.0, 3.0, 2.0)`
cell against a sampled extent of `(8.0, 6.0, 4.0)`. That is the desync described
under **Box geometry**, with the same consequence —
`interpolate_density_periodic` raised `CellError` while `interpolate_density`
returned without complaint. Both isotropic grids measured took this half of the
break and not the other: on 4³ at 1.0 A, `SetSpacing(0.5)` moved the nodes the
way `OEScalarGrid.SetSpacing` did on the same geometry — dims and midpoint kept,
the interval halved — and still left a 4.0 A cell over 2.0 A of nodes.

`OESkewGrid` carries no per-axis setter of any kind, so there is nothing narrower
to reach for. Change a spacing by rebuilding the cell with the **Construction**
recipe instead: on the 4x6x8 grid above, `SetUnitCell`
with the wanted edges followed by `SetMid` gave the requested uniform 0.5 on all
three axes and a cell the grid still tiles.

Three breaks in this section fail without raising: these two spacing calls, and
the unfixed write path above where `OEWriteGrid` returns `True` over a wrong
file, so a migrator who verifies writes by checking that return value gets a
false pass. The other five items in this section are removals, and a removal
announces itself: the method is absent from `OESkewGrid`, and an `OEScalarGrid`
handed to a maptitude function raises `TypeError` at the boundary. Where a
replacement's semantics differ from what it replaces, as under **Containment**
and **Box geometry**, the item says so. Of the 46 public methods the two carriers
share, three are named for a scalar spacing: the two above and `IsSpacingSet()`,
which answered `True` on an anisotropic grid and so will not tell a caller the
map has three. Downstream code reading a spacing off a grid bound for maptitude
must derive three, from `get_grid_params(grid).x_spacing` and its siblings.

### Why recorded values moved

Nine corrections account for every recorded value this release moves. They are
numbered here so the pin table below can cite them.

1. **The carrier switch itself.** On an anisotropic map the scalar path applied
   one spacing to all three axes, resampling values onto a lattice the experiment
   never produced and drifting node positions along the coarser axes. Per-axis
   sampling removes both. Scores on an anisotropic map move; scores on an
   isotropic map do not.
2. **The in-grid domain narrowed.** `grid_contains` is a true fractional-index
   test where bounding-box `IsInGrid` was not, so the outer half-spacing shell is
   now outside the grid, and the six prechecks in `src/Metric.cpp` reject points
   they used to accept.
3. **Interpolation is maptitude's own trilinear kernel, with the node span closed
   on the far face.** `interpolate_density` returns a blend at `f == n - 1` where
   OpenEye returned the caller's default. The blend also runs in `double` where
   the OpenEye call took `float` coordinates and returned a `float`, a
   float-precision residual that does not need an anisotropic map to appear.
   Measured over all 216 `GUARD_SAMPLES` in
   `tests/cpp/test_interpolation.cpp`, on the synthetic 4x4x4 polynomial field
   those samples are defined against: 188 moved and 28 came back bit-identical
   to OpenEye's `float`. The largest residual is `6.9227e-07` taken relative to
   `max(1, abs(expected))`, the denominator that test uses, and `9.4116e-06`
   taken as an absolute difference. That field is not a crystallographic map;
   no equivalent sweep on one is recorded here.
4. **The Q-score radial step** is taken from the smallest of the three node
   intervals rather than from one scalar spacing.
5. **The calculated-density FFT sampling counts** are derived per axis, as
   `n_i = round(edge_i / spacing_i)`.
6. **The calculated-density FFT sampling origin** is the first node rather than
   the bounding-box edge — a half-spacing move.
7. **`wrap_and_pad_grid`'s padding decision is tested against the node bounds**
   `[origin_i, origin_i + (n_i - 1) * spacing_i]` rather than the half-spacing
   box. The old test anchored on `GetXMin()` — the box edge, half a spacing below
   node 0 — and added the node-span width to it, so relative to that test both
   edges move *up* by half a spacing: the tested region is translated by half a
   node interval per axis and keeps its width. `needs_pad` can therefore flip in
   either direction, `false` to `true` at the low edge and `true` to `false` at
   the high edge. **No case in the test suite flipped.** The three tests that
   encode the decision as a literal assertion —
   `WrapAndPadGridNoShiftNeeded` (`false` before, `false` after),
   `WrapAndPadGridCreatesPaddedGrid` (`true`, `true`) and
   `WrapAndPadReturnsNullptrOnlyWhenNoPaddingIsNeeded` (`false`, `false`) — all
   decide as they did before. This item moved no pin: no pin-producing path
   reaches the function.
8. **Both periodic wraps anchor at the node origin rather than the bounding-box
   edge**, a half-spacing shift in where the periodic image lands. The two sites
   are `interpolate_density_periodic` and the padded rebuild inside
   `wrap_and_pad_grid`. The second is the one with a visible consequence: the
   box-edge anchor pushed part of the padded grid outside the source and filled
   it with the default, so a padded grid now carries real density where it
   carried zeros. Measured on the three committed validation structures, read as
   `tests/python/test_validation.py` reads them and then padded by the shipped
   C++ `wrap_and_pad_grid`, **no node of any of the three padded grids is
   default-filled**: 0 values of exactly zero out of 32 422 for `1d26`, 20 592
   for `3q9g` and 66 780 for `340d`. None of the three sources holds an exact
   zero either, so every padded node carries interpolated density. A default fill
   is distinguishable: it is exactly `0.0`, and at a point a quarter of a spacing
   below `1d26`'s first node on all three axes `interpolate_density` returns
   that `0.0` while `interpolate_density_periodic` returns `0.0420`. How many
   nodes that is a change of is not measured: the box-edge anchor is no longer
   in the tree, so no
   shipped entry point produces the value it would have written, and the figures
   this paragraph previously carried came from a Python reimplementation of the
   rebuild rather than from either path as shipped. This item moved no pin
   either, for the same reason as item 7.
9. **`wrap_and_pad_grid` centred the molecule half a spacing off on every axis** —
   the bug fix recorded under **Fixed** above. Unpinned, like items 7 and 8.

**The pins that moved.** Six values in `tests/cpp/pin_values.h` changed, and no
others. The two Q-score pins are read on the in-memory fixture
`generate_pins.cpp` builds as `ObsGrid()` —
`MakeGaussianGrid(0.0, 0.0, 0.0, 1.0, 6.0, 0.5)`, so 25 nodes per axis running
-6.0 to 6.0, exactly isotropic at 0.5 A, centred on the origin — and the
consuming tests at `tests/cpp/test_metric_characterization.cpp:182` and `:190`
call that same builder. That exact isotropy is why items 1, 4 and 5 cannot move
a value there; the size of the move is the float-to-double residual of item 3.
The four `SHELLS4` pins come from a fixture that differs from the single-shell
`FC_ORTHORHOMBIC` case in one argument, `n_scale_shells` 4 against 1, and the
per-shell branch is the only place that samples the observed map from the FFT
origin item 6 moved; the single-shell pins beside it did not move. Item 3 fires
in that same loop — this release also swapped `interpolate_density` for
`interpolate_density_at` there — but it is a float-to-double residual where item
6 is a half-spacing geometric move. The smallest of these four,
`FC_ORTHORHOMBIC_SHELLS4_SUM` at `+0.0031862` relative, is 4 602.5 times the
largest residual the guard sweep under item 3 found, 3.66 orders, and 251.5
times that sweep's worst case against the plain `abs(expected)` denominator,
`1.267e-05`, which is 2.40 orders. That sweep ran on the synthetic guard field
and not on this fixture, so what it fixes is the scale a float-to-double
narrowing works at, not a bound on this loop. A 0.32% move sits 2.40 to 3.66
orders past that scale depending on the convention, so the attribution column
stays at 6.

| pin | from | to | relative change | correction |
|---|---|---|---|---|
| `QSCORE_CARBON_DEFAULT` | `0.96296527998254844` | `0.96296527987682856` | `-1.1e-10` | 3 |
| `QSCORE_CARBON_SIGMA08` | `0.9943922001151142` | `0.99439220042494902` | `+3.1e-10` | 3 |
| `FC_ORTHORHOMBIC_SHELLS4_SUM` | `13893.220623970032` | `13937.487342774868` | `+0.32%` | 6 |
| `FC_ORTHORHOMBIC_SHELLS4_SUM_SQ` | `12354.318438706807` | `12433.449584055335` | `+0.64%` | 6 |
| `FC_ORTHORHOMBIC_SHELLS4_MIN` | `0.82632350921630859` | `0.82111889123916626` | `-0.63%` | 6 |
| `FC_ORTHORHOMBIC_SHELLS4_INDEX_MOMENT` | `-92984.176744340395` | `-115006.82435095034` | `-23.7%` | 6 |

`EDIAM_*`, `COVERAGE_*` and the `FC_*` pins outside the `SHELLS4` group did not
move, so this release does not shift every metric. No Tier 1 analytic test was
repinned: `tests/cpp/test_metric_analytic.cpp` reads no pin, and its assertions
pass unchanged. A value moving below a characterization test's tolerance is not
visible in this table — `tests/cpp/grid_summary.h` compares at `1e-6` — so "no
pin moved" means "no pin moved by more than that", not "nothing changed".

### Performance

Interpolation was characterized on the 1d26 asset in a Release build: mean
**8.283 ms** for the OpenEye kernel against **7.542 ms** for maptitude's, a mean
ratio of **1.0995** over 15 samples, with roughly ±10% run-to-run noise. The
ratio sits inside that noise band, so the measurement supports **no regression,
and roughly 10% faster on this asset** — it does not support a precise
multiplier.

## [0.3.0]

Foundation and safety. This release adds input validation, a typed exception
hierarchy, memory and resource safety fixes, and the project's first C++ test
suite. Measured against the release it follows, `v0.2.4` (`5848b0d`), that is 75
commits touching 48 files as of the final Phase 1 fix round; re-derive with
`git rev-list --count v0.2.4..` and `git diff --name-only v0.2.4.. | wc -l`
rather than trusting the number, which drifts with every later commit. It is not
an accuracy release: no change here was intended to make a metric more correct,
and the six exceptions to that intent are listed under Exceptions to the
neutrality claim below.

### Removed

- **Non-orthorhombic unit cells are no longer supported.** `DensityCalculator`
  now raises `CellError` for monoclinic and triclinic cells. This is a
  capability regression, not a bug fix: those cells previously returned a
  value, but the structure-factor pipeline computes `1/d^2` as
  `(h/a)^2 + (k/b)^2 + (l/c)^2`, which is only correct for an orthorhombic
  lattice, so the returned value was wrong. General lattice support is planned;
  until then an exception is preferable to a plausible wrong answer.
- **`rscc` no longer accepts `AtomRadius::ADAPTIVE`.** It now raises
  `std::invalid_argument` (`RuntimeError` in Python). RSCC never had an adaptive
  radius model: its radius switch listed `FIXED`, `SCALED`, and
  `BINNED: default:`, so `ADAPTIVE` fell through to the binned radius and
  returned the binned score under another name, with nothing on the result
  recording which model had run. The withdrawn pin `RSCC_CARBON_ADAPTIVE` was
  bit-identical to `RSCC_CARBON_BINNED`, which is what its equality had been
  pinning; the characterization test now asserts the rejection instead. The
  switch also lost its `default:` label, so a future enumerator is a compiler
  warning rather than another silent fall-through. `rsr` is unaffected and
  still defaults to `ADAPTIVE`.

### Added

- A typed exception hierarchy: `MaptitudeError` and its four subclasses
  `StructureError`, `GridError`, `SymOpError`, and `CellError`. Previously all
  three C++ exception types surfaced in Python as `ValueError`, so callers had
  to match on message text. See the Exceptions section of the README.
- Unit-cell validation at construction and at every point a cell is consumed.
  Cells with non-positive or non-finite lengths, angles outside `(0, 180)`
  degrees, a non-positive volume radicand, or derived matrices too
  ill-conditioned to round-trip a coordinate are now rejected with `CellError`.
- `QScoreOptions` setters reject invalid values. A non-positive radial step
  previously caused a non-terminating loop.
- A bound on the Q-score radial sweep. The shell count, the per-shell sample
  count, and their product are each capped, and the step and maximum radius are
  each required to be usable, so no combination of options can exhaust memory
  or hang before scoring begins. An unusable sweep raises `GridError` naming
  the option to adjust — under a bound as often as over one, since a step too
  small to advance the accumulator and a radius of zero both fail here. In
  `RadialSampling::ADAPTIVE` the maximum radius is derived from the atom rather
  than from the options, so an atom that carries no radius scores `NaN` and the
  rest of the molecule is still scored. A sweep that *no* atom in the molecule
  can satisfy is a statement about the resolution and the grid spacing instead,
  and raises `GridError`: an adaptive Q-score at a grid spacing of 4.0 A, where
  no atom has a radius large enough to produce even one shell, is a caller error
  rather than a molecule of unscorable atoms.
- Validation of the public numeric arguments that reached unchecked arithmetic.
  `resolution` must be finite and positive at all five entry points; NaN and
  `+inf` both passed the old `resolution <= 0.0` test, and what the argument
  reaches past that test differs by path. `ediam` and
  `DensityCalculator::Calculate` divide by the resolution or its square; `rscc`
  and `rsr` size the scoring radius from it under their default radius models;
  `qscore` derives the adaptive sweep's step from it. `qscore`'s fixed sweep,
  which is the default, never reads the argument at all, and is checked for
  consistency with the adaptive path and with the other four entry points rather
  than because it divides. `n_scale_shells` must be in
  `[1, 1000]` (`MAX_SCALE_SHELLS`); at `UINT_MAX` the `n_scale_shells + 1` that
  sizes the shell-edge table wrapped to zero and the loop filling it could not
  terminate. A grid spacing at or above twice a cell edge is rejected: it rounds
  that FFT dimension to zero, making the Miller-index wrap a division by zero.
  `wrap_and_pad_grid` validates all three cell edges, which it had passed
  straight to `std::fmod` — NaN for a zero divisor, and that NaN reached every
  voxel of the padded grid. These raise `GridError`, except the cell edges,
  which raise `CellError`.
- `CoverageOptions::SetSigma` rejects a non-finite sigma, which made the density
  threshold NaN or infinite. Measured at `v0.2.4`, NaN and `+inf` returned a
  plausible `overall = 0.0`; `-inf`, on a map with nonzero spread, a perfect
  `1.0`. It rejects *only* non-finite values. `QScoreOptions::SetSigma`
  additionally refuses zero and negatives, and the asymmetry is deliberate: the
  two sigmas are different quantities. Q-score's is a Gaussian width, which has
  to be positive; coverage's is a multiplier in the threshold
  `mean + sigma * stddev`, where zero means "threshold at the mean" and a
  negative value means "threshold below the mean", both meaningful requests.
- A bound on the Miller-index box, `MAX_MILLER_BOX_POINTS` (2e8 points). The
  index sweep spans `2*ceil(edge/resolution)+1` points per axis, so its volume
  grows as `(2a/resolution)^3` and a finite positive resolution does not bound
  it: 1e-9 A is an ordinary double, not a subnormal, and the loop does not
  finish. The bound is on the box rather than on the resolution because the two
  depend on the three cell edges jointly. A 200 A cubic cell at 1.0 A resolution
  uses under a third of it. Raises `GridError` naming the resolution, the three
  edges, and the limit.
- The `AtomRadius` and `RadialSampling` setters reject values the enum does not
  declare. Both are scoped enums with underlying type `int`, so
  `SetAtomRadiusMethod(42)` was a valid value of the type and SWIG passed it
  through without a cast. It matched no arm of the RSCC radius switch, which
  then read its radius uninitialized and returned `overall = 1.0` with every
  per-atom score `NaN`; and it took the `else` arm of both of `qscore`'s
  `RadialSampling` tests, so it was validated as adaptive and then executed as
  fixed, with the fixed sweep's parameters never checked. Raises
  `std::invalid_argument` (`RuntimeError` in Python).
- A C++ test suite covering the metrics, the structure-factor pipeline, and
  real CCP4 map I/O. Of its 28 characterization tests, 26 assert pinned values.
  The other two assert rejections and replaced the pins withdrawn by the two
  changes under Removed: five `FC_MONOCLINIC_*` for the monoclinic cell, and
  `RSCC_CARBON_ADAPTIVE` for the adaptive RSCC radius. It runs in CI on pushes
  to `master`, on every pull request, and on manual dispatch.

### Changed

- **Requires OpenEye Toolkits 2026.1 or later** (was 2025.2).
- `coverage(sigma=...)` defaults to `None` rather than `1.0`, meaning "use the
  value carried by `options`". Passing `sigma=1.0` explicitly was previously
  indistinguishable from omitting it, so an explicit `1.0` was silently
  discarded whenever `options` carried a different sigma and the options value
  set the threshold instead. The explicit argument now wins. Calls that pass
  only one of the two are unaffected.
- `rscc`, `rsr`, and `coverage` no longer mutate the options object passed to
  them. Overrides such as `atom_radius` and `sigma` are applied to a copy.
- `combine_maps` and `diff_to_calc` require exactly matching grid geometry --
  dimensions, origin, and spacing -- via `OEGridSameGeometry`. Grids that
  differed only in origin were previously combined element-wise, mixing
  densities from different points in space. Spacing comparison is now exact
  rather than within `1e-6`.
- `wrap_and_pad_grid` raises `StructureError` when the molecule contains no
  heavy atoms. In the C++ API its `nullptr` return now means only "no padding
  was needed". The Python wrapper never exposes that `nullptr`: it substitutes
  the original grid, so it always returns a grid and never `None`.
- Grids returned across the Python boundary are independently owned copies.
  Peak memory doubles for the duration of the copy; an allocation failure
  during the copy raises `MemoryError`, and a grid whose geometry does not
  survive copy-assignment raises `RuntimeError`.
- Passing an object that is not an OpenEye atom predicate as a `mask` raises
  `TypeError` rather than being reinterpreted as a pointer.
- The symmetry-operator parser rejects input it previously accepted or
  mis-parsed, raising `SymOpError` for a repeated axis within one component
  (`x+x,y,z`), consecutive or trailing signs, empty components, zero and
  non-finite denominators, hexadecimal and signed number tokens, translations
  that overflow when summed, and missing separators (`2x`, `xy`, `x.5`,
  `x+2y`).
- `SymOp::ToString` output changed in four ways, which matter to anyone parsing
  the string rather than round-tripping it through the parser: the decimal
  fallback now uses `max_digits10`; a zero component serializes as `0`;
  translations past the `int` fraction range serialize as decimals; and
  translations in `1e-10 < |t| < 5e-9` are suppressed, so `x+1e-9,y,z`
  round-trips to `x,y,z`.
- Nine previously unchecked FFTW allocations and five plan creations now raise
  `GridError` on failure instead of dereferencing null.

### Fixed

- A cross-allocator double free in the grid return path: the SWIG wrapper
  deleted an object allocated inside the OpenEye shared library and swapped in
  its own pointer, so each runtime freed the other's memory.
- The atom-mask argument is type-checked before its `this` pointer is cast
  across the SWIG runtime boundary. A wrong-typed mask previously crashed the
  interpreter.
- FFTW buffers and plans in `DensityCalculator::Calculate` are managed with
  RAII, so an exception mid-pipeline no longer leaks them. Plan creation and
  destruction are serialized: the FFTW planner is not thread-safe.
- Undefined behavior in the symmetry-operator parser: `std::isdigit` was called
  with a possibly-negative `char`.

### Exceptions to the neutrality claim

Every pinned metric value in this release is bit-identical to the value pinned
before any behavior change, verified by regenerating `tests/cpp/pin_values.h`
from the post-fix binary and diffing. The only pins removed are the five
`FC_MONOCLINIC_*` entries withdrawn by exception 1.

Six changes are exceptions to that intent. Four were planned; exceptions 5 and
6 were found by the reconciliation audit and its review, after the work had
landed.

The new rejections listed under Added are not among them. With one measured
exception, none of them refuses input that previously produced a correct value:
the inputs they refuse variously hung, exhausted memory, divided by zero, read
uninitialized memory, or ran to completion and returned NaN or a plausible wrong
number. Those mechanisms illustrate rather than enumerate. What the neutrality
claim rests on is the property, and a mechanism not in that list would not
disturb it.

The measured exception is `qscore`'s fixed radial sweep, which never reads the
`resolution` argument: at `v0.2.4` it returned a bit-identical score at 0.5,
1.0, 2.0, 3.5, 10.0 and `+inf` A, so `qscore(mol, grid, +inf)` did return a
correct Q-score and is now refused. Whether the argument is read can turn on the
configuration rather than on the entry point alone: `rscc` and `rsr` also ignore
it under `AtomRadius::FIXED` and `AtomRadius::SCALED`, at `v0.2.4` as here, since
both models take the radius from the options or from the atom without consulting
it. What distinguishes `qscore`'s fixed sweep is that it is the configuration a
caller reaches by default, `rscc` defaulting to `BINNED` and `rsr` to `ADAPTIVE`,
both of which size the radius from the resolution. The check is kept because the
configurations that do read it raised no error at `+inf` either: finite scores
from `rscc` and `rsr` under those defaults and from `qscore`'s adaptive sweep, a
finite grid from `DensityCalculator::Calculate`, and `NaN` from `ediam`. A
resolution no caller can mean is better reported at the call than honoured in one
configuration and ignored in the next.

The criterion for the list below is narrower than "refuses something", and its
operative half is the input. Exceptions 1, 2, and 4 each refuse a well-formed
request that a caller can legitimately mean and that previously got an answer: a
monoclinic cell, a molecule with no heavy atoms, two grids whose spacings differ
within a tolerance. Whether the answer was right is not the test -- exception 1
is on the list precisely because its answer was wrong. A malformed argument
fails the other half, and that is what keeps the two rejections with the most
plausible before-values under Added rather than here. Measured at `v0.2.4`,
`coverage` returned `overall = 0.0` for a NaN `sigma` and `1.0` for
`sigma = -inf`, both indistinguishable from real scores; `resolution = +inf`
returned the numbers above. Neither a non-finite multiplier nor an infinite
resolution is a request a caller can mean, and the two numbers arise
differently. The `+inf` Q-score is the carve-out above because the path that
produced it never reads the resolution: the argument did not enter the
computation, so the number that came back was the score. The sigma does enter
it, straight into `mean + sigma * stddev`, and what it yields there turns on the
arithmetic of infinity rather than on the map -- `-inf` is below every density
only where the spread is nonzero, and makes the threshold NaN where it is zero.
A `1.0` that holds only for maps with spread is an artifact of the argument
rather than a coverage, so it is not a second carve-out.

The adaptive no-atom-can-sweep `GridError` is the only new rejection that
refuses well-formed input which previously ran to completion without an error.
Measured at `v0.2.4`, an adaptive Q-score on a single carbon at resolution
30.0 A and grid spacing 4.0 A returned `overall = NaN` with
`by_atom = {0: NaN}`: the step is `min(spacing, resolution / 7)` = 4.0 A against
a maximum radius of `2 * 1.7` = 3.4 A, so the sweep spans no shell, only the
replicated centre point is sampled, and `pearson_correlation` reports the
resulting zero variance as NaN rather than dividing by zero. Its input is a
well-formed grid, but a NaN is not an answer, so it fails the second half of the
criterion and is not a seventh exception.

1. **Non-orthorhombic cells raise `CellError`.** A capability regression, as
   described under Removed. This is the only exception that removes a pin.
2. **`wrap_and_pad_grid` raises `StructureError` on an empty heavy-atom set.**
   `nullptr` from that function now means only "no padding needed"; the return
   value is no longer overloaded. Changes rejection behavior, not scores.
3. **`coverage(sigma=...)` no longer discards an explicit `sigma=1.0`.** Scores
   move for exactly one call shape: `coverage(sigma=1.0, options=<options
   carrying a different sigma>)`. Only `1.0` was affected, because the old
   default was `1.0` and the wrapper tested `if sigma != 1.0` to infer whether
   the caller had supplied the argument at all; any other explicit value
   already reached `SetSigma` and already won. No C++ pin exercises that path,
   because the coverage pins are generated from C++ entry points that never
   route through the Python wrapper.
4. **Grid combination requires exact geometry equality.** Near-equal spacings
   and mismatched origins that were previously accepted are now rejected.
   Changes rejection behavior, not scores.
5. **`wrap_and_pad_grid`'s heavy-atom predicate changed** from a hand-rolled
   `GetAtomicNum() == 1` skip to `OEChem::OEIsHeavy()`, in both the centroid
   and the bounding-box loop. The two disagree on Z=0: the old code included
   dummy atoms and virtual sites in the centroid and bounding box, the new code
   excludes them. **Affected input class: molecules carrying Z=0 atoms
   alongside real heavy atoms.** For those, the centroid, the wrap translation
   derived from it, the bounding box, and every metric computed on the
   resulting grid all move. Measured for two carbons at (-1,0,0) and (1,0,0)
   plus one Z=0 atom at (20,20,20), cell 30x30x30, padding 3.0, spacing 1.0:
   the new code returns a 9x7x7 grid centred at (0,0,0), the old behavior a
   28x27x27 grid centred at (9.5,10,10) -- a 46-fold difference in voxel count
   and a displaced centre. No pin and no test covers the mixed case.
6. **`rscc`, `rsr`, and `coverage` no longer mutate the caller's options
   object.** A caller reusing one options object across calls, where an earlier
   call passed `atom_radius` or `sigma`, previously had that setting leak into
   every subsequent call. Removing the leak changes returned scores for those
   callers: one measured pair of calls moved from `0.9999999999999996` to
   `1.0000000000000004`, and a second, on a different grid, from
   `0.9592678096783006` to `0.9503942559650402`. No pin covers this path,
   because the pin corpus is C++-only and never crosses the SWIG boundary.

**What the neutrality claim does and does not bound.** It is a claim about
*pinned* values, not about all values. A bit-identical regeneration bounds
drift over the corpus the pins actually cover and is silent everywhere else --
exception 5 moves real numbers and the regeneration is still identical, because
no pinned structure exercises the mixed Z=0 case. The sharpest limit is not a
coverage percentage but the shape of the corpus: no pin crosses the SWIG
boundary, which is exactly where exceptions 3 and 6 live. For reference, the
last recorded line coverage of `src/` by the whole suite was 84.71%, measured
several dozen commits ago and before much of the new validation code landed;
treat it as a historical snapshot rather than a current floor, and note that it
is suite coverage, not pin coverage.

### Known limitations

- A writable `this` on SWIG proxies defeats cross-runtime pointer validation.
  maptitude's type checks confirm the Python wrapper's real type, but `this` is
  assignable on every SWIG proxy, so an accepted proxy can be made to carry a
  foreign pointer. maptitude cannot validate the pointee: its SWIG runtime
  version differs from the OpenEye modules' and neither exports a joinable type
  table. Passing a deliberately punned proxy is undefined behavior.
- The FFTW planner mutex covers only maptitude's own planning calls. A host
  application that calls FFTW planning routines on another thread remains
  unsynchronized with respect to maptitude, and the library cannot close that
  gap from the inside.
- Several validation guards are pinned only by `EXPECT_THROW`, which is
  structurally blind to a narrowed accepting set and to wrong output on the
  accepting path. `RadialSampling::ADAPTIVE` is the notable case. Three tests
  now assert that an adaptive sweep on valid input returns a real score rather
  than `NaN`, which closes the narrowed-accepting-set half, but no pin records
  an adaptive Q-score: `RSR_CARBON_ADAPTIVE` pins the unrelated
  `AtomRadius::ADAPTIVE`. The accepting path is covered for shape and not for
  value, so a change that moved every adaptive score would leave the suite
  green.
- One branch of `SymOp::ToString` is unpinnable by construction. Suppressing a
  zero-numerator translation makes the distinguishing input for the separator's
  sign decision unreachable, so the two candidate implementations are
  equivalent mutants. The separator itself is pinned; only the source of its
  sign decision is not.
- The FFTW RAII conversion and the nine allocation null checks it added are not
  covered by a regression test. Every test that throws from `Calculate` throws
  before the first `fftw_alloc_complex` — the argument checks and the
  Miller-index bound all fire during setup — so reverting the `FftwBuffer` and
  `FftwPlan` wrappers to raw allocation plus manual `fftw_free` on the success
  path leaves the whole suite green. The evidence for the conversion is a
  one-time manual leak measurement taken during the phase, not something CI
  re-checks. Injecting an FFTW allocation failure from a test would need an
  allocator seam the library does not have.
- The library disagrees with itself about what a heavy atom is. Five sites in
  `DensityCalculator.cpp` and `Metric.cpp` still use a hand-rolled
  `GetAtomicNum() == 1` skip while `Metric.cpp` and `GridOps.cpp` use
  `OEChem::OEIsHeavy()`. The two disagree on Z=0 atoms. Contracting the
  remaining five would move metric values for structures containing Z=0 atoms
  alongside real heavy atoms, and characterization coverage for that case does
  not exist yet, so the change is deferred.
