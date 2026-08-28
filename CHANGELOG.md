# Changelog

All notable changes to this project are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project is pre-1.0: breaking changes may land in a minor release.

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
  `+inf` both passed the old `resolution <= 0.0` test, and every entry point
  divides by the resolution or its square. `n_scale_shells` must be in
  `[1, 1000]` (`MAX_SCALE_SHELLS`); at `UINT_MAX` the `n_scale_shells + 1` that
  sizes the shell-edge table wrapped to zero and the loop filling it could not
  terminate. A grid spacing at or above twice a cell edge is rejected: it rounds
  that FFT dimension to zero, making the Miller-index wrap a division by zero.
  `wrap_and_pad_grid` validates all three cell edges, which it had passed
  straight to `std::fmod` — NaN for a zero divisor, and that NaN reached every
  voxel of the padded grid. These raise `GridError`, except the cell edges,
  which raise `CellError`.
- `CoverageOptions::SetSigma` rejects a non-finite sigma, which made the density
  threshold NaN, every `rho >= threshold` comparison false, and coverage return
  a plausible `0.0`. It rejects *only* non-finite values. `QScoreOptions::SetSigma`
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
correct Q-score and is now refused. The check is kept because every other path
does read the argument, including `qscore`'s own adaptive sweep, and none of
them raised an error at `+inf`: finite scores from `rscc`, `rsr`, and the
adaptive sweep, a finite grid from `DensityCalculator::Calculate`, and `NaN`
from `ediam`. A resolution no caller
can mean is better reported at the call than ignored in one code path and
honoured in another.

The criterion for the list below is narrower than "refuses something", and its
operative half is the input. Exceptions 1, 2, and 4 each refuse a well-formed
request that a caller can legitimately mean and that previously got an answer: a
monoclinic cell, a molecule with no heavy atoms, two grids whose spacings differ
within a tolerance. Whether the answer was right is not the test -- exception 1
is on the list precisely because its answer was wrong. A malformed argument
fails the other half, and that is what keeps the two rejections with the most
plausible before-values under Added rather than here. Measured at `v0.2.4`,
`coverage` with a NaN `sigma` returned `overall = 0.0`, indistinguishable from a
real score; `resolution = +inf` returned the numbers above. A NaN multiplier and
an infinite resolution are not requests, and those numbers were artifacts of an
argument that was never read rather than answers to one.

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
- The checked-in SWIG proxy `python/maptitude/maptitude.py` carries version
  constants and a `__version__` string that lag this release until the next
  build regenerates it. This is a consequence of tracking a generated file, not
  a defect in the release. It is invisible to callers: `maptitude.__version__`
  is assigned by `python/maptitude/__init__.py` and never read from the proxy.
