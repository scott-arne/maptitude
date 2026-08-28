# Changelog

All notable changes to this project are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project is pre-1.0: breaking changes may land in a minor release.

## [0.3.0]

Foundation and safety. This release adds input validation, a typed exception
hierarchy, memory and resource safety fixes, and the project's first C++ test
suite, across 67 commits touching 43 files. It is not an accuracy release: no
change here was intended to make a metric more correct, and the six exceptions
to that intent are listed under Exceptions to the neutrality claim below.

### Removed

- **Non-orthorhombic unit cells are no longer supported.** `DensityCalculator`
  now raises `CellError` for monoclinic and triclinic cells. This is a
  capability regression, not a bug fix: those cells previously returned a
  value, but the structure-factor pipeline computes `1/d^2` as
  `(h/a)^2 + (k/b)^2 + (l/c)^2`, which is only correct for an orthorhombic
  lattice, so the returned value was wrong. General lattice support is planned;
  until then an exception is preferable to a plausible wrong answer.

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
  count, and their product are each capped, so no combination of options can
  exhaust memory or hang before scoring begins. Sweeps over the limit raise
  `GridError` naming the option to adjust.
- A C++ test suite covering the metrics, the structure-factor pipeline, and
  real CCP4 map I/O, of which 28 tests carry pinned values. It runs in CI on
  every push and pull request.

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
  heavy atoms. Its `nullptr` return now means only "no padding was needed".
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

1. **Non-orthorhombic cells raise `CellError`.** A capability regression, as
   described under Removed. This is the only exception that removes a pin.
2. **`wrap_and_pad_grid` raises `StructureError` on an empty heavy-atom set.**
   `nullptr` from that function now means only "no padding needed"; the return
   value is no longer overloaded. Changes rejection behavior, not scores.
3. **`coverage(sigma=...)` no longer discards an explicit sigma.** Scores move
   for exactly one call shape: `coverage(sigma=X, options=<options carrying a
   different sigma>)`. No C++ pin exercises that path, because the coverage
   pins are generated from C++ entry points that never route through the
   Python wrapper.
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
  accepting path. `RadialSampling::ADAPTIVE` is the notable case: its only
  test asserts that a bad input is rejected, and nothing asserts a value from
  an adaptive sweep on valid input.
- One branch of `SymOp::ToString` is unpinnable by construction. Suppressing a
  zero-numerator translation makes the distinguishing input for the separator's
  sign decision unreachable, so the two candidate implementations are
  equivalent mutants. The separator itself is pinned; only the source of its
  sign decision is not.
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
