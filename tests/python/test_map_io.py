"""Python-side tests for the CCP4/MRC map I/O bindings.

The C++ tests pin the values; these check that the SWIG boundary hands them
across intact -- that the grid arrives as a native ``openeye.oegrid.OESkewGrid``
rather than a wrapper, that the symmetry text survives as ``str``, and that the
typed exceptions still land as their Python classes.
"""

from __future__ import annotations

import ctypes
import struct
import subprocess
import sys
from pathlib import Path

import maptitude
import numpy as np
import pytest
from maptitude import (
    CellError,
    GridError,
    MapFile,
    OriginSource,
    SymOpError,
    UnitCell,
    fc_density,
    get_grid_params,
    get_unit_cell,
    parse_symops,
    read_map,
    wrap_and_pad_grid,
)
from openeye import oechem, oegrid

_ASSET_DIR = Path(__file__).resolve().parents[1] / "assets" / "mapq"
_DATA_DIR = Path(__file__).resolve().parents[1] / "data"

_EXPECTED_RSCC = 0.9690309706700857
_EXPECTED_RSR = 0.09864796027088141
_EXPECTED_QSCORE = 0.9052790577235651

_CCP4_HEADER_BYTES = 1024
_SYMOP_RECORD_BYTES = 80

# Every malformed shape the C++ parser pin covers, from
# tests/cpp/test_map_io.cpp:551.
_MALFORMED_SYMOPS = ["not a symop at all", "x,y", "x,y,z,w", "x,y,q*z"]


def test_read_map_returns_a_native_grid_and_symop_text():
    result = read_map(_ASSET_DIR / "1d26_2fofc.ccp4")

    assert isinstance(result, MapFile)
    assert isinstance(result.grid, oegrid.OESkewGrid)
    assert isinstance(result.symops, str)

    assert (result.grid.GetXDim(), result.grid.GetYDim(),
            result.grid.GetZDim()) == (49, 49, 25)
    assert result.grid.GetSpaceGroup() == 96
    assert len(parse_symops(result.symops)) == 8
    assert result.symops.splitlines()[0] == "x,y,z"


def test_map_file_unpacks_as_a_pair():
    grid, symops = read_map(_ASSET_DIR / "340d_2fofc.ccp4")
    assert grid.GetXDim() == 61
    assert len(parse_symops(symops)) == 8


def test_read_map_accepts_a_string_path():
    result = read_map(str(_DATA_DIR / "test_map.ccp4"))
    assert result.grid.GetXDim() == 21
    assert result.symops == ""


def test_read_map_places_the_em_map_at_its_origin_record():
    result = read_map(_ASSET_DIR / "390_emd_30342_A_z4.mrc")
    params = get_grid_params(result.grid)
    assert params.x_origin == pytest.approx(145.825, abs=1e-3)
    assert params.y_origin == pytest.approx(112.825, abs=1e-3)
    assert params.z_origin == pytest.approx(120.517, abs=1e-3)


def test_read_map_carries_the_unit_cell_on_the_grid():
    # 390_emd_30342_A_z4.mrc is the only shipped fixture whose three edge
    # lengths differ, so a defect that copied b from a -- or defaulted any of
    # the three -- goes red here where a tetragonal fixture would not.
    #
    # The three angle assertions are a completeness check, not a
    # discriminating one: read_map rejects a nonorthogonal cell outright, so no
    # fixture can carry a non-90 angle through this path and an implementation
    # that hardcoded 90 would pass here.  That rejection is pinned separately,
    # by test_read_map_rejects_a_skewed_cell_with_a_typed_error below.
    cell = get_unit_cell(read_map(_ASSET_DIR / "390_emd_30342_A_z4.mrc").grid)
    assert cell.a == pytest.approx(99.057, abs=1e-3)
    assert cell.b == pytest.approx(90.153, abs=1e-3)
    assert cell.c == pytest.approx(70.119, abs=1e-3)
    assert cell.alpha == pytest.approx(90.0, abs=1e-4)
    assert cell.beta == pytest.approx(90.0, abs=1e-4)
    assert cell.gamma == pytest.approx(90.0, abs=1e-4)


def test_read_map_rejects_a_skewed_cell_with_a_typed_error(tmp_path):
    # The rejection happens inside the out-typemap's grid-copy validation,
    # after the %exception-protected C++ call has already returned, so the
    # exception class is chosen by hand there.  Before this was fixed the
    # failure arrived as builtins.RuntimeError and `except MaptitudeError`
    # walked straight past it.
    #
    # The control half is what makes this a measurement rather than a
    # tautology: the same six-float patch with every angle left at 90 reads
    # back exactly, so the rejection is attributable to beta and not to the
    # patching.
    raw = bytearray((_ASSET_DIR / "1d26_2fofc.ccp4").read_bytes())
    struct.pack_into("<6f", raw, (11 - 1) * 4, 30.0, 40.0, 50.0, 90.0, 105.0, 90.0)
    skewed = tmp_path / "beta_105.ccp4"
    skewed.write_bytes(bytes(raw))

    with pytest.raises(CellError) as caught:
        read_map(skewed)
    assert "axis-aligned" in str(caught.value)
    assert issubclass(CellError, maptitude.MaptitudeError)

    struct.pack_into("<6f", raw, (11 - 1) * 4, 30.0, 40.0, 50.0, 90.0, 90.0, 90.0)
    control = tmp_path / "beta_90.ccp4"
    control.write_bytes(bytes(raw))
    cell = get_unit_cell(read_map(control).grid)
    assert cell.a == pytest.approx(30.0, abs=1e-4)
    assert cell.b == pytest.approx(40.0, abs=1e-4)
    assert cell.c == pytest.approx(50.0, abs=1e-4)


def test_read_map_reports_a_failed_geometry_copy_as_grid_error(tmp_path):
    # The other arm of the branch the skewed-cell test above drives.  Round 4
    # made the copy-validation block choose between CellError and GridError;
    # this pins the choice rather than one side of it.  Before round 4 both
    # arms arrived as builtins.RuntimeError.
    #
    # A zero cell edge is what reaches this arm: the copy succeeds, then
    # same_grid_geometry cannot derive a node interval from the degenerate
    # axis.  OpenEye prints "SetXMid unable to handle NaN" warnings to stderr
    # on this input; they come from the toolkit's own setters and are not a
    # failure.
    raw = bytearray((_ASSET_DIR / "1d26_2fofc.ccp4").read_bytes())
    struct.pack_into("<f", raw, (11 - 1) * 4, 0.0)
    degenerate = tmp_path / "a_zero.ccp4"
    degenerate.write_bytes(bytes(raw))

    with pytest.raises(GridError) as caught:
        read_map(degenerate)
    assert "Grid geometry copy failed" in str(caught.value)
    # The discriminating assertions.  CellError is a sibling of GridError, not
    # a parent, so a change that routed every copy failure to CellError would
    # still satisfy `pytest.raises(MaptitudeError)` -- but not this.  And
    # GridError is not a RuntimeError subclass, so this also fails if the
    # pre-round-4 behaviour comes back.
    assert not isinstance(caught.value, CellError)
    assert not isinstance(caught.value, RuntimeError)


def test_tiebreak_selects_which_record_wins(tmp_path):
    # Both records must be live for the tiebreak to have anything to break.
    # test_map.ccp4 carries NxSTART (-10, -10, -10), which the reader places at
    # node 0 = (-5, -5, -5); patching header words 50-52 adds a nonzero ORIGIN.
    # Mirrors the C++ pin at tests/cpp/test_map_io.cpp:456.
    assert int(OriginSource.ORIGIN_RECORD) != int(OriginSource.NXSTART)

    raw = bytearray((_DATA_DIR / "test_map.ccp4").read_bytes())
    struct.pack_into("<3f", raw, (50 - 1) * 4, 7.0, 8.0, 9.0)
    variant = tmp_path / "both_records_set.ccp4"
    variant.write_bytes(raw)

    by_origin = get_grid_params(read_map(variant, OriginSource.ORIGIN_RECORD).grid)
    assert by_origin.x_origin == pytest.approx(7.0, abs=1e-6)
    assert by_origin.y_origin == pytest.approx(8.0, abs=1e-6)
    assert by_origin.z_origin == pytest.approx(9.0, abs=1e-6)

    by_nxstart = get_grid_params(read_map(variant, OriginSource.NXSTART).grid)
    assert by_nxstart.x_origin == pytest.approx(-5.0, abs=1e-6)
    assert by_nxstart.y_origin == pytest.approx(-5.0, abs=1e-6)
    assert by_nxstart.z_origin == pytest.approx(-5.0, abs=1e-6)

    # The default must be ORIGIN_RECORD, not merely a value that parses.
    assert get_grid_params(read_map(variant).grid).x_origin == pytest.approx(
        7.0, abs=1e-6)


def test_read_map_raises_grid_error_for_a_missing_file():
    with pytest.raises(GridError):
        read_map(_DATA_DIR / "no_such_map.ccp4")


@pytest.mark.parametrize(
    "bad", [None, 3, ["m.ccp4"], b"tests/data/test_map.ccp4"])
def test_read_map_rejects_arguments_that_are_not_paths(bad):
    # str(path) turned each of these into a filename and reported a missing
    # file: None became "None", and the bytes case names a file that really
    # exists but stringifies to "b'tests/data/test_map.ccp4'".  GridError is
    # not a TypeError, so this goes red against the old body.  os.fspath
    # rejects the first three; the std::string typemap rejects the bytes.
    with pytest.raises(TypeError):
        read_map(bad)


def test_read_map_grid_survives_repeated_attribute_access():
    # .grid is tuple slot 0, so both reads return the same object by
    # construction.  What this pins is narrower: that the object the typemap
    # built is a live, usable grid and that MapFile keeps it alive for the
    # caller.  It says nothing about the released source grid's lifetime --
    # _maptitude_wrap_as_oe_skew_grid copy-assigns into a separate
    # Python-owned grid, so no Python-visible behaviour distinguishes a leaked
    # source from a freed one.  Ownership is established by reading that
    # helper, not here.
    result = read_map(_DATA_DIR / "test_map.ccp4")
    first = result.grid
    second = result.grid
    assert first is second
    assert first.GetSize() == 9261


@pytest.mark.parametrize("record", _MALFORMED_SYMOPS)
def test_read_map_raises_the_typed_symop_error(tmp_path, record):
    # The hierarchy check alone cannot show that SymOpError survives the SWIG
    # boundary: it holds whether or not read_map can raise and marshal one.
    # Overwriting record 1 of 1d26's eight-record block leaves the file length
    # and every header word untouched, so the symop parser is the only thing
    # that can reject the file -- and record 0 stays valid, so this pins the
    # parser rejecting a record rather than rejecting the block wholesale.
    assert issubclass(SymOpError, maptitude.MaptitudeError)

    raw = bytearray((_ASSET_DIR / "1d26_2fofc.ccp4").read_bytes())
    start = _CCP4_HEADER_BYTES + _SYMOP_RECORD_BYTES
    raw[start:start + _SYMOP_RECORD_BYTES] = record.encode().ljust(
        _SYMOP_RECORD_BYTES, b" ")
    variant = tmp_path / "malformed_symop.ccp4"
    variant.write_bytes(bytes(raw))

    with pytest.raises(SymOpError):
        read_map(variant)


def test_write_map_round_trips_a_crystallographic_asset(tmp_path):
    source = read_map(_ASSET_DIR / "1d26_2fofc.ccp4")
    out = tmp_path / "round_trip.ccp4"

    maptitude.write_map(out, source.grid, source.symops)

    back = read_map(out)
    assert (back.grid.GetXDim(), back.grid.GetYDim(),
            back.grid.GetZDim()) == (49, 49, 25)
    assert back.grid.GetSpaceGroup() == 96
    assert back.symops == source.symops

    before = get_grid_params(source.grid)
    after = get_grid_params(back.grid)
    assert after.x_spacing == pytest.approx(before.x_spacing, rel=1e-6)
    assert after.y_spacing == pytest.approx(before.y_spacing, rel=1e-6)
    assert after.z_spacing == pytest.approx(before.z_spacing, rel=1e-6)
    assert after.x_origin == pytest.approx(before.x_origin, abs=1e-6)

    original = source.grid.GetValues()
    written = back.grid.GetValues()
    assert [written[i] for i in range(back.grid.GetSize())] == [
        original[i] for i in range(source.grid.GetSize())
    ]


# Deliberately a child script rather than an in-process loop: the leak this
# pins is a process address, which is stable within one process, so a
# single-process comparison passes while the defect is present.
_REPRODUCIBILITY_CHILD = """
import sys

import maptitude

bundle = maptitude.read_map(sys.argv[1])
maptitude.write_map(sys.argv[2], bundle.grid, bundle.symops)
"""


def test_write_map_is_byte_reproducible_across_processes(tmp_path):
    """Two processes writing one grid must produce the same bytes.

    OEWriteGrid leaves SKWTRN, header words 35-37, holding whatever was in the
    memory behind them, and write_map hands it a copied grid. Word 35 then
    carried the low half of a heap pointer, which ASLR moves per process: the
    same grid written twice gave different files, so content hashing and diff
    failed spuriously, and four bytes of the writer's address space rode along
    in the file. PatchHeaderRecords zeroes those words.
    """
    source = _ASSET_DIR / "1d26_2fofc.ccp4"
    written = []
    for tag in ("first", "second"):
        out = tmp_path / f"{tag}.ccp4"
        result = subprocess.run(
            [sys.executable, "-c", _REPRODUCIBILITY_CHILD, str(source), str(out)],
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0, (
            f"{tag} child exited {result.returncode}:\n{result.stderr}"
        )
        written.append(out.read_bytes())

    first, second = written
    assert len(first) == len(second)

    # Report the offsets rather than the buffers: these files are a quarter of
    # a megabyte, and the header word is the diagnostic.
    differing = [index for index, (a, b) in enumerate(zip(first, second)) if a != b]
    header_words = sorted(
        {offset // 4 + 1 for offset in differing if offset < _CCP4_HEADER_BYTES}
    )
    past_header = sum(1 for offset in differing if offset >= _CCP4_HEADER_BYTES)
    assert not differing, (
        "two processes wrote different bytes for the same grid: "
        f"{len(differing)} differing bytes, 1-based header words "
        f"{header_words}, {past_header} bytes past the header"
    )


def test_write_map_carries_the_em_origin(tmp_path):
    source = read_map(_ASSET_DIR / "390_emd_30342_A_z4.mrc")
    out = tmp_path / "em.mrc"

    maptitude.write_map(out, source.grid)

    after = get_grid_params(read_map(out).grid)
    assert after.x_origin == pytest.approx(145.825, abs=1e-3)
    assert after.y_origin == pytest.approx(112.825, abs=1e-3)
    assert after.z_origin == pytest.approx(120.517, abs=1e-3)


def test_write_map_accepts_the_ccp4_and_mrc_spellings(tmp_path):
    source = read_map(_DATA_DIR / "test_map.ccp4")
    for name in ("out.ccp4", "out.mrc"):
        out = tmp_path / name
        maptitude.write_map(out, source.grid)
        assert read_map(out).grid.GetSize() == 9261


def test_write_map_defaults_an_absent_space_group_to_p1(tmp_path):
    source = read_map(_ASSET_DIR / "390_emd_30342_A_z4.mrc")
    assert source.grid.GetSpaceGroup() == 0
    out = tmp_path / "p1.mrc"
    maptitude.write_map(out, source.grid)
    assert read_map(out).grid.GetSpaceGroup() == 1


def test_write_map_raises_on_an_unwritable_extension(tmp_path):
    source = read_map(_DATA_DIR / "test_map.ccp4")
    out = tmp_path / "density.dat"
    with pytest.raises(GridError):
        maptitude.write_map(out, source.grid)
    assert not out.exists()


def test_write_map_raises_on_unparseable_symops(tmp_path):
    source = read_map(_DATA_DIR / "test_map.ccp4")
    out = tmp_path / "bad_symops.ccp4"
    with pytest.raises(SymOpError):
        maptitude.write_map(out, source.grid, "not a symop")
    assert not out.exists()


# pdb, resolution, node count. The resolutions are the ones
# tests/python/test_validation.py already pins in _PDB_RESOLUTIONS:113-117; the
# node counts are what read_map returns for the observed map, which fc_density
# copies onto its output.
_FC_CASES = [
    ("1d26", 2.12, (49, 49, 25)),
    ("3q9g", 2.05, (37, 37, 61)),
    ("340d", 1.60, (61, 61, 37)),
]


@pytest.mark.parametrize(("pdb", "resolution", "dim"), _FC_CASES)
def test_write_map_round_trips_fc_density_output(pdb, resolution, dim, tmp_path):
    """Section 6's fc_density case, on all three crystallographic assets.

    Distinct from the read-fixture round trips above in what produced the
    carrier: fc_density built this grid and computed its payload, where those
    cases write back a grid read_map constructed from a header. The geometry is
    inherited from the observed map, which is what keeps it writable -- it is
    the contrast case for the sub-box refusal below, since fc_density leaves the
    observed map's cell-to-spacing ratio alone and wrap_and_pad_grid does not.
    It is also the shape a caller most wants to write: a model map they just
    computed.
    """
    mol = oechem.OEGraphMol()
    ifs = oechem.oemolistream(str(_ASSET_DIR / f"{pdb}.cif"))
    assert oechem.OEReadMolecule(ifs, mol)
    ifs.close()

    source = read_map(_ASSET_DIR / f"{pdb}_2fofc.ccp4")
    cell = get_unit_cell(source.grid)
    # 90-degree angles, as test_validation.py:395 does for the same three
    # structures. DensityCalculator refuses any other angle outright
    # (src/DensityCalculator.cpp:43), so this is not a rounding convenience.
    fc = fc_density(
        mol, source.grid, resolution,
        UnitCell(cell.a, cell.b, cell.c, 90.0, 90.0, 90.0),
        symops=parse_symops(source.symops) or None,
    )

    out = tmp_path / f"{pdb}_fc.ccp4"
    maptitude.write_map(out, fc, source.symops)

    back = read_map(out)
    assert (back.grid.GetXDim(), back.grid.GetYDim(),
            back.grid.GetZDim()) == dim
    assert back.grid.GetSpaceGroup() == fc.GetSpaceGroup()
    assert back.symops == source.symops

    before = get_grid_params(fc)
    after = get_grid_params(back.grid)
    assert after.x_spacing == pytest.approx(before.x_spacing, rel=1e-6)
    assert after.y_spacing == pytest.approx(before.y_spacing, rel=1e-6)
    assert after.z_spacing == pytest.approx(before.z_spacing, rel=1e-6)
    assert after.x_origin == pytest.approx(before.x_origin, abs=1e-6)

    original = fc.GetValues()
    written = back.grid.GetValues()
    assert [written[i] for i in range(back.grid.GetSize())] == [
        original[i] for i in range(fc.GetSize())
    ]


def _padded_sub_box(source):
    """The wrap_and_pad_grid box around 1d26's ligand: the one refused shape.

    Lives in Python because the C++ suite has no molecule-from-file loader.
    Its C++ counterpart reaches the same branch with a hand-built grid, which
    reproduces the geometric condition but is not the call section 2.4 measured
    the refusal on -- so this is where that measurement is actually pinned.
    """
    mol = oechem.OEGraphMol()
    ifs = oechem.oemolistream(str(_ASSET_DIR / "1d26.cif"))
    assert oechem.OEReadMolecule(ifs, mol)
    ifs.close()

    cell = get_unit_cell(source.grid)
    padded = wrap_and_pad_grid(source.grid, mol, cell.a, cell.b, cell.c)
    assert padded is not source.grid, (
        "wrap_and_pad_grid returned the grid it was given, so the ligand already "
        "fits and there is no sub-box to refuse"
    )
    return padded


def test_write_map_refuses_a_wrap_and_pad_sub_box(tmp_path):
    source = read_map(_ASSET_DIR / "1d26_2fofc.ccp4")
    out = tmp_path / "sub_box.ccp4"

    with pytest.raises(GridError) as excinfo:
        maptitude.write_map(out, _padded_sub_box(source))

    assert "dimensions differ" in str(excinfo.value)
    assert not out.exists()


def test_a_refused_write_leaves_the_destination_alone(tmp_path):
    source = read_map(_ASSET_DIR / "1d26_2fofc.ccp4")
    out = tmp_path / "existing.ccp4"
    maptitude.write_map(out, source.grid, source.symops)
    before = out.read_bytes()

    with pytest.raises(GridError):
        maptitude.write_map(out, _padded_sub_box(source))

    assert out.read_bytes() == before


def test_write_map_accepts_a_string_path(tmp_path):
    # Every other case that writes successfully passes a pathlib.Path, so str
    # is the spelling nothing else covers.
    source = read_map(_DATA_DIR / "test_map.ccp4")
    out = tmp_path / "from_str.ccp4"
    maptitude.write_map(str(out), source.grid)
    assert read_map(out).grid.GetSize() == 9261


def test_write_map_refuses_a_path_with_an_embedded_nul(tmp_path):
    # os.fspath passes an embedded NUL through unchanged, and the std::string
    # typemap carries it across with its length, so the whole name reaches C++.
    # There a NUL splits it in two: std::filesystem reports the extension of
    # "victim.dat\0.ccp4" as ".ccp4", which the extension gate admits, while
    # every filesystem call downstream goes through c_str() and stops at the
    # NUL. Measured before the guard, this returned successfully having
    # overwritten victim.dat with a 38068-byte map, so the bytes are what this
    # pins; the refusal is caught rather than required, so that check runs
    # either way.
    source = read_map(_DATA_DIR / "test_map.ccp4")
    victim = tmp_path / "victim.dat"
    victim.write_bytes(b"ORIGINAL CONTENTS\n")

    refusal = None
    try:
        maptitude.write_map(str(victim) + "\0.ccp4", source.grid)
    except GridError as error:
        refusal = error

    assert victim.read_bytes() == b"ORIGINAL CONTENTS\n"
    assert refusal is not None


@pytest.mark.parametrize(
    "bad", [None, 3, ["m.ccp4"], b"tests/data/test_map.ccp4"])
def test_write_map_rejects_arguments_that_are_not_paths(bad, tmp_path,
                                                        monkeypatch):
    # Task 3's read_map case in write form, where a str(path) body costs more
    # than a misleading error: str(["m.ccp4"]) is "['m.ccp4']", whose extension
    # std::filesystem reports as ".ccp4']", which the extension gate admits, so
    # under a str(path) body that argument writes a real map under a name the
    # caller never gave. The chdir puts that name inside tmp_path, which is what
    # gives the directory check a place to look, and the refusal is caught
    # rather than wrapped in pytest.raises so the directory check runs even for
    # the list argument, which is the one a str(path) body does not refuse.
    # os.fspath rejects the first three; the std::string typemap rejects the
    # bytes, which os.fspath passes through unchanged.
    source = read_map(_DATA_DIR / "test_map.ccp4")
    monkeypatch.chdir(tmp_path)

    refusal = None
    try:
        maptitude.write_map(bad, source.grid)
    except TypeError as error:
        refusal = error

    assert list(tmp_path.iterdir()) == []
    assert refusal is not None


def test_read_map_scores_the_same_as_the_retired_loader():
    """RSCC, RSR and Q on one asset must not move when the loader changes.

    First pinned against values computed through the loader being retired, so
    a difference here meant the swap had moved a score rather than a call
    site; the swap moved nothing. RSCC and RSR were re-pinned once, when
    fc_density gained its own radius preparation: on a file-fresh molecule
    they had been 1.9e-4 and 4.5e-4 lower, the solvent mask's 1.7 A fallback
    standing in for every atom's radius, and they now equal what a
    scorer-prepared molecule always gave. Q does not use fc_density and did
    not move.

    The tolerance is the floor the C++ characterisation pins use for a pinned
    value below 1 (grid_summary.h: absolute 1e-6), not a same-build tolerance.
    The values were measured on macOS arm64; the Windows x64 and Linux aarch64
    wheel jobs computed RSCC 6.0e-8 and 3.7e-8 below the earlier pin from the
    same source, the cross-build drift that floor exists for.
    """
    mol = oechem.OEGraphMol()
    ifs = oechem.oemolistream(str(_ASSET_DIR / "340d.cif"))
    assert oechem.OEReadMolecule(ifs, mol)
    ifs.close()

    result = read_map(_ASSET_DIR / "340d_2fofc.ccp4")
    cell = get_unit_cell(result.grid)
    symops = parse_symops(result.symops) if result.symops else None
    grid = maptitude.wrap_and_pad_grid(result.grid, mol, cell.a, cell.b, cell.c)

    calc = maptitude.fc_density(
        mol, grid, 1.60,
        maptitude.UnitCell(cell.a, cell.b, cell.c, 90.0, 90.0, 90.0),
        symops=symops)

    assert maptitude.rscc(mol, grid, 1.60, calc_grid=calc).overall == (
        pytest.approx(_EXPECTED_RSCC, abs=1e-6))
    assert maptitude.rsr(mol, grid, 1.60, calc_grid=calc).overall == (
        pytest.approx(_EXPECTED_RSR, abs=1e-6))
    assert maptitude.qscore(mol, grid, 1.60).overall == (
        pytest.approx(_EXPECTED_QSCORE, abs=1e-6))


def test_combine_maps_rejects_an_out_of_range_op():
    """An out-of-range op used to return an all-zero grid with no error."""
    source = read_map(_DATA_DIR / "test_map.ccp4")
    other = read_map(_DATA_DIR / "test_map.ccp4")
    with pytest.raises(ValueError, match="MapOp value"):
        maptitude.combine_maps(source.grid, other.grid, 99)


def test_combine_maps_rejects_a_non_integer_op():
    source = read_map(_DATA_DIR / "test_map.ccp4")
    other = read_map(_DATA_DIR / "test_map.ccp4")
    with pytest.raises(TypeError, match="op must be one of"):
        maptitude.combine_maps(source.grid, other.grid, "ADD")


def test_combine_maps_rejects_a_bool_op():
    """True would otherwise be SUBTRACT and False would be ADD, silently."""
    source = read_map(_DATA_DIR / "test_map.ccp4")
    for value in (True, False):
        with pytest.raises(TypeError, match="op must be one of"):
            maptitude.combine_maps(source.grid, source.grid, value)


def test_combine_maps_still_accepts_every_valid_op():
    source = read_map(_DATA_DIR / "test_map.ccp4")
    other = read_map(_DATA_DIR / "test_map.ccp4")
    for op in (maptitude.MapOp.ADD, maptitude.MapOp.SUBTRACT,
               maptitude.MapOp.MIN, maptitude.MapOp.MAX):
        assert maptitude.combine_maps(source.grid, other.grid, op) is not None


def test_read_map_rejects_an_out_of_range_tiebreak():
    with pytest.raises(ValueError, match="OriginSource value"):
        read_map(_DATA_DIR / "test_map.ccp4", 99)


def test_read_map_still_accepts_both_tiebreaks():
    for tiebreak in (OriginSource.ORIGIN_RECORD, OriginSource.NXSTART):
        assert read_map(_DATA_DIR / "test_map.ccp4", tiebreak).grid is not None


def test_read_map_rejects_a_non_integer_tiebreak():
    with pytest.raises(TypeError, match="tiebreak must be"):
        read_map(_DATA_DIR / "test_map.ccp4", "ORIGIN_RECORD")


def test_read_map_rejects_a_bool_tiebreak():
    for value in (True, False):
        with pytest.raises(TypeError, match="tiebreak must be"):
            read_map(_DATA_DIR / "test_map.ccp4", value)


def test_the_enum_boundary_is_the_width_of_a_c_long():
    """The docstrings promise ValueError inside a C long and OverflowError past it.

    The limits are derived rather than written out: a C ``long`` is 64-bit on LP64
    platforms and 32-bit on Windows, and the documented claim is about the type, not
    about either width.
    """
    long_max = 2 ** (8 * ctypes.sizeof(ctypes.c_long) - 1) - 1
    long_min = -long_max - 1
    source = read_map(_DATA_DIR / "test_map.ccp4")
    for value in (long_max, long_min):
        with pytest.raises(ValueError, match="MapOp value"):
            maptitude.combine_maps(source.grid, source.grid, value)
        with pytest.raises(ValueError, match="OriginSource value"):
            read_map(_DATA_DIR / "test_map.ccp4", value)
    for value in (long_max + 1, long_min - 1):
        with pytest.raises(OverflowError):
            maptitude.combine_maps(source.grid, source.grid, value)
        with pytest.raises(OverflowError):
            read_map(_DATA_DIR / "test_map.ccp4", value)


def test_the_enum_boundary_rejects_minus_one_as_a_value_not_a_sentinel():
    """-1 is what PyLong_AsLong returns on failure, so it is the one colliding input.

    The typemaps separate the two meanings with ``&& PyErr_Occurred()``. Dropping that
    conjunct would send a real -1 to SWIG_fail with no exception set; this test is what
    notices. Its pair is the OverflowError assertion above, which is what notices if the
    sentinel check is dropped altogether.

    The patterns are anchored and carry the trailing text because ``match=`` is an
    unanchored ``re.search``: a bare ``MapOp value -1`` also matches the message for
    ``-10``, ``-100`` and every other value whose digits begin with ``1``.
    """
    source = read_map(_DATA_DIR / "test_map.ccp4")
    with pytest.raises(ValueError, match=r"^MapOp value -1 is out of range;"):
        maptitude.combine_maps(source.grid, source.grid, -1)
    with pytest.raises(ValueError, match=r"^OriginSource value -1 is out of range;"):
        read_map(_DATA_DIR / "test_map.ccp4", -1)


@pytest.mark.parametrize(
    "scalar",
    [np.int8(0), np.int32(0), np.int64(0), np.uint64(0), np.longlong(0)],
)
def test_the_enum_boundary_rejects_numpy_integer_scalars(scalar):
    """The docstrings promise the check is on the Python type, not on integer-ness."""
    source = read_map(_DATA_DIR / "test_map.ccp4")
    with pytest.raises(TypeError, match="op must be one of"):
        maptitude.combine_maps(source.grid, source.grid, scalar)
    with pytest.raises(TypeError, match="tiebreak must be"):
        read_map(_DATA_DIR / "test_map.ccp4", scalar)
