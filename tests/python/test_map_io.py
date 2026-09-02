"""Python-side tests for the CCP4/MRC map I/O bindings.

The C++ tests pin the values; these check that the SWIG boundary hands them
across intact -- that the grid arrives as a native ``openeye.oegrid.OESkewGrid``
rather than a wrapper, that the symmetry text survives as ``str``, and that the
typed exceptions still land as their Python classes.
"""

from __future__ import annotations

import struct
from pathlib import Path

import maptitude
import pytest
from maptitude import (
    CellError,
    GridError,
    MapFile,
    OriginSource,
    SymOpError,
    get_grid_params,
    get_unit_cell,
    parse_symops,
    read_map,
)
from openeye import oegrid

_ASSET_DIR = Path(__file__).resolve().parents[1] / "assets" / "mapq"
_DATA_DIR = Path(__file__).resolve().parents[1] / "data"

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
