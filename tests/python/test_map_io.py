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
    cell = get_unit_cell(read_map(_ASSET_DIR / "3q9g_2fofc.ccp4").grid)
    assert cell.a == pytest.approx(32.867, abs=1e-3)
    assert cell.c == pytest.approx(55.413, abs=1e-3)
    assert cell.alpha == pytest.approx(90.0, abs=1e-4)


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


def test_read_map_grid_survives_repeated_attribute_access():
    # The out-typemap releases the heap grid into a Python-owned object exactly
    # once, at call time. Two attribute reads return that same object; this
    # pins that it is live and usable rather than a corpse.
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
