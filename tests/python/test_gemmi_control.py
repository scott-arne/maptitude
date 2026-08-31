"""Independent index-correspondence control for the map readers.

maptitude reads its grids through OpenEye. A check that compares the result
against OpenEye, or against values pinned from an earlier OpenEye read, cannot
see a misreading of CCP4/MRC that maptitude and OpenEye share. gemmi is an
unrelated implementation of the same formats, so putting the two decoded grids
side by side can.

What agreement here rules out, on the four assets in CASES:

* A format misreading maptitude and OpenEye share but gemmi does not, insofar
  as it shows in the quantities compared -- the dims, the unit cell, and the
  node values on the lattice the value test walks.
* An off-by-one at the closing plane. The dims relation is asserted between the
  two readers, and every position of all three appended planes is compared to
  the plane it copies.
* A node origin maptitude puts where gemmi's node 0 does not, to 1e-3 A, on
  390_emd_30342_A_z4.mrc -- the only asset the origin test reads.
* An x/y transposition of the grid *dims*, on 390_emd_30342_A_z4.mrc alone: at
  90 x 82 x 64 it is the only asset whose x and y dims differ.

What it does not rule out, and what a fifth asset should therefore target:

* An x/y transposition of the grid *spacings*. 1d26, 340d and 3q9g have
  x_spacing == y_spacing exactly (43.3/48, 43.18/60 and 32.867/36, the same on
  both axes); 390_emd_30342_A_z4.mrc has 99.05701/89 against 90.15301/81, which
  differ by 1.0e-08 relative -- a hundredfold under the rel=1e-6 the spacing
  test asserts at. So the two assertions that would catch such a swap compare
  quantities nothing here can tell apart, and the swap passes the whole module.
  An asset whose x and y spacings differ measurably is what closes this.
* A misplacement both readers share. The MRC ORIGIN record is exactly that
  case; test_neither_reader_applies_the_origin_record records it.
* Whatever hides at the payload nodes the value test skips. It walks a
  stride-3/5/7 lattice, roughly one node in a hundred.

tests/data/test_map.ccp4 is deliberately not covered here. Its header NX (21)
and MX (42) differ and its NxSTART is (-10, -10, -10): it is a sub-block of a
larger sampled cell, so the dims relation asserted below does not apply to it.
It is the isotropic-invariance fixture of the carrier phase, not a faithfulness
fixture.
"""

import pathlib

import gemmi
import pytest
from openeye import oegrid

import maptitude

_ASSET_DIR = pathlib.Path(__file__).parent.parent / "assets" / "mapq"

# (asset, gemmi dims, OpenEye dims, cell a/b/c/alpha/beta/gamma). OpenEye
# appends a closing plane on every axis, so its dims are gemmi's plus one on
# each; the appended plane is a byte copy of plane 0, the periodic image.
CASES = [
    ("1d26_2fofc.ccp4", (48, 48, 24), (49, 49, 25),
     (43.3, 43.3, 24.52, 90.0, 90.0, 90.0)),
    ("340d_2fofc.ccp4", (60, 60, 36), (61, 61, 37),
     (43.18, 43.18, 25.12, 90.0, 90.0, 90.0)),
    ("3q9g_2fofc.ccp4", (36, 36, 60), (37, 37, 61),
     (32.867, 32.867, 55.413, 90.0, 90.0, 90.0)),
    ("390_emd_30342_A_z4.mrc", (89, 81, 63), (90, 82, 64),
     (99.05701, 90.15301, 70.119, 90.0, 90.0, 90.0)),
]


def _read_both(name):
    path = str(_ASSET_DIR / name)
    g = gemmi.read_ccp4_map(path)
    oe = oegrid.OESkewGrid()
    assert oegrid.OEReadGrid(path, oe), f"OEReadGrid failed on {name}"
    return g, oe


@pytest.mark.parametrize("name,gemmi_dims,oe_dims,cell", CASES)
def test_dimensions_differ_by_exactly_the_closing_plane(name, gemmi_dims, oe_dims, cell):
    """KNOWN GAP, closed by the map I/O phase.

    OpenEye returns one node more per axis than the payload stores. Asserting
    the difference exactly, rather than asserting equality and tolerating a
    failure, makes this an executable statement of what the map I/O phase owes.
    That phase replaces the reader and this assertion is expected to change
    with it.
    """
    g, oe = _read_both(name)
    g_dims = (g.grid.nu, g.grid.nv, g.grid.nw)
    o_dims = (oe.GetXDim(), oe.GetYDim(), oe.GetZDim())
    assert g_dims == gemmi_dims
    assert o_dims == oe_dims

    # The relation has to hold between the two readers, not between two columns
    # of the table above, or deleting a reader would not change its outcome.
    assert o_dims == tuple(n + 1 for n in g_dims)

    # Grid.h:232 makes these three fields the basis of maptitude's element
    # index, so the node-for-node test below only means what it says if they
    # are the dims the reader actually returned.
    gp = maptitude.get_grid_params(oe)
    assert (gp.x_dim, gp.y_dim, gp.z_dim) == o_dims


@pytest.mark.parametrize("name,gemmi_dims,oe_dims,cell", CASES)
def test_per_axis_spacing_matches_gemmi(name, gemmi_dims, oe_dims, cell):
    """maptitude's walked spacing is gemmi's cell edge over its node count.

    A scalar carrier could not represent 1d26's (0.902083, 0.902083, 1.021667)
    without resampling an axis, which is why the carrier holds three spacings
    rather than one. This separates z from x and y on 1d26, 340d and 3q9g;
    390_emd_30342_A_z4.mrc holds its three spacings within 1.2e-07 relative of
    one another, a factor of nine under the rel=1e-6 asserted here, so it
    separates no pair of axes; and no asset here separates x from y (see the
    module docstring).
    """
    g, oe = _read_both(name)
    gp = maptitude.get_grid_params(oe)
    uc = g.grid.unit_cell

    assert gp.x_spacing == pytest.approx(uc.a / g.grid.nu, rel=1e-6)
    assert gp.y_spacing == pytest.approx(uc.b / g.grid.nv, rel=1e-6)
    assert gp.z_spacing == pytest.approx(uc.c / g.grid.nw, rel=1e-6)


@pytest.mark.parametrize("name,gemmi_dims,oe_dims,cell", CASES)
def test_unit_cell_matches_gemmi(name, gemmi_dims, oe_dims, cell):
    """The two readers resolve the same cell out of the same header words.

    The cell divided by the node count is the spacing, so a disagreement here
    stretches the node lattice. maptitude is therefore compared to gemmi
    directly: routing that comparison through the literals would bound the pair
    only by the sum of the two literal tolerances, which is far looser than the
    readers actually agree. The literals stay behind as an anchor that both
    readers opened the intended file.
    """
    g, oe = _read_both(name)
    uc = g.grid.unit_cell
    assert (uc.a, uc.b, uc.c) == pytest.approx(cell[:3], rel=1e-6)
    assert (uc.alpha, uc.beta, uc.gamma) == pytest.approx(cell[3:], rel=1e-6)

    mine = maptitude.get_unit_cell(oe)
    assert (mine.a, mine.b, mine.c) == pytest.approx(cell[:3], rel=1e-5)
    assert (mine.alpha, mine.beta, mine.gamma) == pytest.approx(cell[3:], rel=1e-5)

    # On all twelve lengths measured here maptitude returns the header's
    # float32 word widened to double and gemmi returns the decimal that word
    # rounds from, so the entire disagreement is that one round trip: 7.1e-09
    # to 4.7e-08 relative, against the 2**-24 = 6.0e-08 ceiling on a float32
    # rounding. 1e-6 clears that ceiling by more than an order of magnitude, so
    # a failure here is a real disagreement and not a rounding one. The angles
    # agree exactly, but every asset here is 90/90/90 and 90.0 is exact in
    # float32; they get the same tolerance for the first asset that is not.
    assert (mine.a, mine.b, mine.c) == pytest.approx((uc.a, uc.b, uc.c), rel=1e-6)
    assert (mine.alpha, mine.beta, mine.gamma) == pytest.approx(
        (uc.alpha, uc.beta, uc.gamma), rel=1e-6
    )


@pytest.mark.parametrize("name,gemmi_dims,oe_dims,cell", CASES)
def test_values_agree_node_for_node(name, gemmi_dims, oe_dims, cell):
    """The two readers land the same density on the same element index.

    Dims and cell can both be right while the payload is transposed, shifted by
    a plane, or walked in the wrong axis order, which no comparison of dims or
    cell can show. This walks a strided lattice of the shared payload block and
    demands bit equality at every node it visits -- roughly one node in a
    hundred; the module docstring says what the rest leave open.
    """
    g, oe = _read_both(name)
    values = oe.GetValues()
    nx, ny, _ = oe_dims

    # Exact float equality: both readers decode the same IEEE-754 words, so any
    # difference is a decode or an index error, not rounding. The comparison
    # samples the shared payload block [0,NX) x [0,NY) x [0,NZ) and indexes
    # nowhere else -- OpenEye's appended planes have no gemmi counterpart. The
    # strides are literals so a failure is reproducible.
    for iw in range(0, gemmi_dims[2], 7):
        for iv in range(0, gemmi_dims[1], 5):
            for iu in range(0, gemmi_dims[0], 3):
                el = iw * nx * ny + iv * nx + iu
                assert values[el] == g.grid.get_value(iu, iv, iw), (
                    f"{name} disagrees at ({iu}, {iv}, {iw})"
                )


@pytest.mark.parametrize("name,gemmi_dims,oe_dims,cell", CASES)
def test_the_closing_plane_repeats_the_first_on_every_axis(name, gemmi_dims, oe_dims, cell):
    """A statement about OpenEye's reader convention, not about density.

    The appended planes have no gemmi counterpart, so they cannot be compared
    across readers; what can be checked is that each is the periodic image of
    plane 0 on its axis. The sweep is exhaustive -- every position of all three
    appended planes, some 37k across the four assets, which costs a few
    milliseconds.
    """
    _, oe = _read_both(name)
    values = oe.GetValues()
    nx, ny, nz = oe_dims

    for iz in range(nz):
        for iy in range(ny):
            base = iz * nx * ny + iy * nx
            assert values[base + nx - 1] == values[base], (
                f"{name} x closing plane differs at ({iy}, {iz})"
            )
    for iz in range(nz):
        for ix in range(nx):
            plane = iz * nx * ny
            assert values[plane + (ny - 1) * nx + ix] == values[plane + ix], (
                f"{name} y closing plane differs at ({ix}, {iz})"
            )
    for iy in range(ny):
        for ix in range(nx):
            row = iy * nx + ix
            assert values[(nz - 1) * nx * ny + row] == values[row], (
                f"{name} z closing plane differs at ({ix}, {iy})"
            )


def test_neither_reader_applies_the_origin_record():
    """KNOWN GAP, closed by the map I/O phase.

    390_emd_30342_A_z4.mrc is the only asset carrying a nonzero ORIGIN record
    (header words 50-52). Measured against gemmi 0.7.5 and OpenEye 2026.1.0,
    BOTH readers ignore it and place node 0 at the coordinate origin, so the
    two agree with each other while both sit below the map's true position by
    exactly ORIGIN. Agreement here is therefore not evidence of correctness --
    it is a shared gap, recorded so the map I/O phase has an executable
    statement of what it must change.
    """
    origin = (145.825, 112.825, 120.517)
    g, oe = _read_both("390_emd_30342_A_z4.mrc")

    assert [g.header_float(i) for i in (50, 51, 52)] == pytest.approx(origin, abs=1e-3)

    gemmi_node0 = g.grid.get_position(0, 0, 0)
    assert (gemmi_node0.x, gemmi_node0.y, gemmi_node0.z) == pytest.approx(
        (0.0, 0.0, 0.0), abs=1e-6
    )

    gp = maptitude.get_grid_params(oe)
    assert (gp.x_origin, gp.y_origin, gp.z_origin) == pytest.approx(
        (gemmi_node0.x, gemmi_node0.y, gemmi_node0.z), abs=1e-3
    )
