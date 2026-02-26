"""Shared utilities for maptitude vs bms-bio benchmarks."""

import pathlib
import struct
import time

from openeye import oechem, oegrid

# ---------------------------------------------------------------------------
# Paths and constants
# ---------------------------------------------------------------------------

ASSET_DIR = pathlib.Path(__file__).resolve().parent.parent / "tests" / "assets" / "mapq"

RESOLUTIONS = {"1d26": 2.12, "3q9g": 2.05, "340d": 1.60, "7cec": 3.9}

# Published per-ligand RSCC/RSR from wwPDB validation pipeline
# Source: Smart et al., Acta Cryst D, 2018 (ba5278sup4.csv.gz)
PUBLISHED_RSCC_RSR = {
    "1d26": [
        {"resname": "G31", "chain": "A", "resnum": 5,
         "rsr": 0.062, "rscc": 0.987},
    ],
    "3q9g": [
        {"resname": "4BF", "chain": "A", "resnum": 5,
         "rsr": 0.077, "rscc": 0.981},
        {"resname": "ORN", "chain": "A", "resnum": 6,
         "rsr": 0.090, "rscc": 0.953},
        {"resname": "HAO", "chain": "A", "resnum": 7,
         "rsr": 0.079, "rscc": 0.955},
        {"resname": "ORN", "chain": "A", "resnum": 10,
         "rsr": 0.081, "rscc": 0.965},
    ],
    "340d": [
        {"resname": "5CM", "chain": "A", "resnum": 2,
         "rsr": 0.066, "rscc": 0.973},
        {"resname": "5CM", "chain": "A", "resnum": 4,
         "rsr": 0.060, "rscc": 0.968},
    ],
}

# Per-residue Q-scores from mapq (sigma=0.6, resolution=3.9 A)
# Source: 390_7cec_A100__Q__390_emd_30342_A_z4_All.txt
MAPQ_QSCORES = {
    1: ("PHE", 0.322624),
    2: ("ASN", 0.531169),
    3: ("LEU", 0.411691),
    4: ("ASP", 0.509571),
    5: ("THR", 0.566151),
    6: ("ARG", 0.533167),
    7: ("GLU", 0.605880),
    8: ("ASP", 0.564404),
    9: ("ASN", 0.471037),
    10: ("VAL", 0.608591),
    11: ("ILE", 0.592637),
    12: ("ARG", 0.604452),
}


# ---------------------------------------------------------------------------
# Loading helpers
# ---------------------------------------------------------------------------

def load_mol(path: pathlib.Path):
    """Load a CIF/PDB file into an OEGraphMol.

    :param path: Path to the structure file.
    :returns: Loaded molecule.
    """
    mol = oechem.OEGraphMol()
    ifs = oechem.oemolistream(str(path))
    oechem.OEReadMolecule(ifs, mol)
    ifs.close()
    return mol


def load_mrc_grid(mrc_path: pathlib.Path):
    """Load an MRC grid and fix the origin from header bytes 196-207.

    :param mrc_path: Path to the MRC file.
    :returns: OEScalarGrid with corrected origin.
    """
    with open(mrc_path, "rb") as f:
        header = f.read(1024)
    origin_x, origin_y, origin_z = struct.unpack_from("3f", header, 196)

    grid = oegrid.OEScalarGrid()
    ifs = oechem.oeifstream(str(mrc_path))
    oegrid.OEReadGrid(ifs, grid, oegrid.OEGridFileType_CCP4)
    ifs.close()

    sp = grid.GetSpacing()
    grid.SetXMid(origin_x + (grid.GetXDim() - 1) * sp / 2.0)
    grid.SetYMid(origin_y + (grid.GetYDim() - 1) * sp / 2.0)
    grid.SetZMid(origin_z + (grid.GetZDim() - 1) * sp / 2.0)
    return grid


def load_ccp4_grid(ccp4_path: pathlib.Path):
    """Load a CCP4 map with cell dimensions and symmetry operators.

    Splits the extended header into 80-char records (the CCP4 format for
    symmetry operators) before returning.

    :param ccp4_path: Path to the CCP4 map file.
    :returns: Tuple of (grid, (a, b, c) cell dimensions, symops_text).
    """
    with open(ccp4_path, "rb") as f:
        header = f.read(1024)
    a, b, c = struct.unpack_from("3f", header, 40)
    nsymbt = struct.unpack_from("<i", header, 92)[0]

    symops_text = ""
    if nsymbt > 0:
        with open(ccp4_path, "rb") as f:
            f.seek(1024)
            raw = f.read(nsymbt).decode("ascii", errors="ignore")
        # CCP4 stores symops as 80-char fixed-width records
        records = [raw[i:i + 80].strip() for i in range(0, len(raw), 80)]
        symops_text = "\n".join(r for r in records if r)

    grid = oegrid.OEScalarGrid()
    ifs = oechem.oeifstream(str(ccp4_path))
    oegrid.OEReadGrid(ifs, grid, oegrid.OEGridFileType_CCP4)
    ifs.close()

    return grid, (a, b, c), symops_text


def wrap_and_pad(grid, mol, cell, padding: float = 3.0):
    """Translate molecule into the unit cell and pad the grid if needed.

    The molecule is modified in-place (coordinates shifted).

    :param grid: Input OEScalarGrid (one unit cell).
    :param mol: Molecule to shift.
    :param cell: Tuple (a, b, c) of cell dimensions.
    :param padding: Padding in Angstroms around the molecule.
    :returns: Original or padded grid covering the molecule.
    """
    a, b, c = cell
    sp = grid.GetSpacing()
    coords = oechem.OEFloatArray(3)

    # Compute heavy-atom centroid
    cx, cy, cz, n = 0.0, 0.0, 0.0, 0
    for atom in mol.GetAtoms(oechem.OEIsHeavy()):
        mol.GetCoords(atom, coords)
        cx += coords[0]; cy += coords[1]; cz += coords[2]; n += 1
    if n == 0:
        return grid
    cx /= n; cy /= n; cz /= n

    # Shift centroid to grid centre using integer unit-cell vectors
    gxm = grid.GetXMin() + (grid.GetXDim() - 1) * sp / 2.0
    gym = grid.GetYMin() + (grid.GetYDim() - 1) * sp / 2.0
    gzm = grid.GetZMin() + (grid.GetZDim() - 1) * sp / 2.0
    sx = round((gxm - cx) / a) * a if a > 0 else 0.0
    sy = round((gym - cy) / b) * b if b > 0 else 0.0
    sz = round((gzm - cz) / c) * c if c > 0 else 0.0

    if abs(sx) > 0.01 or abs(sy) > 0.01 or abs(sz) > 0.01:
        for atom in mol.GetAtoms():
            mol.GetCoords(atom, coords)
            mol.SetCoords(atom, oechem.OEFloatArray(
                [coords[0] + sx, coords[1] + sy, coords[2] + sz]))

    # Check whether all atoms fall inside the grid (with padding)
    xs, ys, zs = [], [], []
    for atom in mol.GetAtoms(oechem.OEIsHeavy()):
        mol.GetCoords(atom, coords)
        xs.append(float(coords[0]))
        ys.append(float(coords[1]))
        zs.append(float(coords[2]))

    gxmin, gymin, gzmin = grid.GetXMin(), grid.GetYMin(), grid.GetZMin()
    gxmax = gxmin + (grid.GetXDim() - 1) * sp
    gymax = gymin + (grid.GetYDim() - 1) * sp
    gzmax = gzmin + (grid.GetZDim() - 1) * sp

    if (min(xs) - padding < gxmin or max(xs) + padding > gxmax
            or min(ys) - padding < gymin or max(ys) + padding > gymax
            or min(zs) - padding < gzmin or max(zs) + padding > gzmax):
        minmax = oechem.OEDoubleArray([
            min(xs) - padding, min(ys) - padding, min(zs) - padding,
            max(xs) + padding, max(ys) + padding, max(zs) + padding])
        padded = oegrid.OEScalarGrid(minmax, sp)
        for i in range(padded.GetSize()):
            x, y, z = padded.ElementToSpatialCoord(i)
            wx = gxmin + ((x - gxmin) % a)
            wy = gymin + ((y - gymin) % b)
            wz = gzmin + ((z - gzmin) % c)
            padded.SetValue(i, float(
                oechem.OEFloatGridLinearInterpolate(grid, wx, wy, wz, 0.0)))
        return padded
    return grid


# ---------------------------------------------------------------------------
# Atom predicates
# ---------------------------------------------------------------------------

class LigandPred(oechem.OEUnaryAtomPred):
    """Atom predicate matching a specific residue by name, chain, and number."""

    def __init__(self, res_name: str, chain_id: str, res_num: int):
        super().__init__()
        self._name = res_name
        self._chain = chain_id
        self._num = res_num

    def __call__(self, atom):
        res = oechem.OEAtomGetResidue(atom)
        return (res.GetName().strip() == self._name
                and res.GetChainID().strip() == self._chain
                and res.GetResidueNumber() == self._num)

    def CreateCopy(self):
        return LigandPred(self._name, self._chain, self._num)


def ligand_mask(res_name: str, chain_id: str, res_num: int):
    """Build an atom predicate selecting a single ligand residue.

    :param res_name: Residue name (e.g. "G31").
    :param chain_id: Chain identifier (e.g. "A").
    :param res_num: Residue number.
    :returns: LigandPred instance.
    """
    return LigandPred(res_name, chain_id, res_num)


# ---------------------------------------------------------------------------
# Timing utility
# ---------------------------------------------------------------------------

def bench(func, n: int = 5, warmup: int = 1):
    """Time a function, returning (min_time_seconds, last_result).

    :param func: Callable to benchmark.
    :param n: Number of timed iterations.
    :param warmup: Number of warmup iterations.
    :returns: Tuple of (min_elapsed_seconds, result).
    """
    result = None
    for _ in range(warmup):
        result = func()
    times = []
    for _ in range(n):
        t0 = time.perf_counter()
        result = func()
        times.append(time.perf_counter() - t0)
    return min(times), result
