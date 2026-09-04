"""Shared utilities for maptitude vs bms-bio benchmarks."""

import pathlib
import time

import maptitude
from openeye import oechem

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


def wrap_and_pad(grid, mol, cell, padding: float = 3.0):
    """Translate molecule into the unit cell and pad the grid if needed.

    An adapter over :func:`maptitude.wrap_and_pad_grid`, which takes the three
    cell edges separately where the loaders here return them as a tuple. The
    benchmarks measure the shipped padding path, so this must not grow a second
    implementation of it: an earlier Python copy that filled the padded grid by
    sampling the non-periodic entry point produced the better wwPDB agreement
    figures this suite once reported.

    The molecule is modified in-place (coordinates shifted).

    :param grid: Input OESkewGrid (one unit cell).
    :param mol: Molecule to shift.
    :param cell: Tuple (a, b, c) of cell dimensions.
    :param padding: Padding in Angstroms around the molecule.
    :returns: Padded grid, or the original when no padding is needed.
    :raises StructureError: If the molecule has no heavy atoms.
    :raises CellError: If an edge is zero, negative, or non-finite, or does not
        round to ``n`` or ``n - 1`` of that axis's node intervals to within the
        allowance made for float node coordinates.
    """
    a, b, c = cell
    return maptitude.wrap_and_pad_grid(grid, mol, a, b, c, padding)


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
