"""
Validation tests for maptitude scoring metrics against reference data.

Tests RSCC/RSR against published wwPDB validation pipeline values (Smart et al.,
Acta Cryst D, 2018) and Q-Score against mapq reference values (Pintilie et al.,
2020).  All tests use the maptitude C++ bindings, verifying that the C++
implementations produce results consistent with the reference data.
"""

import math
import pathlib
import struct

import pytest

np = pytest.importorskip("numpy", reason="numpy required for validation tests")

_ASSET_DIR = pathlib.Path(__file__).parent.parent / "assets" / "mapq"


def _has_openeye():
    """Check if OpenEye toolkits are available."""
    try:
        from openeye import oechem, oegrid
        return True
    except ImportError:
        return False


pytestmark = pytest.mark.skipif(
    not _has_openeye(),
    reason="OpenEye toolkits not available",
)


# ---------------------------------------------------------------------------
# Reference data
# ---------------------------------------------------------------------------

# Per-residue Q-scores from mapq (sigma=0.6, resolution=3.9 A)
# Source: 390_7cec_A100__Q__390_emd_30342_A_z4_All.txt
_MAPQ_RESIDUE_QSCORES = {
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

# Per-atom Q-scores from mapq annotated CIF
# Source: 390_7cec_A100__Q__390_emd_30342_A_z4.cif (_atom_site.Q_score column)
_MAPQ_ATOM_QSCORES = {
    (1, "N"): 0.373, (1, "CA"): 0.541, (1, "C"): 0.591, (1, "O"): 0.393,
    (1, "CB"): 0.418, (1, "CG"): 0.236, (1, "CD1"): -0.069, (1, "CD2"): 0.262,
    (1, "CE1"): 0.124, (1, "CE2"): 0.293, (1, "CZ"): 0.386,
    (2, "N"): 0.779, (2, "CA"): 0.719, (2, "C"): 0.517, (2, "O"): 0.486,
    (2, "CB"): 0.552, (2, "CG"): 0.481, (2, "OD1"): 0.456, (2, "ND2"): 0.261,
    (3, "N"): 0.350, (3, "CA"): 0.023, (3, "C"): 0.466, (3, "O"): 0.393,
    (3, "CB"): 0.496, (3, "CG"): 0.564, (3, "CD1"): 0.646, (3, "CD2"): 0.356,
    (4, "N"): 0.633, (4, "CA"): 0.796, (4, "C"): 0.734, (4, "O"): 0.330,
    (4, "CB"): 0.758, (4, "CG"): 0.515, (4, "OD1"): 0.046, (4, "OD2"): 0.265,
    (5, "N"): 0.598, (5, "CA"): 0.703, (5, "C"): 0.569, (5, "O"): 0.556,
    (5, "CB"): 0.686, (5, "OG1"): 0.313, (5, "CG2"): 0.537,
    (6, "N"): 0.505, (6, "CA"): 0.736, (6, "C"): 0.628, (6, "O"): 0.283,
    (6, "CB"): 0.771, (6, "CG"): 0.311, (6, "CD"): 0.438, (6, "NE"): 0.453,
    (6, "CZ"): 0.635, (6, "NH1"): 0.511, (6, "NH2"): 0.594,
    (7, "N"): 0.575, (7, "CA"): 0.651, (7, "C"): 0.679, (7, "O"): 0.564,
    (7, "CB"): 0.594, (7, "CG"): 0.416, (7, "CD"): 0.721, (7, "OE1"): 0.529,
    (7, "OE2"): 0.724,
    (8, "N"): 0.736, (8, "CA"): 0.648, (8, "C"): 0.527, (8, "O"): 0.485,
    (8, "CB"): 0.436, (8, "CG"): 0.478, (8, "OD1"): 0.528, (8, "OD2"): 0.677,
    (9, "N"): 0.461, (9, "CA"): 0.619, (9, "C"): 0.522, (9, "O"): 0.529,
    (9, "CB"): 0.369, (9, "CG"): 0.447, (9, "OD1"): 0.379, (9, "ND2"): 0.443,
    (10, "N"): 0.705, (10, "CA"): 0.768, (10, "C"): 0.622, (10, "O"): 0.357,
    (10, "CB"): 0.683, (10, "CG1"): 0.549, (10, "CG2"): 0.577,
    (11, "N"): 0.672, (11, "CA"): 0.759, (11, "C"): 0.775, (11, "O"): 0.350,
    (11, "CB"): 0.737, (11, "CG1"): 0.495, (11, "CG2"): 0.506, (11, "CD1"): 0.447,
    (12, "N"): 0.783, (12, "CA"): 0.811, (12, "C"): 0.380, (12, "O"): 0.336,
    (12, "CB"): 0.743, (12, "CG"): 0.710, (12, "CD"): 0.468,
}

# Per-ligand RSR/RSCC from wwPDB validation pipeline
# Source: Smart et al., Acta Cryst D, 2018 (ba5278sup4.csv.gz)
_PUBLISHED_LIGAND_METRICS = {
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

# Resolutions from PDB headers
_PDB_RESOLUTIONS = {
    "1d26": 2.12,
    "3q9g": 2.05,
    "340d": 1.60,
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_mrc_grid(mrc_path: pathlib.Path):
    """Load an MRC/CCP4 grid and fix the origin from the MRC header.

    OpenEye's OEReadGrid does not apply the ORIGIN field (words 50-52) from
    the MRC header, so the grid ends up at the wrong spatial location. This
    helper reads the origin from the binary header and shifts the grid
    mid-coordinates accordingly.
    """
    from openeye import oechem, oegrid

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


def _load_ccp4_grid(ccp4_path: pathlib.Path):
    """Load a CCP4 map from PDBe EDS including symmetry operators.

    :returns: Tuple of (grid, (a, b, c) cell dimensions, symops_text).
    """
    from openeye import oechem, oegrid
    from maptitude import parse_symops

    with open(ccp4_path, "rb") as f:
        header = f.read(1024)
    a, b, c = struct.unpack_from("3f", header, 40)
    nsymbt = struct.unpack_from("<i", header, 92)[0]

    symops_list = []
    if nsymbt > 0:
        with open(ccp4_path, "rb") as f:
            f.seek(1024)
            sym_bytes = f.read(nsymbt)
        symops_list = parse_symops(sym_bytes.decode("ascii", errors="ignore"))

    grid = oegrid.OEScalarGrid()
    ifs = oechem.oeifstream(str(ccp4_path))
    oegrid.OEReadGrid(ifs, grid, oegrid.OEGridFileType_CCP4)
    ifs.close()

    return grid, (a, b, c), symops_list


def _wrap_and_pad_grid(grid, mol, cell, padding: float = 3.0):
    """Translate molecule into the unit cell and pad the grid if needed.

    Crystallographic CCP4 maps cover one unit cell.  Deposited coordinates
    may extend slightly beyond the cell boundary.  This helper:

    1. Applies a rigid-body shift (integer multiples of cell vectors) to
       bring the molecular centroid near the grid center.
    2. If any atoms still fall outside the grid (plus *padding*), creates a
       new grid covering the full coordinate range and fills it by sampling
       the original grid with periodic wrapping.

    The molecule is modified in-place (coordinates shifted).
    """
    from openeye import oechem, oegrid

    a, b, c = cell
    sp = grid.GetSpacing()
    coords = oechem.OEFloatArray(3)

    # Compute heavy-atom centroid
    cx, cy, cz, n = 0.0, 0.0, 0.0, 0
    for atom in mol.GetAtoms(oechem.OEIsHeavy()):
        mol.GetCoords(atom, coords)
        cx += coords[0]
        cy += coords[1]
        cz += coords[2]
        n += 1
    cx /= n
    cy /= n
    cz /= n

    # Shift centroid to grid centre using integer unit-cell vectors
    grid_xmid = grid.GetXMin() + (grid.GetXDim() - 1) * sp / 2.0
    grid_ymid = grid.GetYMin() + (grid.GetYDim() - 1) * sp / 2.0
    grid_zmid = grid.GetZMin() + (grid.GetZDim() - 1) * sp / 2.0

    shift_x = round((grid_xmid - cx) / a) * a if a > 0 else 0.0
    shift_y = round((grid_ymid - cy) / b) * b if b > 0 else 0.0
    shift_z = round((grid_zmid - cz) / c) * c if c > 0 else 0.0

    if abs(shift_x) > 0.01 or abs(shift_y) > 0.01 or abs(shift_z) > 0.01:
        for atom in mol.GetAtoms():
            mol.GetCoords(atom, coords)
            mol.SetCoords(
                atom,
                oechem.OEFloatArray(
                    [coords[0] + shift_x, coords[1] + shift_y, coords[2] + shift_z]
                ),
            )

    # Check whether all atoms fall inside the grid (with padding)
    xs, ys, zs = [], [], []
    for atom in mol.GetAtoms(oechem.OEIsHeavy()):
        mol.GetCoords(atom, coords)
        xs.append(float(coords[0]))
        ys.append(float(coords[1]))
        zs.append(float(coords[2]))

    grid_xmin = grid.GetXMin()
    grid_ymin = grid.GetYMin()
    grid_zmin = grid.GetZMin()
    grid_xmax = grid_xmin + (grid.GetXDim() - 1) * sp
    grid_ymax = grid_ymin + (grid.GetYDim() - 1) * sp
    grid_zmax = grid_zmin + (grid.GetZDim() - 1) * sp

    needs_pad = (
        min(xs) - padding < grid_xmin
        or max(xs) + padding > grid_xmax
        or min(ys) - padding < grid_ymin
        or max(ys) + padding > grid_ymax
        or min(zs) - padding < grid_zmin
        or max(zs) + padding > grid_zmax
    )

    if not needs_pad:
        return grid

    # Build a padded grid using periodic wrapping of the unit-cell density
    minmax = oechem.OEDoubleArray(
        [
            min(xs) - padding, min(ys) - padding, min(zs) - padding,
            max(xs) + padding, max(ys) + padding, max(zs) + padding,
        ]
    )
    padded = oegrid.OEScalarGrid(minmax, sp)
    orig_xmin = grid.GetXMin()
    orig_ymin = grid.GetYMin()
    orig_zmin = grid.GetZMin()
    for i in range(padded.GetSize()):
        x, y, z = padded.ElementToSpatialCoord(i)
        wx = orig_xmin + ((x - orig_xmin) % a)
        wy = orig_ymin + ((y - orig_ymin) % b)
        wz = orig_zmin + ((z - orig_zmin) % c)
        padded.SetValue(
            i, float(oechem.OEFloatGridLinearInterpolate(grid, wx, wy, wz, 0.0))
        )
    return padded


class _HasResidueName:
    """Atom predicate matching a specific residue name (for mask building)."""

    def __init__(self, name: str):
        from openeye import oechem
        self._pred = oechem.OEUnaryAtomPred()
        self._name = name

    def __call__(self, atom):
        from openeye import oechem
        res = oechem.OEAtomGetResidue(atom)
        return res.GetName().strip() == self._name


def _ligand_mask(res_name: str, chain_id: str, res_num: int):
    """Build an atom predicate selecting a single ligand residue."""
    from openeye import oechem

    return oechem.OEAndAtom(
        oechem.OEAndAtom(
            oechem.OEHasChainID(chain_id),
            oechem.OEHasResidueNumber(res_num),
        ),
        oechem.OEHasResName(res_name),
    )


# =====================================================================
# Q-Score validation against mapq reference (7CEC / EMD-30342)
# =====================================================================

class TestQScoreMapqComparison:
    """Compare maptitude Q-score against reference values from mapq.

    Uses the 7CEC structure (12 residues, 100 heavy atoms, chain A) and the
    EMD-30342 density map at 3.9 A resolution.  Reference Q-scores were
    computed by the mapq package (Pintilie et al., 2020) with sigma=0.6.
    """

    _RESOLUTION = 3.9

    @classmethod
    def setup_class(cls):
        from openeye import oechem
        from maptitude import qscore

        cls.mol = oechem.OEGraphMol()
        ifs = oechem.oemolistream(str(_ASSET_DIR / "390_7cec_A100.cif"))
        oechem.OEReadMolecule(ifs, cls.mol)
        ifs.close()

        cls.grid = _load_mrc_grid(_ASSET_DIR / "390_emd_30342_A_z4.mrc")

        cls.result = qscore(cls.mol, cls.grid, resolution=cls._RESOLUTION)

    def test_loads_correct_atom_count(self):
        assert self.mol.NumAtoms() == 100

    def test_grid_covers_structure(self):
        from openeye import oechem

        coords = oechem.OEFloatArray(3)
        for atom in self.mol.GetAtoms():
            self.mol.GetCoords(atom, coords)
            assert self.grid.IsInGrid(coords[0], coords[1], coords[2]), (
                f"Atom {atom.GetIdx()} at ({coords[0]:.1f}, {coords[1]:.1f}, "
                f"{coords[2]:.1f}) is outside grid bounds"
            )

    def test_overall_qscore_not_nan(self):
        assert not math.isnan(self.result.overall)

    def test_overall_qscore_positive(self):
        assert self.result.overall > 0.0

    def test_correct_residue_count(self):
        assert len(self.result.by_residue) == 12

    def test_correct_atom_count(self):
        assert len(self.result.by_atom) == 100

    def test_per_residue_correlation_with_reference(self):
        """Per-residue Q-scores should strongly correlate with mapq reference.

        With the default options matching the mapq algorithm (fixed sigma=0.6,
        map normalization, point isolation), the per-residue correlation should
        be high.  Minor differences remain due to interpolation method and
        sphere-point algorithm details.
        """
        ours = []
        refs = []
        for res, score in self.result.by_residue.items():
            if res.number in _MAPQ_RESIDUE_QSCORES and not math.isnan(score):
                ours.append(score)
                refs.append(_MAPQ_RESIDUE_QSCORES[res.number][1])

        assert len(ours) > 0, "No residues matched reference data"
        corr = float(np.corrcoef(ours, refs)[0, 1])
        assert corr > 0.7, (
            f"Per-residue correlation with mapq reference is {corr:.3f} "
            f"(expected > 0.7)"
        )

    def test_per_atom_correlation_with_reference(self):
        """Per-atom Q-scores should strongly correlate with mapq reference."""
        from openeye import oechem

        ours = []
        refs = []
        for atom in self.mol.GetAtoms():
            res = oechem.OEAtomGetResidue(atom)
            key = (res.GetResidueNumber(), atom.GetName().strip())
            if key in _MAPQ_ATOM_QSCORES and atom.GetIdx() in self.result.by_atom:
                our_val = self.result.by_atom[atom.GetIdx()]
                if not math.isnan(our_val):
                    ours.append(our_val)
                    refs.append(_MAPQ_ATOM_QSCORES[key])

        assert len(ours) > 0, "No atoms matched reference data"
        corr = float(np.corrcoef(ours, refs)[0, 1])
        assert corr > 0.7, (
            f"Per-atom correlation with mapq reference is {corr:.3f} "
            f"(expected > 0.7)"
        )

    def test_overall_qscore_close_to_reference(self):
        """Overall Q-score should be close to the mapq reference value of ~0.52."""
        assert abs(self.result.overall - 0.52) < 0.15, (
            f"Overall Q-score={self.result.overall:.3f} "
            f"(expected ~0.52 +/- 0.15)"
        )


# =====================================================================
# Integration tests: all 4 metrics on real 7CEC structure
# =====================================================================

class TestIntegration:
    """Integration tests using a real structure with real EM density.

    Uses the 7CEC fragment (12 residues, 100 heavy atoms, chain A) and the
    EMD-30342 density map at 3.9 A resolution.
    """

    _RESOLUTION = 3.9

    @classmethod
    def setup_class(cls):
        from openeye import oechem
        from maptitude import rscc, rsr, qscore, ediam

        cls.mol = oechem.OEGraphMol()
        ifs = oechem.oemolistream(str(_ASSET_DIR / "390_7cec_A100.cif"))
        oechem.OEReadMolecule(ifs, cls.mol)
        ifs.close()

        cls.grid = _load_mrc_grid(_ASSET_DIR / "390_emd_30342_A_z4.mrc")

        cls.rscc_result = rscc(cls.mol, cls.grid, resolution=cls._RESOLUTION)
        cls.rsr_result = rsr(cls.mol, cls.grid, resolution=cls._RESOLUTION)
        cls.qscore_result = qscore(cls.mol, cls.grid, resolution=cls._RESOLUTION)
        cls.ediam_result = ediam(cls.mol, cls.grid, resolution=cls._RESOLUTION)

    def test_rscc_on_real_structure(self):
        assert not math.isnan(self.rscc_result.overall)
        assert len(self.rscc_result.by_residue) > 0
        assert len(self.rscc_result.by_atom) > 0

    def test_rsr_on_real_structure(self):
        assert not math.isnan(self.rsr_result.overall)
        assert self.rsr_result.overall >= 0.0

    def test_qscore_on_real_structure(self):
        assert not math.isnan(self.qscore_result.overall)

    def test_ediam_on_real_structure(self):
        assert not math.isnan(self.ediam_result.overall)
        assert self.ediam_result.overall >= 0.0
        assert self.ediam_result.overall <= 1.0

    def test_rscc_and_rsr_return_same_atom_count(self):
        """RSCC and RSR should identify the same set of atoms."""
        assert set(self.rscc_result.by_atom.keys()) == set(
            self.rsr_result.by_atom.keys()
        )

    def test_qscore_and_ediam_return_same_atom_count(self):
        """Q-score and EDIAm should cover the same atoms."""
        assert set(self.qscore_result.by_atom.keys()) == set(
            self.ediam_result.by_atom.keys()
        )


# =====================================================================
# RSCC / RSR benchmark against published wwPDB validation values
# =====================================================================

class TestRSCCRSRBenchmark:
    """Benchmark RSCC/RSR against published wwPDB validation values.

    Uses three small X-ray structures with 2Fo-Fc maps from PDBe EDS
    and structure-factor (Fc) model density computed via Cromer-Mann
    scattering factors + inverse FFT with flat bulk solvent correction
    (Jiang & Brunger, 1994).  Crystal symmetry operators from the CCP4
    extended header are applied so the Fc calculation covers the full
    unit cell.  Riding hydrogens at idealized positions are included in
    the Fc calculation (include_h=True).

    Published per-ligand values come from the wwPDB validation pipeline
    (Smart et al., Acta Cryst D, 2018) which uses REFMAC5 Fc maps.

    - RSCC values agree within ~0.10 of published because RSCC is
      insensitive to the overall density scale.
    - RSR values are close to published thanks to crystal symmetry
      expansion and the bulk solvent correction.
    """

    @classmethod
    def setup_class(cls):
        from openeye import oechem, oegrid
        from maptitude import (
            rscc, rsr, fc_density, UnitCell, parse_symops,
        )

        cls.mols = {}
        cls.grids = {}
        cls.cells = {}
        cls.symops = {}
        cls.fc_grids = {}
        cls.rscc_results = {}
        cls.rsr_results = {}
        cls.ligand_rscc = {}
        cls.ligand_rsr = {}

        for pdb in ["1d26", "3q9g", "340d"]:
            mol = oechem.OEGraphMol()
            ifs = oechem.oemolistream(str(_ASSET_DIR / f"{pdb}.cif"))
            oechem.OEReadMolecule(ifs, mol)
            ifs.close()

            grid, cell_dims, symops_list = _load_ccp4_grid(
                _ASSET_DIR / f"{pdb}_2fofc.ccp4"
            )
            grid = _wrap_and_pad_grid(grid, mol, cell_dims)

            cls.mols[pdb] = mol
            cls.grids[pdb] = grid
            cls.cells[pdb] = cell_dims
            cls.symops[pdb] = symops_list

            # Build UnitCell for fc_density — use 90-degree angles for
            # these orthorhombic-like structures
            a, b, c = cell_dims
            cell = UnitCell(a, b, c, 90.0, 90.0, 90.0)

            # Compute Fc model density with bulk solvent correction +
            # symmetry + riding hydrogens
            fc_grid = fc_density(
                mol, grid, _PDB_RESOLUTIONS[pdb], cell,
                symops=symops_list or None,
                include_h=True,
            )
            cls.fc_grids[pdb] = fc_grid

            res_val = _PDB_RESOLUTIONS[pdb]
            cls.rscc_results[pdb] = rscc(
                mol, grid, resolution=res_val, calc_grid=fc_grid,
            )
            cls.rsr_results[pdb] = rsr(
                mol, grid, resolution=res_val, calc_grid=fc_grid,
            )

            # Per-ligand results
            cls.ligand_rscc[pdb] = {}
            cls.ligand_rsr[pdb] = {}
            for entry in _PUBLISHED_LIGAND_METRICS[pdb]:
                mask = _ligand_mask(
                    entry["resname"], entry["chain"], entry["resnum"]
                )
                key = (entry["resname"], entry["resnum"])
                cls.ligand_rscc[pdb][key] = rscc(
                    mol, grid, resolution=res_val, mask=mask, calc_grid=fc_grid,
                )
                cls.ligand_rsr[pdb][key] = rsr(
                    mol, grid, resolution=res_val, mask=mask, calc_grid=fc_grid,
                )

    # -- Sanity bounds -----------------------------------------------

    def test_grid_covers_all_atoms(self):
        """All heavy atoms should be inside the grid for each structure."""
        from openeye import oechem

        coords = oechem.OEFloatArray(3)
        for pdb in ["1d26", "3q9g", "340d"]:
            mol = self.mols[pdb]
            grid = self.grids[pdb]
            for atom in mol.GetAtoms(oechem.OEIsHeavy()):
                mol.GetCoords(atom, coords)
                assert grid.IsInGrid(coords[0], coords[1], coords[2]), (
                    f"{pdb}: atom {atom.GetIdx()} at "
                    f"({coords[0]:.1f},{coords[1]:.1f},{coords[2]:.1f}) "
                    f"is outside grid"
                )

    def test_rscc_in_valid_range(self):
        """All RSCC values should be in [-1, 1]."""
        for pdb in ["1d26", "3q9g", "340d"]:
            result = self.rscc_results[pdb]
            assert -1.0 <= result.overall <= 1.0, (
                f"{pdb} overall RSCC={result.overall:.3f}"
            )
            for res, val in result.by_residue.items():
                if not math.isnan(val):
                    assert -1.0 <= val <= 1.0, f"{pdb} {res} RSCC={val:.3f}"

    def test_rsr_in_valid_range(self):
        """All RSR values should be in [0, 1]."""
        for pdb in ["1d26", "3q9g", "340d"]:
            result = self.rsr_results[pdb]
            assert 0.0 <= result.overall <= 1.0, (
                f"{pdb} overall RSR={result.overall:.3f}"
            )
            for res, val in result.by_residue.items():
                if not math.isnan(val):
                    assert 0.0 <= val <= 1.0, f"{pdb} {res} RSR={val:.3f}"

    def test_ligand_results_not_nan(self):
        """All ligand RSCC/RSR values should be finite (not NaN)."""
        for pdb, entries in _PUBLISHED_LIGAND_METRICS.items():
            for entry in entries:
                key = (entry["resname"], entry["resnum"])
                rscc_val = self.ligand_rscc[pdb][key].overall
                rsr_val = self.ligand_rsr[pdb][key].overall
                assert not math.isnan(rscc_val), (
                    f"{pdb} {entry['resname']}:{entry['resnum']} RSCC is NaN"
                )
                assert not math.isnan(rsr_val), (
                    f"{pdb} {entry['resname']}:{entry['resnum']} RSR is NaN"
                )

    def test_per_residue_coverage(self):
        """Each structure should have per-residue scores for most residues."""
        for pdb in ["1d26", "3q9g", "340d"]:
            n_rscc = sum(
                1 for v in self.rscc_results[pdb].by_residue.values()
                if not math.isnan(v)
            )
            assert n_rscc > 0, (
                f"{pdb} has no non-NaN per-residue RSCC values"
            )

    # -- Ordering ----------------------------------------------------

    def test_ligand_rscc_ordering_3q9g(self):
        """Our RSCC should preserve published ordering: 4BF > worst ORN."""
        rscc_4bf = self.ligand_rscc["3q9g"][("4BF", 5)].overall
        rscc_orn6 = self.ligand_rscc["3q9g"][("ORN", 6)].overall
        rscc_orn10 = self.ligand_rscc["3q9g"][("ORN", 10)].overall
        assert rscc_4bf > min(rscc_orn6, rscc_orn10), (
            f"4BF RSCC ({rscc_4bf:.3f}) should exceed worst ORN "
            f"({min(rscc_orn6, rscc_orn10):.3f})"
        )

    # -- Absolute proximity to published RSCC ------------------------

    def test_ligand_rscc_proximity(self):
        """Each ligand RSCC should be within 0.10 of the published value."""
        for pdb, entries in _PUBLISHED_LIGAND_METRICS.items():
            for entry in entries:
                key = (entry["resname"], entry["resnum"])
                our_val = self.ligand_rscc[pdb][key].overall
                pub_val = entry["rscc"]
                delta = abs(our_val - pub_val)
                assert delta <= 0.10, (
                    f"{pdb} {entry['resname']}:{entry['resnum']} "
                    f"RSCC={our_val:.3f} vs published={pub_val:.3f} "
                    f"(|delta|={delta:.3f}, tolerance=0.10)"
                )

    def test_overall_rscc_high(self):
        """Overall RSCC (all atoms) should exceed 0.85 for each structure."""
        for pdb in ["1d26", "3q9g", "340d"]:
            val = self.rscc_results[pdb].overall
            assert val > 0.85, (
                f"{pdb} overall RSCC={val:.3f} (expected > 0.85)"
            )

    # -- RSR functional check ----------------------------------------

    def test_ligand_rsr_bounded(self):
        """Each ligand RSR should be in (0, 0.20) with Fc model density."""
        for pdb, entries in _PUBLISHED_LIGAND_METRICS.items():
            for entry in entries:
                key = (entry["resname"], entry["resnum"])
                our_val = self.ligand_rsr[pdb][key].overall
                assert our_val > 0.0, (
                    f"{pdb} {entry['resname']}:{entry['resnum']} "
                    f"RSR={our_val:.3f} (expected > 0.0)"
                )
                assert our_val < 0.20, (
                    f"{pdb} {entry['resname']}:{entry['resnum']} "
                    f"RSR={our_val:.3f} (expected < 0.20)"
                )

    def test_overall_rsr_bounded(self):
        """Overall RSR should be in (0, 0.25) for each structure."""
        for pdb in ["1d26", "3q9g", "340d"]:
            val = self.rsr_results[pdb].overall
            assert 0.0 < val < 0.25, (
                f"{pdb} overall RSR={val:.3f} (expected in (0, 0.25))"
            )

    def test_ligand_rsr_proximity(self):
        """Each ligand RSR should be within 0.10 of the published value."""
        for pdb, entries in _PUBLISHED_LIGAND_METRICS.items():
            for entry in entries:
                key = (entry["resname"], entry["resnum"])
                our_val = self.ligand_rsr[pdb][key].overall
                pub_val = entry["rsr"]
                delta = abs(our_val - pub_val)
                assert delta <= 0.10, (
                    f"{pdb} {entry['resname']}:{entry['resnum']} "
                    f"RSR={our_val:.3f} vs published={pub_val:.3f} "
                    f"(|delta|={delta:.3f}, tolerance=0.10)"
                )
