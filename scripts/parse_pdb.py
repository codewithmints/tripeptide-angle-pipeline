import argparse
import csv
import logging
import multiprocessing as mp
import os
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")

# ── Biopython (heavy; imported once per worker) ───────────────────────────────
try:
    from Bio import PDB
    from Bio.PDB import DSSP
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    warnings.filterwarnings("ignore", category=PDBConstructionWarning)
except ImportError:
    sys.exit("[ERROR] Biopython is not installed. Run: pip install biopython")

# ── Logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# ── Constants ─────────────────────────────────────────────────────────────────

# Standard three-letter → one-letter amino acid mapping
THREE_TO_ONE: dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# Left-residue size categories (from assignment specification)
AA_CATEGORY: dict[str, str] = {
    "K": "Large",        "R": "Bulky",
    "D": "Intermediate", "E": "Large",
    "G": "Tiny",         "A": "Tiny",
    "V": "Small",        "L": "Intermediate",
    "I": "Intermediate", "M": "Large",
    "P": "Small",        "S": "Small",
    "T": "Small",        "N": "Intermediate",
    "Q": "Large",        "C": "Small",
    "F": "Bulky",        "Y": "Bulky",
    "W": "Bulky",        "H": "Large",
}

# DSSP secondary-structure codes that represent helices
HELIX_CODES: frozenset[str] = frozenset({"H", "G", "I"})
#   H = alpha helix, G = 3-10 helix, I = pi helix

# Output CSV columns
CSV_COLUMNS = ["pdb_id", "left_res", "central_res", "right_res", "angle", "category"]

# ── Global config (set once in initialiser) ───────────────────────────────────
_DSSP_EXE: str = "mkdssp"


def _worker_init(dssp_exe: str) -> None:
    """Initialise per-process globals (avoids re-passing via args)."""
    global _DSSP_EXE
    _DSSP_EXE = dssp_exe


# ── Angle computation (professor-provided, verbatim logic) ────────────────────

def signed_angle_3d(
    v1: np.ndarray,
    v2: np.ndarray,
    axis: np.ndarray,
    degrees: bool = True,
) -> float:
    """
    Signed angle from v1 to v2 around *axis*.

    Parameters
    ----------
    v1, v2 : array-like, shape (3,)
        Vectors to measure the angle between.
    axis : array-like, shape (3,)
        Reference axis that determines the sign.
    degrees : bool
        Return degrees (True) or radians (False).

    Returns
    -------
    float
        Signed angle in the range (–180, 180] degrees (or radians).
    """
    v1   = np.array(v1,   dtype=float)
    v2   = np.array(v2,   dtype=float)
    axis = np.array(axis, dtype=float)

    v1   /= np.linalg.norm(v1)
    v2   /= np.linalg.norm(v2)
    axis /= np.linalg.norm(axis)

    cross = np.cross(v1, v2)
    dot   = np.dot(v1, v2)
    angle = np.arctan2(np.dot(cross, axis), dot)

    return float(np.degrees(angle)) if degrees else float(angle)


# ── Structure helpers ─────────────────────────────────────────────────────────

def _ca_coords(residue: "PDB.Residue.Residue") -> np.ndarray | None:
    """Return the (3,) coordinates of the Cα atom, or None if absent."""
    if "CA" not in residue:
        return None
    return residue["CA"].get_vector().get_array().copy()


def _centroid(residue: "PDB.Residue.Residue") -> np.ndarray:
    """Return the (3,) centroid of *all* atoms in the residue."""
    coords = np.array(
        [a.get_vector().get_array() for a in residue.get_atoms()],
        dtype=float,
    )
    return coords.mean(axis=0)


def _run_dssp(model: "PDB.Model.Model", pdb_path: str) -> dict | None:
    """
    Attempt to run DSSP (mkdssp first, then 'dssp' as fallback).
    Returns a dict keyed by (chain_id, res_id) or None on failure.
    """
    for exe in (_DSSP_EXE, "dssp"):
        try:
            dssp_obj = DSSP(model, pdb_path, dssp=exe)
            return dict(dssp_obj)
        except Exception:
            continue
    return None


# ── Per-file processing ───────────────────────────────────────────────────────

def _process_single_pdb(pdb_path: str) -> list[dict]:
    """
    Parse one PDB file and return a list of angle records.
    Returns [] on any failure (file is skipped gracefully).
    """
    records: list[dict] = []
    parser  = PDB.PDBParser(QUIET=True)

    # ── 1. Parse structure ────────────────────────────────────────────────
    try:
        structure = parser.get_structure("prot", pdb_path)
    except Exception as exc:
        logger.warning("Parse error  [%s]: %s", os.path.basename(pdb_path), exc)
        return records

    pdb_id = os.path.splitext(os.path.basename(pdb_path))[0]

    # Only analyse the first MODEL (NMR ensembles: model 0 only)
    try:
        model = next(iter(structure))
    except StopIteration:
        return records

    # ── 2. DSSP ──────────────────────────────────────────────────────────
    dssp_dict = _run_dssp(model, pdb_path)
    if dssp_dict is None:
        logger.debug("DSSP failed   [%s]: skipping file", pdb_id)
        return records

    # ── 3. Iterate chains → tripeptides ──────────────────────────────────
    for chain in model:
        chain_id = chain.get_id()

        # Collect standard amino-acid residues only (skip HETATM / water)
        aa_residues = [
            r for r in chain.get_residues()
            if PDB.is_aa(r, standard=True)
        ]

        for i in range(1, len(aa_residues) - 1):
            left_res   = aa_residues[i - 1]
            center_res = aa_residues[i]
            right_res  = aa_residues[i + 1]

            # ── 3a. Central residue must be Arginine ──────────────────
            if center_res.get_resname() != "ARG":
                continue

            # ── 3b. Central residue must be in a helix (DSSP) ─────────
            dssp_key = (chain_id, center_res.get_id())
            if dssp_key not in dssp_dict:
                continue
            # dssp tuple: (index, aa, ss, rasa, phi, psi, ...)
            ss_code = dssp_dict[dssp_key][2]
            if ss_code not in HELIX_CODES:
                continue

            # ── 3c. Resolve one-letter codes ──────────────────────────
            left_one   = THREE_TO_ONE.get(left_res.get_resname())
            right_one  = THREE_TO_ONE.get(right_res.get_resname(), "X")

            if left_one is None or left_one not in AA_CATEGORY:
                continue  # unknown / non-standard left residue

            # ── 3d. Cα coordinates ────────────────────────────────────
            ca_left   = _ca_coords(left_res)
            ca_center = _ca_coords(center_res)

            if ca_left is None or ca_center is None:
                continue

            # ── 3e. CA→centroid vectors ───────────────────────────────
            v_left   = _centroid(left_res)   - ca_left
            v_center = _centroid(center_res) - ca_center

            if np.linalg.norm(v_left) < 1e-6 or np.linalg.norm(v_center) < 1e-6:
                continue  # degenerate (e.g., single-atom residue)

            # ── 3f. Axis for signed angle: CA_left → CA_center ────────
            axis = ca_center - ca_left
            if np.linalg.norm(axis) < 1e-6:
                continue  # coincident Cα atoms

            # ── 3g. Compute angle ─────────────────────────────────────
            angle = signed_angle_3d(v_left, v_center, axis, degrees=True)

            records.append(
                {
                    "pdb_id":     pdb_id,
                    "left_res":   left_one,
                    "central_res": "R",
                    "right_res":  right_one,
                    "angle":      round(angle, 4),
                    "category":   AA_CATEGORY[left_one],
                }
            )

    return records


def _safe_process(pdb_path: str) -> list[dict]:
    """Wrapper that catches *any* unexpected exception so the pool keeps running."""
    try:
        return _process_single_pdb(pdb_path)
    except Exception as exc:
        logger.warning("Unexpected error [%s]: %s", os.path.basename(pdb_path), exc)
        return []


# ── Main ──────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Parse a batch of PDB files and compute side-chain angles."
    )
    p.add_argument(
        "--manifest", required=True,
        help="Text file listing one PDB path per line.",
    )
    p.add_argument(
        "--output", required=True,
        help="Output CSV file path.",
    )
    p.add_argument(
        "--workers", type=int, default=4,
        help="Number of parallel worker processes (default: 4).",
    )
    p.add_argument(
        "--dssp", default="mkdssp",
        help="Name or path of the DSSP binary (default: mkdssp).",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    # ── Read manifest ─────────────────────────────────────────────────────
    if not os.path.isfile(args.manifest):
        sys.exit(f"[ERROR] Manifest not found: {args.manifest}")

    with open(args.manifest) as fh:
        pdb_paths = [ln.strip() for ln in fh if ln.strip()]

    if not pdb_paths:
        logger.warning("Empty manifest: %s – writing empty CSV.", args.manifest)
        _write_csv([], args.output)
        return

    logger.info(
        "Batch: %d PDB files | workers: %d | dssp: %s",
        len(pdb_paths), args.workers, args.dssp,
    )

    # ── Ensure output directory exists ────────────────────────────────────
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    # ── Parallel processing ───────────────────────────────────────────────
    all_records: list[dict] = []
    n_workers = max(1, args.workers)

    with mp.Pool(
        processes=n_workers,
        initializer=_worker_init,
        initargs=(args.dssp,),
    ) as pool:
        for file_records in pool.imap_unordered(
            _safe_process, pdb_paths, chunksize=10
        ):
            all_records.extend(file_records)

    logger.info("Extracted %d angle records from %d files.", len(all_records), len(pdb_paths))
    _write_csv(all_records, args.output)


def _write_csv(records: list[dict], path: str) -> None:
    """Write records (list of dicts) to *path* as a CSV."""
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(records)
    logger.info("Written → %s  (%d rows)", path, len(records))


if __name__ == "__main__":
    main()
