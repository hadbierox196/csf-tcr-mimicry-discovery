"""Extraction of interface confidence and geometry metrics.

From a predicted TCR-pMHC complex structure (PDB format).
Implements a minimal, dependency-free PDB ATOM-record reader rather than
depending on BioPython for this narrow task. Follows the standard PDB
fixed-column convention. AlphaFold-family models (and PyMOL/ChimeraX
exports generally) write per-residue pLDDT into the standard PDB
B-factor column, which is what mean_plddt reads back out.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

_CONTACT_DISTANCE_CUTOFF = 5.0  # Angstroms, standard interface-contact cutoff


@dataclass(frozen=True)
class PDBAtom:
    """One parsed ATOM record from a PDB file."""

    serial: int
    name: str
    res_name: str
    chain_id: str
    res_seq: int
    coord: tuple[float, float, float]
    b_factor: float


def parse_pdb_atoms(pdb_path: Path | str) -> list[PDBAtom]:
    """Parse ATOM/HETATM records from a PDB file using fixed-width columns."""
    pdb_path = Path(pdb_path)
    if not pdb_path.exists():
        raise FileNotFoundError(f"Structure file not found: {pdb_path}")

    atoms: list[PDBAtom] = []
    with pdb_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atoms.append(
                PDBAtom(
                    serial=int(line[6:11]),
                    name=line[12:16].strip(),
                    res_name=line[17:20].strip(),
                    chain_id=line[21:22].strip(),
                    res_seq=int(line[22:26]),
                    coord=(float(line[30:38]), float(line[38:46]), float(line[46:54])),
                    b_factor=float(line[60:66]),
                )
            )
    return atoms


def _chain_ca_atoms(atoms: list[PDBAtom], chain_id: str) -> list[PDBAtom]:
    """Return CA atoms belonging to a given chain, in file order.

    Args:
        atoms: Parsed atoms, e.g. from :func:`parse_pdb_atoms`.
        chain_id: Chain identifier to filter to.

    Returns:
        The subset of ``atoms`` with matching ``chain_id`` and CA name.
    """
    return [a for a in atoms if a.chain_id == chain_id and a.name == "CA"]


def interface_residue_count(
    atoms: list[PDBAtom],
    chain_a_id: str,
    chain_b_id: str,
    distance_cutoff: float = _CONTACT_DISTANCE_CUTOFF,
) -> int:
    """Count CA atoms in chain_a within distance_cutoff of chain_b."""
    chain_a = _chain_ca_atoms(atoms, chain_a_id)
    chain_b = _chain_ca_atoms(atoms, chain_b_id)
    if not chain_a or not chain_b:
        return 0
    b_coords: NDArray[np.float64] = np.array([a.coord for a in chain_b])
    n_contacts = 0
    for atom in chain_a:
        distances = np.linalg.norm(b_coords - np.array(atom.coord), axis=1)
        if np.any(distances <= distance_cutoff):
            n_contacts += 1
    return n_contacts


def mean_plddt(atoms: list[PDBAtom], chain_id: str | None = None) -> float:
    """Compute mean per-atom pLDDT stored in the B-factor column."""
    selected = [a for a in atoms if chain_id is None or a.chain_id == chain_id]
    if not selected:
        raise ValueError(f"No atoms found{f' for chain {chain_id}' if chain_id else ''}.")
    return float(np.mean([a.b_factor for a in selected]))
