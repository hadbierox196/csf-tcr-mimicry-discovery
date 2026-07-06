"""Unit tests for mimicry_discovery.structure.interface_metrics."""

from __future__ import annotations

from pathlib import Path

import pytest

from mimicry_discovery.structure.interface_metrics import (
    interface_residue_count,
    mean_plddt,
    parse_pdb_atoms,
)


def _fmt_atom_line(
    serial: int,
    name: str,
    res_name: str,
    chain_id: str,
    res_seq: int,
    x: float,
    y: float,
    z: float,
    occupancy: float,
    b_factor: float,
    element: str,
) -> str:
    """Format a single fixed-width PDB ATOM record for test fixtures."""
    return (
        f"ATOM  {serial:>5} {name:<4}{'':1}{res_name:<3} {chain_id}{res_seq:>4}{'':1}   "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{occupancy:>6.2f}{b_factor:>6.2f}"
        f"{'':10}{element:>2}"
    )


@pytest.fixture
def mini_pdb(tmp_path: Path) -> Path:
    """Write a minimal two-chain synthetic PDB file for interface-metric tests.

    Chain A has two CA atoms; chain B has two CA atoms, one of which
    (B-res1) sits close to both chain-A atoms, and one (B-res2) sits far
    from everything.
    """
    lines = [
        _fmt_atom_line(1, "CA", "ALA", "A", 1, 10.0, 10.0, 10.0, 1.0, 90.0, "C"),
        _fmt_atom_line(2, "CA", "GLY", "A", 2, 11.5, 10.0, 10.0, 1.0, 85.0, "C"),
        _fmt_atom_line(3, "CA", "SER", "B", 1, 13.5, 10.0, 10.0, 1.0, 70.0, "C"),
        _fmt_atom_line(4, "CA", "VAL", "B", 2, 30.0, 30.0, 30.0, 1.0, 60.0, "C"),
    ]
    pdb_path = tmp_path / "mini.pdb"
    pdb_path.write_text("\n".join(lines) + "\nEND\n")
    return pdb_path


def test_parse_pdb_atoms_reads_expected_fields(mini_pdb: Path) -> None:
    """Serial, chain, residue number, and coordinates parse correctly."""
    atoms = parse_pdb_atoms(mini_pdb)
    assert len(atoms) == 4
    assert atoms[0].chain_id == "A"
    assert atoms[0].res_seq == 1
    assert atoms[0].coord == pytest.approx((10.0, 10.0, 10.0))
    assert atoms[2].chain_id == "B"


def test_parse_pdb_atoms_missing_file_raises(tmp_path: Path) -> None:
    """A missing PDB path raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        parse_pdb_atoms(tmp_path / "missing.pdb")


def test_interface_residue_count_within_cutoff(mini_pdb: Path) -> None:
    """Both chain-A residues are within the 5A cutoff of chain-B residue 1."""
    atoms = parse_pdb_atoms(mini_pdb)
    n_contacts = interface_residue_count(atoms, chain_a_id="A", chain_b_id="B")
    assert n_contacts == 2


def test_interface_residue_count_missing_chain_returns_zero(mini_pdb: Path) -> None:
    """A chain ID with no atoms yields zero contacts rather than an error."""
    atoms = parse_pdb_atoms(mini_pdb)
    assert interface_residue_count(atoms, chain_a_id="A", chain_b_id="Z") == 0


def test_mean_plddt_overall_and_per_chain(mini_pdb: Path) -> None:
    """Mean pLDDT is correctly averaged overall and when restricted to a chain."""
    atoms = parse_pdb_atoms(mini_pdb)

    overall = mean_plddt(atoms)
    assert overall == pytest.approx((90 + 85 + 70 + 60) / 4)

    chain_a_only = mean_plddt(atoms, chain_id="A")
    assert chain_a_only == pytest.approx((90 + 85) / 2)


def test_mean_plddt_no_matching_atoms_raises(mini_pdb: Path) -> None:
    """Requesting a chain with no atoms raises ValueError, not a silent NaN."""
    atoms = parse_pdb_atoms(mini_pdb)
    with pytest.raises(ValueError, match="No atoms found"):
        mean_plddt(atoms, chain_id="Z")
