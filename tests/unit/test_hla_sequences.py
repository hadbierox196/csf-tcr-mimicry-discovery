"""Unit tests for mimicry_discovery.io.hla_sequences.

The header format used in the fixtures below was verified against a
real fetched IMGT/HLA FASTA file (see hla_sequences.py's module
docstring) -- the sequence content itself is a short synthetic stand-in,
not a real allele's actual sequence, since these tests only need to
exercise the parser.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mimicry_discovery.io.hla_sequences import (
    load_hla_protein_sequences,
    to_two_field_lookup,
)

_SYNTHETIC_FASTA = """\
>HLA:HLA00001 A*02:01:01:01 11 aa
MAVMAPRTLLL
>HLA:HLA00002 A*02:01:01:02 11 aa
MAVMAPRTLLL
>HLA:HLA00003 B*07:02:01:01 9 aa
GSHSMRYFY
"""


def test_load_hla_protein_sequences_parses_multi_record_fasta(tmp_path: Path) -> None:
    """Multiple records, each header/sequence pair parsed independently."""
    fasta_path = tmp_path / "test_prot.fasta"
    fasta_path.write_text(_SYNTHETIC_FASTA)

    sequences = load_hla_protein_sequences(fasta_path)

    assert sequences == {
        "A*02:01:01:01": "MAVMAPRTLLL",
        "A*02:01:01:02": "MAVMAPRTLLL",
        "B*07:02:01:01": "GSHSMRYFY",
    }


def test_load_hla_protein_sequences_missing_file_raises(tmp_path: Path) -> None:
    """A missing FASTA path raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        load_hla_protein_sequences(tmp_path / "missing.fasta")


def test_load_hla_protein_sequences_rejects_malformed_header(tmp_path: Path) -> None:
    """A header not matching the IMGT/HLA convention raises a clear ValueError."""
    bad_path = tmp_path / "bad.fasta"
    bad_path.write_text(">not_the_right_format\nMAVMAPRTLLL\n")
    with pytest.raises(ValueError, match="Unrecognized IMGT/HLA FASTA header"):
        load_hla_protein_sequences(bad_path)


def test_load_hla_protein_sequences_joins_multiline_sequences(tmp_path: Path) -> None:
    """A sequence wrapped across multiple lines is joined into one string."""
    fasta_path = tmp_path / "wrapped.fasta"
    fasta_path.write_text(">HLA:HLA00001 A*02:01:01:01 11 aa\nMAVM\nAPRT\nLLL\n")

    sequences = load_hla_protein_sequences(fasta_path)

    assert sequences["A*02:01:01:01"] == "MAVMAPRTLLL"


def test_to_two_field_lookup_collapses_and_keeps_first(tmp_path: Path) -> None:
    """Full-resolution alleles sharing a two-field prefix collapse, first wins."""
    fasta_path = tmp_path / "test_prot.fasta"
    fasta_path.write_text(_SYNTHETIC_FASTA)
    sequences = load_hla_protein_sequences(fasta_path)

    two_field = to_two_field_lookup(sequences)

    assert set(two_field.keys()) == {"A*02:01", "B*07:02"}
    assert two_field["A*02:01"] == "MAVMAPRTLLL"  # first (...01:01) wins
    assert two_field["B*07:02"] == "GSHSMRYFY"
