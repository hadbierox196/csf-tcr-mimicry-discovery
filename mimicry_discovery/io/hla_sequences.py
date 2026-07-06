"""Loads IPD-IMGT/HLA protein sequences from the official reference FASTA.

Source: https://github.com/ANHIG/IMGTHLA -- specifically the per-locus
``X_prot.fasta`` files (or the combined ``hla_prot.fasta``). The header
format below was verified against a real fetched file, not assumed from
memory (``raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/fasta/C_nuc.fasta``,
checked 2026-07)::

    >HLA:HLA00401 C*01:02:01:01 1101 bp
    ATGCGGGTCATGGCGCCCCGAACC...

That confirms the accession/allele-name/length header convention for the
nucleotide FASTA variant. The ANHIG/IMGTHLA README documents the same
convention for ``X_prot.fasta`` files; this environment couldn't fetch an
actual protein FASTA to independently confirm the length-unit suffix
reads ``"aa"`` there instead of ``"bp"`` (the parser below doesn't care
what that suffix says, only that a length token follows the allele name,
so this doesn't affect correctness -- noted for anyone auditing this
against the real file for the first time).
"""

from __future__ import annotations

import re
from pathlib import Path

_HEADER_PATTERN = re.compile(r"^>HLA:\S+\s+(?P<allele>\S+)\s+\d+\s*\S*")


def load_hla_protein_sequences(fasta_path: Path | str) -> dict[str, str]:
    """Parse an IPD-IMGT/HLA protein FASTA file into an allele->sequence map.

    Args:
        fasta_path: Path to ``hla_prot.fasta`` or a per-locus
            ``X_prot.fasta`` file downloaded from
            https://github.com/ANHIG/IMGTHLA (see ``scripts/download_reference_data.sh``).

    Returns:
        Dict mapping the full-resolution allele name exactly as it
        appears in the FASTA header (e.g. ``"A*02:01:01:01"`` -- without
        the ``"HLA-"`` prefix used elsewhere in this codebase) to amino
        acid sequence.

    Raises:
        FileNotFoundError: If ``fasta_path`` does not exist.
        ValueError: If a header line doesn't match the expected
            ``>HLA:<accession> <allele> <length> <unit>`` format.
    """
    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"HLA protein FASTA not found: {fasta_path}")

    sequences: dict[str, str] = {}
    current_allele: str | None = None
    current_seq: list[str] = []

    def _flush() -> None:
        """Commit the in-progress allele's accumulated sequence lines."""
        if current_allele is not None:
            sequences[current_allele] = "".join(current_seq)

    for line in fasta_path.read_text().splitlines():
        if line.startswith(">"):
            _flush()
            match = _HEADER_PATTERN.match(line)
            if match is None:
                raise ValueError(f"Unrecognized IMGT/HLA FASTA header: {line!r}")
            current_allele = match.group("allele")
            current_seq = []
        elif line.strip():
            current_seq.append(line.strip())
    _flush()
    return sequences


def to_two_field_lookup(sequences: dict[str, str]) -> dict[str, str]:
    """Collapse a full-resolution allele->sequence map to two-field keys.

    Most binding-relevant comparisons only need two-field resolution
    (e.g. ``"A*02:01"``), but the reference FASTA is keyed by
    full-resolution names (e.g. ``"A*02:01:01:01"``). Where multiple
    full-resolution alleles share a two-field prefix, the first one
    encountered wins -- a simplification worth knowing about (some
    alleles synonymous at two-field resolution do differ at higher
    resolution), documented here rather than silently applied.

    Args:
        sequences: Output of :func:`load_hla_protein_sequences`.

    Returns:
        Dict keyed by two-field allele name, e.g. ``"A*02:01"``.
    """
    two_field: dict[str, str] = {}
    for allele, sequence in sequences.items():
        fields = allele.split(":")
        key = ":".join(fields[:2])
        two_field.setdefault(key, sequence)
    return two_field
