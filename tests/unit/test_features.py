"""Unit tests for mimicry_discovery.scoring.features."""

from __future__ import annotations

import pytest

from mimicry_discovery.scoring.features import anchor_conservation, peptide_similarity


def test_peptide_similarity_identical_sequences_is_one() -> None:
    """Comparing a peptide to itself yields the maximum similarity score."""
    assert peptide_similarity("SLLMWITQV", "SLLMWITQV") == pytest.approx(1.0)


def test_peptide_similarity_rejects_length_mismatch() -> None:
    """Peptides of different lengths are rejected rather than silently compared."""
    with pytest.raises(ValueError, match="equal-length"):
        peptide_similarity("SLLMWITQV", "SLLMW")


def test_peptide_similarity_ranks_conservative_above_nonconservative_substitution() -> None:
    """A same-property-class substitution scores higher than a cross-class one."""
    conservative = peptide_similarity("AAAAAAAAA", "VVVVVVVVV")  # both nonpolar
    nonconservative = peptide_similarity("AAAAAAAAA", "DDDDDDDDD")  # nonpolar vs. negative
    assert conservative > nonconservative


def test_peptide_similarity_property_match_credit_is_configurable() -> None:
    """A property_match_credit of 0 gives conservative substitutions no credit."""
    score = peptide_similarity("AAAAAAAAA", "VVVVVVVVV", property_match_credit=0.0)
    assert score == pytest.approx(0.0)


def test_anchor_conservation_all_positions_match() -> None:
    """A peptide pair sharing both anchor residues scores 1.0."""
    # Position 2 (index 1) and position 9 (index 8) match; the rest differ.
    score = anchor_conservation("SLLMWITQV", "XLXXXXXXV", anchor_positions=(2, 9))
    assert score == pytest.approx(1.0)


def test_anchor_conservation_no_positions_match() -> None:
    """A peptide pair matching at neither anchor position scores 0.0."""
    score = anchor_conservation("SLLMWITQV", "XXXXXXXXX", anchor_positions=(2, 9))
    assert score == pytest.approx(0.0)


def test_anchor_conservation_supports_negative_indexing() -> None:
    """A negative anchor position is interpreted relative to the C-terminus."""
    # position 2 (index 1) matches; position -1 (index 4, C-term) does not.
    score = anchor_conservation("ABCDE", "XBCDX", anchor_positions=(2, -1))
    assert score == pytest.approx(0.5)


def test_anchor_conservation_out_of_range_position_raises() -> None:
    """An anchor position beyond the peptide length raises ValueError."""
    with pytest.raises(ValueError, match="out of range"):
        anchor_conservation("SLL", "SLL", anchor_positions=(9,))
