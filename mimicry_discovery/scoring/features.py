"""Feature engineering for the molecular-mimicry risk score.

Combines peptide sequence similarity, HLA anchor-residue conservation,
and (elsewhere) structural confidence into the feature vector consumed
by :mod:`mimicry_discovery.scoring.risk_model`.

``peptide_similarity`` below uses a simple identity + physicochemical
property-group scheme rather than a BLOSUM-style substitution matrix.
TODO(scoring): swap in a real substitution matrix (e.g. BLOSUM62 via
Biopython's ``Bio.Align.substitution_matrices``, listed in
``environment.yml``) for biochemically-weighted scoring once that
dependency is available in your runtime -- the property-group scheme
here is a deliberately simple, dependency-free baseline in the
meantime, not a claim that it's the more rigorous option.
"""

from __future__ import annotations

# Standard four-way amino acid physicochemical classification.
_PROPERTY_GROUPS: dict[str, str] = {
    **{aa: "nonpolar" for aa in "GAVLIMPFW"},
    **{aa: "polar" for aa in "STCYNQ"},
    **{aa: "positive" for aa in "KRH"},
    **{aa: "negative" for aa in "DE"},
}


def _validate_equal_length(peptide_a: str, peptide_b: str, fn_name: str) -> None:
    """Raise ValueError if the two peptides differ in length.

    Args:
        peptide_a: First peptide sequence.
        peptide_b: Second peptide sequence.
        fn_name: Name of the calling function, used in the error
            message.

    Raises:
        ValueError: If the lengths differ.
    """
    if len(peptide_a) != len(peptide_b):
        raise ValueError(
            f"{fn_name} requires equal-length peptides, got "
            f"{len(peptide_a)} and {len(peptide_b)}."
        )


def peptide_similarity(
    peptide_a: str, peptide_b: str, property_match_credit: float = 0.5
) -> float:
    """Score similarity between two equal-length peptides.

    Used to compare a tumor-derived peptide against a candidate
    neuronal self-peptide presented on a compatible HLA allele -- the
    core "molecular mimicry" signal. Each position contributes 1.0 if
    the residues are identical, ``property_match_credit`` if they
    differ but share a physicochemical class (see
    :data:`_PROPERTY_GROUPS`), and 0.0 otherwise.

    Args:
        peptide_a: First peptide sequence.
        peptide_b: Second peptide sequence, same length as
            ``peptide_a``.
        property_match_credit: Partial credit for a conservative
            (same-property-class) substitution. Must be in ``[0, 1]``.

    Returns:
        A similarity score in ``[0, 1]``, where 1.0 means identical
        sequences.

    Raises:
        ValueError: If the two peptides differ in length, or
            ``property_match_credit`` is outside ``[0, 1]``.
    """
    _validate_equal_length(peptide_a, peptide_b, "peptide_similarity")
    if not 0.0 <= property_match_credit <= 1.0:
        raise ValueError("property_match_credit must be in [0, 1].")

    total = 0.0
    for a, b in zip(peptide_a.upper(), peptide_b.upper()):
        if a == b:
            total += 1.0
        elif _PROPERTY_GROUPS.get(a) == _PROPERTY_GROUPS.get(b) and a in _PROPERTY_GROUPS:
            total += property_match_credit
    return total / len(peptide_a)


def anchor_conservation(
    peptide_a: str,
    peptide_b: str,
    anchor_positions: tuple[int, ...] = (2, 9),
) -> float:
    """Fraction of HLA anchor positions that are identical between peptides.

    HLA class I binding is dominated by a small number of anchor
    positions (commonly position 2 and the C-terminus); peptides that
    diverge elsewhere but match at anchors are more likely to be
    presented in a similar conformation, which matters for whether the
    same TCR could plausibly engage both.

    Args:
        peptide_a: First peptide sequence.
        peptide_b: Second peptide sequence, same length as
            ``peptide_a``.
        anchor_positions: 1-indexed anchor positions to compare (a
            negative value counts from the C-terminus, e.g. ``-1`` is
            the last residue). Defaults to the canonical HLA class I
            P2/PΩ pair for a 9-mer. TODO(scoring): parameterize this
            per-allele using a real binding-motif reference rather than
            a fixed default.

    Returns:
        Fraction (0-1) of the given anchor positions that match
        exactly.

    Raises:
        ValueError: If the two peptides differ in length, or an anchor
            position falls outside the peptide.
    """
    _validate_equal_length(peptide_a, peptide_b, "anchor_conservation")
    length = len(peptide_a)
    matches = 0
    for pos in anchor_positions:
        index = pos - 1 if pos > 0 else length + pos
        if not 0 <= index < length:
            raise ValueError(
                f"Anchor position {pos} is out of range for a peptide of "
                f"length {length}."
            )
        if peptide_a[index] == peptide_b[index]:
            matches += 1
    return matches / len(anchor_positions)
