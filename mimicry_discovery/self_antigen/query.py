"""Queries the neuronal self-peptide reference for HLA/motif matches.

Depends on :func:`mimicry_discovery.self_antigen.build_reference.
build_neuronal_self_peptidome` having been run first to produce a
reference list. The filtering logic itself needs no domain-specific
model call, so unlike its sibling module this one is fully implemented
rather than stubbed.
"""

from __future__ import annotations

from mimicry_discovery.self_antigen.build_reference import SelfPeptide


def find_self_peptides_for_allele(
    hla_allele: str, reference: list[SelfPeptide]
) -> list[SelfPeptide]:
    """Return self-peptides presented on a given HLA allele.

    Args:
        hla_allele: Two-field IPD-IMGT/HLA allele string to match.
        reference: The compiled self-peptidome reference, typically
            loaded from the output of
            :func:`~mimicry_discovery.self_antigen.build_reference.build_neuronal_self_peptidome`.

    Returns:
        The subset of ``reference`` presented on ``hla_allele``, in the
        same relative order.
    """
    return [entry for entry in reference if entry.hla_allele == hla_allele]


def find_self_peptides_by_length(
    peptide_length: int, reference: list[SelfPeptide]
) -> list[SelfPeptide]:
    """Return self-peptides of a given length.

    Useful for restricting candidates to the length of a specific tumor
    peptide before running
    :func:`mimicry_discovery.scoring.features.peptide_similarity`, which
    requires equal-length inputs.

    Args:
        peptide_length: Exact peptide length to match.
        reference: The compiled self-peptidome reference.

    Returns:
        The subset of ``reference`` whose peptide has length
        ``peptide_length``.

    Raises:
        ValueError: If ``peptide_length`` is not positive.
    """
    if peptide_length <= 0:
        raise ValueError("peptide_length must be positive.")
    return [entry for entry in reference if len(entry.peptide) == peptide_length]
