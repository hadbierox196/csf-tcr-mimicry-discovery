"""Builds the neuronal self-peptide / HLA-ligandome reference set.

STUB. This module should compile a reference of HLA class I/II ligands
from neuronal-tissue-expressed proteins, sourced from a combination of
in-house immunopeptidomics MS data (if available) and predicted ligands
for neuronal-expressed transcripts (e.g. filtered to the Human Protein
Atlas "neuronal" expression tier), for
:mod:`mimicry_discovery.self_antigen.query` to search against.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class SelfPeptide:
    """A candidate neuronal self-peptide presented on a given HLA allele.

    Attributes:
        peptide: Self-peptide sequence.
        hla_allele: Presenting HLA allele, two-field IPD-IMGT/HLA form.
        source_gene: Gene of origin.
        evidence: Provenance of this call, e.g. ``"ms_confirmed"`` or
            ``"predicted"``.
    """

    peptide: str
    hla_allele: str
    source_gene: str
    evidence: str


def build_neuronal_self_peptidome(
    neuronal_expression_reference: Path | str,
    output_path: Path | str,
) -> list[SelfPeptide]:
    """Compile the neuronal self-peptide/HLA reference set.

    Args:
        neuronal_expression_reference: Path to a neuronal-tissue gene
            expression reference (e.g. an HPA export) used to restrict
            candidate source proteins to those actually expressed in
            CNS tissue.
        output_path: Where to write the compiled reference.

    Returns:
        The compiled list of :class:`SelfPeptide` entries.

    Raises:
        NotImplementedError: Always -- TODO(self-antigen): source and
            wire in either in-house immunopeptidomics MS data, a public
            immunopeptidomics resource, or an HLA-binding-prediction
            fallback over neuronal-expressed transcripts.
    """
    raise NotImplementedError("build_neuronal_self_peptidome is a stub -- see module docstring.")
