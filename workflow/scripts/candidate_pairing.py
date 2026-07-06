"""Snakemake script: generate TCR x tumor-neoantigen x self-peptide candidates.

Same testable-``main()`` pattern as ``ingest_qc.py`` -- see that file's
docstring.

Self-peptide matching depends on
:func:`mimicry_discovery.self_antigen.build_reference.build_neuronal_self_peptidome`,
which is still a documented stub (roadmap Sprint 3). Rather than fail this
whole rule until that lands, this script checks whether a compiled
reference file exists yet: if not, it emits tumor-reactivity candidates
with zero self-peptide matches (an honest "no self-antigen reference
exists yet" result, not a crash) so the rest of the DAG stays runnable.
Once Sprint 3 lands and a reference file exists at
``<reference_dir>/neuronal_self_peptidome.json``, this script picks it up
automatically -- no rule change needed.
"""

from __future__ import annotations

import itertools
import json
import logging
from pathlib import Path

from mimicry_discovery.self_antigen.build_reference import SelfPeptide
from mimicry_discovery.self_antigen.query import (
    find_self_peptides_by_length,
    find_self_peptides_for_allele,
)

logger = logging.getLogger(__name__)


def _load_self_peptide_reference(reference_dir: str) -> list[SelfPeptide]:
    """Load the compiled neuronal self-peptidome reference, if it exists.

    Args:
        reference_dir: The configured reference data directory
            (``config.paths.reference_dir``).

    Returns:
        The compiled reference, or an empty list with a logged warning if
        ``build_neuronal_self_peptidome`` hasn't been run yet (roadmap
        Sprint 3).
    """
    reference_path = Path(reference_dir) / "neuronal_self_peptidome.json"
    if not reference_path.exists():
        logger.warning(
            "No neuronal self-peptidome reference found at %s -- "
            "self_antigen.build_reference is still a stub (see "
            "docs/roadmap.md Sprint 3). Emitting zero self-peptide matches.",
            reference_path,
        )
        return []
    raw = json.loads(reference_path.read_text())
    return [SelfPeptide(**entry) for entry in raw]


def main(
    clonotypes_path: str,
    neoantigens_path: str,
    reference_dir: str,
    candidates_out: str,
) -> None:
    """Generate TCR x neoantigen x self-peptide candidate triples.

    Args:
        clonotypes_path: Path to the ``clonotypes.json`` written by the
            ``ingest_qc`` rule.
        neoantigens_path: Path to the ``neoantigens.json`` written by the
            ``ingest_qc`` rule.
        reference_dir: Configured reference data directory, used to look
            for a compiled self-peptidome reference.
        candidates_out: Where to write the generated candidates (JSON).
    """
    clonotypes = json.loads(Path(clonotypes_path).read_text())
    neoantigens = json.loads(Path(neoantigens_path).read_text())
    self_peptide_reference = _load_self_peptide_reference(reference_dir)

    candidates = []
    for clonotype, neoantigen in itertools.product(clonotypes, neoantigens):
        if clonotype["cdr3_beta"] is None:
            # Every backend needs at least the beta chain; TCRdock further
            # requires alpha too (checked at the structure_prediction
            # stage, not here, so alpha-only-missing candidates still
            # reach scoring via a beta-only-capable backend like ESMFold).
            continue
        same_allele = find_self_peptides_for_allele(
            neoantigen["hla_allele"], self_peptide_reference
        )
        same_length = find_self_peptides_by_length(
            len(neoantigen["peptide"]), same_allele
        )
        for self_peptide in same_length:
            candidates.append(
                {
                    "sample_id": neoantigen["sample_id"],
                    "tcr_cdr3_beta": clonotype["cdr3_beta"],
                    "tcr_v_gene": clonotype["v_gene_beta"],
                    "tcr_j_gene": clonotype["j_gene_beta"],
                    "tcr_cdr3_alpha": clonotype["cdr3_alpha"],
                    "tcr_v_gene_alpha": clonotype["v_gene_alpha"],
                    "tcr_j_gene_alpha": clonotype["j_gene_alpha"],
                    "tumor_peptide": neoantigen["peptide"],
                    "self_peptide": self_peptide.peptide,
                    "hla_allele": neoantigen["hla_allele"],
                    "source_gene": self_peptide.source_gene,
                }
            )

    Path(candidates_out).parent.mkdir(parents=True, exist_ok=True)
    Path(candidates_out).write_text(json.dumps(candidates, indent=2))
    logger.info(
        "%s: %d clonotypes x %d neoantigens -> %d candidates (%d self-peptides in reference)",
        neoantigens[0]["sample_id"] if neoantigens else "unknown",
        len(clonotypes),
        len(neoantigens),
        len(candidates),
        len(self_peptide_reference),
    )


if "snakemake" in globals():
    logging.basicConfig(level=logging.INFO)
    main(
        clonotypes_path=snakemake.input.clonotypes,  # noqa: F821
        neoantigens_path=snakemake.input.neoantigens,  # noqa: F821
        reference_dir=snakemake.params.reference_dir,  # noqa: F821
        candidates_out=snakemake.output.candidates,  # noqa: F821
    )
