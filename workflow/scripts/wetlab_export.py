"""Snakemake script: export the top-ranked candidates for tetramer ordering.

Formats the highest-risk candidates as a CSV a wet-lab teammate can hand
directly to a tetramer vendor / core facility, and that
``wetlab/validation_reports/build_concordance_report.py`` (roadmap
Sprint 4, not yet built) will eventually join real assay results back
against.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

_EXPORT_COLUMNS = [
    "sample_id",
    "tcr_cdr3_beta",
    "tcr_v_gene",
    "tcr_j_gene",
    "self_peptide",
    "hla_allele",
    "source_gene",
    "risk_score",
]


def main(risk_scores_path: str, export_out: str, top_n: int = 10) -> None:
    """Write the top-N highest-risk candidates as a tetramer-order CSV.

    Args:
        risk_scores_path: Path to ``risk_scores.json`` from the
            ``risk_scoring`` rule (already sorted descending by
            ``risk_score``).
        export_out: Where to write the export CSV.
        top_n: How many top-ranked candidates to include.

    Raises:
        ValueError: If ``top_n`` is not positive.
    """
    if top_n <= 0:
        raise ValueError("top_n must be positive.")

    candidates = json.loads(Path(risk_scores_path).read_text())
    top_candidates = candidates[:top_n]

    Path(export_out).parent.mkdir(parents=True, exist_ok=True)
    with Path(export_out).open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_EXPORT_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for candidate in top_candidates:
            writer.writerow(candidate)


if "snakemake" in globals():
    main(
        risk_scores_path=snakemake.input.risk_scores,  # noqa: F821
        export_out=snakemake.output.export,  # noqa: F821
        top_n=snakemake.params.top_n,  # noqa: F821
    )
