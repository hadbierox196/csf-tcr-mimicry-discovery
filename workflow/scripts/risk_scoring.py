"""Snakemake script: compute the composite mimicry risk score per candidate.

Unlike ingest_qc/candidate_pairing/structure_prediction, every function
this script calls (``mimicry_discovery.scoring``) is fully implemented,
not scaffolded -- this is the one rule in the DAG that needs no smoke-test
mode or NotImplementedError caveat.
"""

from __future__ import annotations

import json
from pathlib import Path

from mimicry_discovery.config import ScoringConfig
from mimicry_discovery.scoring.features import anchor_conservation, peptide_similarity
from mimicry_discovery.scoring.risk_model import MimicryRiskInput, compute_risk_score


def main(
    structure_results_path: str,
    scoring_config: ScoringConfig,
    risk_scores_out: str,
) -> None:
    """Score every candidate and write ranked results.

    Args:
        structure_results_path: Path to ``structure_results.json`` from
            the ``structure_prediction`` rule.
        scoring_config: Feature weights and the high-risk threshold.
        risk_scores_out: Where to write ranked risk scores (JSON).
    """
    candidates = json.loads(Path(structure_results_path).read_text())

    scored = []
    for candidate in candidates:
        similarity = peptide_similarity(candidate["tumor_peptide"], candidate["self_peptide"])
        anchors = anchor_conservation(candidate["tumor_peptide"], candidate["self_peptide"])
        features = MimicryRiskInput(
            peptide_similarity=similarity,
            anchor_conservation=anchors,
            structural_confidence=candidate["structural_confidence"],
        )
        result = compute_risk_score(features, scoring_config)
        scored.append(
            {
                **candidate,
                "peptide_similarity": similarity,
                "anchor_conservation": anchors,
                "risk_score": result.risk_score,
                "is_high_risk": result.is_high_risk,
            }
        )

    scored.sort(key=lambda c: c["risk_score"], reverse=True)

    Path(risk_scores_out).parent.mkdir(parents=True, exist_ok=True)
    Path(risk_scores_out).write_text(json.dumps(scored, indent=2))


if "snakemake" in globals():
    main(
        structure_results_path=snakemake.input.results,  # noqa: F821
        scoring_config=snakemake.params.scoring_config,  # noqa: F821
        risk_scores_out=snakemake.output.risk_scores,  # noqa: F821
    )
