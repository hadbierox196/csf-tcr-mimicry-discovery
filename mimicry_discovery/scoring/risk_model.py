"""Composite molecular-mimicry risk score.

Combines peptide similarity, structural interface confidence, and
anchor conservation into a single risk score via a transparent,
config-weighted linear combination. This rule-based model is the
default. TODO(scoring): once tetramer-validated labels exist (see
``wetlab/validation_reports/``), swap in a learned model here (e.g.
gradient-boosted trees) calibrated against those labels -- keep the
same ``MimicryRiskInput``/``MimicryRiskResult`` contract so downstream
Snakemake rules and reports don't need to change.
"""

from __future__ import annotations

from dataclasses import dataclass

from mimicry_discovery.config import ScoringConfig


@dataclass(frozen=True)
class MimicryRiskInput:
    """Feature values for one TCR x tumor-peptide x self-peptide triple.

    Attributes:
        peptide_similarity: Output of
            :func:`mimicry_discovery.scoring.features.peptide_similarity`,
            comparing the tumor peptide to the candidate self-peptide.
        anchor_conservation: Output of
            :func:`mimicry_discovery.scoring.features.anchor_conservation`
            for the same peptide pair.
        structural_confidence: Normalized (0-1) interface confidence
            for the TCR:self-peptide-HLA complex, typically the
            :func:`mimicry_discovery.structure.interface_metrics.mean_plddt`
            interface value divided by 100.
    """

    peptide_similarity: float
    anchor_conservation: float
    structural_confidence: float

    def __post_init__(self) -> None:
        """Validate all three features are normalized to [0, 1].

        Raises:
            ValueError: If any feature value is outside ``[0, 1]``.
        """
        for name in ("peptide_similarity", "anchor_conservation", "structural_confidence"):
            value = getattr(self, name)
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be in [0, 1], got {value}.")


@dataclass(frozen=True)
class MimicryRiskResult:
    """Output of the mimicry risk-scoring model.

    Attributes:
        risk_score: Composite score in ``[0, 1]``; higher indicates
            stronger predicted cross-reactivity.
        is_high_risk: Whether ``risk_score`` clears the configured
            ``high_risk_threshold``.
    """

    risk_score: float
    is_high_risk: bool


def compute_risk_score(
    features: MimicryRiskInput, scoring_config: ScoringConfig
) -> MimicryRiskResult:
    """Compute the composite molecular-mimicry risk score.

    Args:
        features: Feature values for one candidate triple.
        scoring_config: Weights and the high-risk threshold, typically
            loaded from ``config/config.yaml`` via
            :func:`mimicry_discovery.config.load_config`.

    Returns:
        A :class:`MimicryRiskResult` with the composite score and a
        boolean high-risk flag.
    """
    weighted_sum = (
        features.peptide_similarity * scoring_config.peptide_similarity_weight
        + features.structural_confidence * scoring_config.structural_confidence_weight
        + features.anchor_conservation * scoring_config.anchor_conservation_weight
    )
    total_weight = (
        scoring_config.peptide_similarity_weight
        + scoring_config.structural_confidence_weight
        + scoring_config.anchor_conservation_weight
    )
    risk_score = weighted_sum / total_weight if total_weight > 0 else 0.0
    risk_score = min(1.0, max(0.0, risk_score))
    return MimicryRiskResult(
        risk_score=risk_score,
        is_high_risk=risk_score >= scoring_config.high_risk_threshold,
    )


def rank_candidates(
    results: list[tuple[str, MimicryRiskResult]],
) -> list[tuple[str, MimicryRiskResult]]:
    """Sort candidate results by descending risk score.

    Args:
        results: ``(candidate_id, result)`` pairs, e.g. keyed by a
            ``f"{sample_id}:{cdr3}:{peptide}"``-style identifier.

    Returns:
        The same pairs, sorted by descending ``risk_score``.
    """
    return sorted(results, key=lambda pair: pair[1].risk_score, reverse=True)
