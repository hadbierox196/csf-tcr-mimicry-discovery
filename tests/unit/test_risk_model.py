"""Unit tests for mimicry_discovery.scoring.risk_model."""

from __future__ import annotations

import pytest

from mimicry_discovery.config import ScoringConfig
from mimicry_discovery.scoring.risk_model import (
    MimicryRiskInput,
    compute_risk_score,
    rank_candidates,
)


def test_compute_risk_score_strong_features_score_high_risk() -> None:
    """Strong similarity/confidence/conservation yields a high risk score."""
    features = MimicryRiskInput(
        peptide_similarity=0.95, anchor_conservation=1.0, structural_confidence=0.9
    )
    result = compute_risk_score(features, ScoringConfig())
    assert result.risk_score == pytest.approx(0.94)
    assert result.is_high_risk is True


def test_compute_risk_score_weak_features_score_low_risk() -> None:
    """Weak features across the board yield a low score, not flagged high risk."""
    features = MimicryRiskInput(
        peptide_similarity=0.1, anchor_conservation=0.0, structural_confidence=0.2
    )
    result = compute_risk_score(features, ScoringConfig())
    assert result.risk_score < 0.5
    assert result.is_high_risk is False


def test_mimicry_risk_input_rejects_out_of_range_feature() -> None:
    """A feature value outside [0, 1] is rejected at construction."""
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        MimicryRiskInput(peptide_similarity=1.5, anchor_conservation=0.5, structural_confidence=0.5)


def test_rank_candidates_sorts_by_descending_risk_score() -> None:
    """rank_candidates orders results from highest to lowest risk score."""
    config = ScoringConfig()
    low = (
        "cand-low",
        compute_risk_score(
            MimicryRiskInput(
                peptide_similarity=0.1, anchor_conservation=0.1, structural_confidence=0.1
            ),
            config,
        ),
    )
    high = (
        "cand-high",
        compute_risk_score(
            MimicryRiskInput(
                peptide_similarity=0.9, anchor_conservation=0.9, structural_confidence=0.9
            ),
            config,
        ),
    )

    ranked = rank_candidates([low, high])

    assert [candidate_id for candidate_id, _ in ranked] == ["cand-high", "cand-low"]
