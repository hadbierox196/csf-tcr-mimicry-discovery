"""Unit tests for mimicry_discovery.scoring.calibration."""

from __future__ import annotations

import math

import pytest

from mimicry_discovery.scoring.calibration import evaluate_against_ground_truth


def test_evaluate_against_ground_truth_perfect_separation() -> None:
    """Perfectly separated scores/labels yield AUROC, precision, recall of 1.0."""
    metrics = evaluate_against_ground_truth(
        risk_scores=[0.1, 0.2, 0.8, 0.9],
        tetramer_positive=[False, False, True, True],
        threshold=0.5,
    )
    assert metrics["auroc"] == pytest.approx(1.0)
    assert metrics["precision"] == pytest.approx(1.0)
    assert metrics["recall"] == pytest.approx(1.0)


def test_evaluate_against_ground_truth_single_class_auroc_is_nan() -> None:
    """AUROC is undefined (NaN) when all ground-truth labels are the same class."""
    metrics = evaluate_against_ground_truth(
        risk_scores=[0.1, 0.2], tetramer_positive=[False, False], threshold=0.5
    )
    assert math.isnan(metrics["auroc"])


def test_evaluate_against_ground_truth_mismatched_lengths_raises() -> None:
    """Mismatched-length inputs raise ValueError rather than a cryptic numpy error."""
    with pytest.raises(ValueError, match="same length"):
        evaluate_against_ground_truth([0.1, 0.2], [True], threshold=0.5)


def test_evaluate_against_ground_truth_empty_input_raises() -> None:
    """An empty candidate set raises ValueError rather than a downstream metric error."""
    with pytest.raises(ValueError, match="empty"):
        evaluate_against_ground_truth([], [], threshold=0.5)
