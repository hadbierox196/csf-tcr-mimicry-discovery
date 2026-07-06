"""Calibrates the mimicry risk score against tetramer-confirmed ground truth.

Consumes the join produced by
``wetlab/validation_reports/build_concordance_report.py`` once real
assay results exist. Until then, these functions operate on any
scores/labels arrays, so they can be exercised in tests using synthetic
data.
"""

from __future__ import annotations

import numpy as np
from sklearn.metrics import precision_score, recall_score, roc_auc_score


def evaluate_against_ground_truth(
    risk_scores: list[float],
    tetramer_positive: list[bool],
    threshold: float,
) -> dict[str, float]:
    """Compute concordance metrics between risk scores and tetramer results.

    Args:
        risk_scores: Predicted risk scores in ``[0, 1]``, one per
            tested candidate.
        tetramer_positive: Ground-truth tetramer-binding call for the
            same candidates, in matching order.
        threshold: Score threshold above which a prediction counts as
            "positive" for precision/recall (typically
            ``ScoringConfig.high_risk_threshold``).

    Returns:
        A dict with ``auroc``, ``precision``, and ``recall``. ``auroc``
        is ``float("nan")`` if all labels are the same class
        (undefined in that case).

    Raises:
        ValueError: If the two input lists have different lengths, or
            either is empty.
    """
    if len(risk_scores) != len(tetramer_positive):
        raise ValueError("risk_scores and tetramer_positive must be the same length.")
    if not risk_scores:
        raise ValueError("Cannot evaluate against an empty set of candidates.")

    scores = np.asarray(risk_scores)
    labels = np.asarray(tetramer_positive, dtype=bool)
    predictions = scores >= threshold

    n_classes = len(set(labels.tolist()))
    auroc = float(roc_auc_score(labels, scores)) if n_classes > 1 else float("nan")

    return {
        "auroc": auroc,
        "precision": float(precision_score(labels, predictions, zero_division=0)),
        "recall": float(recall_score(labels, predictions, zero_division=0)),
    }
