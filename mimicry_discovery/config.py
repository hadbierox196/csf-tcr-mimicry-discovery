"""Typed, validated configuration for the mimicry-discovery pipeline.

Configuration is authored as YAML under ``config/config.yaml`` and loaded
through :func:`load_config`. We use stdlib ``dataclasses`` with
``__post_init__`` validation rather than Pydantic: it keeps the runtime
dependency footprint smaller, matches the plain-dataclass style used
throughout AlphaFold's own codebase, and Snakemake already owns config
loading/override semantics via its own profile system (see
``docs/adr/0002-dataclasses-over-pydantic.md``), so a second, heavier
validation framework would duplicate responsibility without adding much.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields, is_dataclass
from enum import Enum
from pathlib import Path
from typing import Any, TypeVar, get_type_hints

import yaml

T = TypeVar("T")


class StructureBackend(str, Enum):
    """Supported TCR-pMHC structure prediction backends."""

    TCRDOCK = "tcrdock"
    ESMFOLD = "esmfold"
    AF_MULTIMER = "af_multimer"


@dataclass(frozen=True)
class PathsConfig:
    """Filesystem locations the pipeline reads from and writes to.

    All paths are read from config rather than hardcoded so the same
    rule code runs unmodified under the SLURM and cloud Snakemake
    profiles -- only ``config/profiles/{slurm,cloud}`` changes.

    Attributes:
        samples_tsv: Path to the sample manifest.
        reference_dir: Directory holding DVC-tracked reference data.
        models_dir: Directory holding DVC-tracked model weights.
        output_dir: Directory pipeline outputs are written to.
    """

    samples_tsv: Path
    reference_dir: Path
    models_dir: Path
    output_dir: Path

    def __post_init__(self) -> None:
        """Coerce string paths to Path objects."""
        for f in fields(self):
            value = getattr(self, f.name)
            if not isinstance(value, Path):
                object.__setattr__(self, f.name, Path(value))


@dataclass(frozen=True)
class QcThresholds:
    """Quality-control cutoffs applied during TCR-seq ingestion.

    Attributes:
        min_umi_count: Minimum summed UMI support to keep a clonotype.
        min_clonal_frequency: Minimum within-sample clonal frequency
            (0-1) to keep a clonotype.
        require_productive: Drop non-productive rearrangements.
        require_high_confidence: Drop low-confidence contig calls.
    """

    min_umi_count: int = 2
    min_clonal_frequency: float = 0.0
    require_productive: bool = True
    require_high_confidence: bool = True

    def __post_init__(self) -> None:
        """Validate threshold ranges.

        Raises:
            ValueError: If ``min_umi_count`` is negative, or
                ``min_clonal_frequency`` is outside ``[0, 1]``.
        """
        if self.min_umi_count < 0:
            raise ValueError("min_umi_count must be >= 0.")
        if not 0.0 <= self.min_clonal_frequency <= 1.0:
            raise ValueError("min_clonal_frequency must be in [0, 1].")


@dataclass(frozen=True)
class StructurePredictionConfig:
    """Settings controlling the structure-prediction stage.

    Attributes:
        backend: Which structure-prediction backend to use.
        num_recycles: Number of structure-refinement recycles.
        device: Torch device string, e.g. ``"cuda"`` or ``"cpu"``.
    """

    backend: StructureBackend = StructureBackend.TCRDOCK
    num_recycles: int = 3
    device: str = "cuda"

    def __post_init__(self) -> None:
        """Coerce a raw string backend value to StructureBackend.

        Raises:
            ValueError: If ``num_recycles`` is negative, or ``backend``
                is not a recognized value.
        """
        if not isinstance(self.backend, StructureBackend):
            object.__setattr__(self, "backend", StructureBackend(self.backend))
        if self.num_recycles < 0:
            raise ValueError("num_recycles must be >= 0.")


@dataclass(frozen=True)
class ScoringConfig:
    """Weights and thresholds for the composite mimicry risk score.

    Attributes:
        peptide_similarity_weight: Weight on the peptide-similarity
            feature.
        structural_confidence_weight: Weight on the structural
            interface-confidence feature.
        anchor_conservation_weight: Weight on the HLA anchor-residue
            conservation feature.
        high_risk_threshold: Score at/above which a candidate is
            flagged high-risk.
    """

    peptide_similarity_weight: float = 0.4
    structural_confidence_weight: float = 0.4
    anchor_conservation_weight: float = 0.2
    high_risk_threshold: float = 0.75

    def __post_init__(self) -> None:
        """Validate weights and threshold are in sensible ranges.

        Raises:
            ValueError: If any weight is negative, or
                ``high_risk_threshold`` is exactly 0.0 or 1.0 (which
                would flag either everything or nothing as high risk).
        """
        for name in (
            "peptide_similarity_weight",
            "structural_confidence_weight",
            "anchor_conservation_weight",
        ):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} must be >= 0.")
        if self.high_risk_threshold in (0.0, 1.0):
            raise ValueError(
                "high_risk_threshold of exactly 0.0 or 1.0 flags all-or-"
                "nothing results; choose a value strictly between them."
            )


@dataclass(frozen=True)
class PipelineConfig:
    """Root configuration object for the mimicry-discovery pipeline.

    Attributes:
        paths: Filesystem locations (see :class:`PathsConfig`).
        qc: TCR-seq QC thresholds.
        structure_prediction: Structure-prediction backend settings.
        scoring: Mimicry risk-score weights and threshold.
    """

    paths: PathsConfig
    qc: QcThresholds = field(default_factory=QcThresholds)
    structure_prediction: StructurePredictionConfig = field(
        default_factory=StructurePredictionConfig
    )
    scoring: ScoringConfig = field(default_factory=ScoringConfig)


def _dataclass_from_dict(cls: type[T], data: dict[str, Any]) -> T:
    """Recursively build a (possibly nested) dataclass from a plain dict.

    Args:
        cls: The dataclass type to construct.
        data: A dict of field values, as loaded from YAML.

    Returns:
        An instance of ``cls``.
    """
    resolved_types = get_type_hints(cls)
    kwargs: dict[str, Any] = {}
    for f in fields(cls):  # type: ignore[arg-type]
        if f.name not in data:
            continue
        value = data[f.name]
        field_type = resolved_types[f.name]
        if is_dataclass(field_type) and isinstance(value, dict):
            kwargs[f.name] = _dataclass_from_dict(field_type, value)  # type: ignore[arg-type]
        else:
            kwargs[f.name] = value
    return cls(**kwargs)  # type: ignore[return-value]


def load_config(config_path: Path | str) -> PipelineConfig:
    """Load and validate a pipeline configuration from a YAML file.

    Args:
        config_path: Path to a ``config.yaml`` file. Typically the same
            file passed to Snakemake via ``--configfile``.

    Returns:
        A validated :class:`PipelineConfig` instance.

    Raises:
        FileNotFoundError: If ``config_path`` does not exist.
        ValueError: If the YAML content fails validation (e.g. an
            out-of-range weight or a missing required field).
    """
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    with config_path.open("r", encoding="utf-8") as fh:
        raw = yaml.safe_load(fh) or {}
    try:
        return _dataclass_from_dict(PipelineConfig, raw)
    except TypeError as exc:
        raise ValueError(f"Invalid config in {config_path}: {exc}") from exc
