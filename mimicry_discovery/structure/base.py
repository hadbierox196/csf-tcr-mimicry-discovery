"""Abstract interface all TCR-pMHC structure predictors implement.

Concrete backends (TCRdock, ESMFold, AlphaFold-Multimer) are selected at
runtime via ``config.structure_prediction.backend``. Every backend
consumes the same request type and returns the same result type, so the
scoring stage never needs to know which backend produced a given
structure.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass


@dataclass(frozen=True)
class StructurePredictionRequest:
    """Inputs needed to model one candidate TCR-pMHC complex.

    Beta-chain fields are required; alpha-chain fields are optional,
    since not every clonotype recovers both chains (see
    :class:`mimicry_discovery.io.tcr_parsers.TCRClonotype`). Backends
    that need the full paired receptor (TCRdock) should reject a request
    with no alpha chain rather than silently modeling beta alone -- see
    :class:`~mimicry_discovery.structure.tcrdock_adapter.TCRDockPredictor`.

    Attributes:
        tcr_cdr3_beta: Beta-chain CDR3 amino acid sequence.
        tcr_v_gene: Beta-chain V-gene segment (IMGT nomenclature).
        tcr_j_gene: Beta-chain J-gene segment (IMGT nomenclature).
        peptide: Candidate presented peptide sequence.
        hla_allele: Presenting HLA allele, two-field IPD-IMGT/HLA form.
        tcr_cdr3_alpha: Alpha-chain CDR3 amino acid sequence, if the
            clonotype's alpha chain was recovered.
        tcr_v_gene_alpha: Alpha-chain V-gene segment, if recovered.
        tcr_j_gene_alpha: Alpha-chain J-gene segment, if recovered.
        organism: ``"human"`` or ``"mouse"`` -- TCRdock needs this to
            select the right germline gene database.
        mhc_class: ``1`` or ``2``.
    """

    tcr_cdr3_beta: str
    tcr_v_gene: str
    tcr_j_gene: str
    peptide: str
    hla_allele: str
    tcr_cdr3_alpha: str | None = None
    tcr_v_gene_alpha: str | None = None
    tcr_j_gene_alpha: str | None = None
    organism: str = "human"
    mhc_class: int = 1

    @property
    def is_paired(self) -> bool:
        """Whether both alpha and beta chain info are present."""
        return self.tcr_cdr3_alpha is not None


@dataclass(frozen=True)
class StructurePredictionResult:
    """Output of a structure-prediction backend for one complex.

    Attributes:
        pdb_path: Path to the predicted complex structure (PDB).
        backend: Name of the backend that produced this result.
        mean_plddt: Mean per-residue confidence (0-100) across the
            modeled complex, when the backend provides one (ESMFold,
            AlphaFold-Multimer); ``None`` for backends that don't.
        interface_plddt: Mean confidence restricted to interface
            residues, if computed at prediction time; otherwise this is
            left to
            :mod:`mimicry_discovery.structure.interface_metrics` to
            compute post hoc.
    """

    pdb_path: str
    backend: str
    mean_plddt: float | None = None
    interface_plddt: float | None = None

    def __post_init__(self) -> None:
        """Validate that any provided confidence values are in [0, 100].

        Raises:
            ValueError: If ``mean_plddt`` or ``interface_plddt`` is
                provided but outside ``[0, 100]``.
        """
        for name in ("mean_plddt", "interface_plddt"):
            value = getattr(self, name)
            if value is not None and not 0.0 <= value <= 100.0:
                raise ValueError(f"{name} must be in [0, 100] if provided, got {value}.")


class StructurePredictor(ABC):
    """Common interface for TCR-pMHC complex structure predictors."""

    @abstractmethod
    def predict(self, request: StructurePredictionRequest) -> StructurePredictionResult:
        """Predict the 3D structure of a TCR-pMHC complex.

        Args:
            request: The TCR/peptide/HLA combination to model.

        Returns:
            A :class:`StructurePredictionResult` describing the modeled
            complex and its confidence metrics.

        Raises:
            RuntimeError: If the underlying model fails to produce a
                structure (implementation-specific).
        """
        raise NotImplementedError
