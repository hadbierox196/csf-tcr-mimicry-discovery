"""Unit tests for mimicry_discovery.structure.esmfold_adapter."""

from __future__ import annotations

from pathlib import Path

import pytest

from mimicry_discovery.structure.base import StructurePredictionRequest
from mimicry_discovery.structure.esmfold_adapter import ESMFoldPredictor

_REQUEST = StructurePredictionRequest(
    tcr_cdr3_beta="CASSIRSSYEQYF", tcr_v_gene="TRBV19", tcr_j_gene="TRBJ2-7",
    peptide="SLLMWITQV", hla_allele="HLA-A*02:01",
)


def test_predict_without_hla_lookup_raises_key_error(tmp_path: Path) -> None:
    """No HLA reference loaded -> a clear KeyError, checked before model
    loading is even attempted (fail on missing input, not a heavy model
    that isn't wired in yet either)."""
    predictor = ESMFoldPredictor(output_dir=tmp_path)
    with pytest.raises(KeyError, match="No reference sequence"):
        predictor.predict(_REQUEST)


def test_predict_with_hla_lookup_reaches_model_loading_not_implemented(tmp_path: Path) -> None:
    """Once the HLA lookup covers the requested allele, the request passes
    that check and reaches the (still-scaffolded) model-loading step."""
    predictor = ESMFoldPredictor(
        output_dir=tmp_path, hla_sequence_lookup={"A*02:01": "MAVMAPRTLLL"}
    )
    with pytest.raises(NotImplementedError):
        predictor.predict(_REQUEST)


def test_build_complex_sequence_joins_chains_with_colons(tmp_path: Path) -> None:
    """The complex sequence is beta-CDR3:peptide:HLA-sequence, colon-joined."""
    predictor = ESMFoldPredictor(
        output_dir=tmp_path, hla_sequence_lookup={"A*02:01": "MAVMAPRTLLL"}
    )
    sequence = predictor._build_complex_sequence(_REQUEST)
    assert sequence == "CASSIRSSYEQYF:SLLMWITQV:MAVMAPRTLLL"


def test_from_reference_fasta_loads_and_collapses_to_two_field(tmp_path: Path) -> None:
    """from_reference_fasta wires a real FASTA straight into the lookup."""
    fasta_path = tmp_path / "A_prot.fasta"
    fasta_path.write_text(">HLA:HLA00001 A*02:01:01:01 11 aa\nMAVMAPRTLLL\n")

    predictor = ESMFoldPredictor.from_reference_fasta(
        output_dir=tmp_path, hla_fasta_path=fasta_path
    )

    assert predictor.hla_sequence_lookup == {"A*02:01": "MAVMAPRTLLL"}
    with pytest.raises(NotImplementedError):
        predictor.predict(_REQUEST)  # lookup now covers this allele
