"""Unit tests for mimicry_discovery.structure.base and its adapters.

The three real backends (TCRdock, ESMFold, AlphaFold-Multimer) are
intentionally-scaffolded stubs -- see their module docstrings. These
tests assert the *documented* stub behavior (a clear NotImplementedError)
plus the shared ABC contract, using a minimal fake predictor to exercise
the contract itself.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mimicry_discovery.structure.af_multimer_adapter import AlphaFoldMultimerPredictor
from mimicry_discovery.structure.base import (
    StructurePredictionRequest,
    StructurePredictionResult,
    StructurePredictor,
)
from mimicry_discovery.structure.esmfold_adapter import ESMFoldPredictor
from mimicry_discovery.structure.tcrdock_adapter import TCRDockPredictor

_SAMPLE_REQUEST = StructurePredictionRequest(
    tcr_cdr3_beta="CASSLGQGNTIYF",
    tcr_v_gene="TRBV20-1",
    tcr_j_gene="TRBJ2-1",
    peptide="SLLMWITQV",
    hla_allele="HLA-A*02:01",
)


class _FakePredictor(StructurePredictor):
    """Minimal concrete StructurePredictor used only to test the ABC contract."""

    def predict(self, request: StructurePredictionRequest) -> StructurePredictionResult:
        """Return a fixed, fake result without doing any real modeling."""
        return StructurePredictionResult(pdb_path="fake.pdb", backend="fake", mean_plddt=87.5)


def test_structure_predictor_cannot_be_instantiated_directly() -> None:
    """StructurePredictor is an ABC and cannot be instantiated on its own."""
    with pytest.raises(TypeError):
        StructurePredictor()  # type: ignore[abstract]


def test_fake_predictor_satisfies_the_contract() -> None:
    """A concrete subclass implementing predict() works as expected."""
    result = _FakePredictor().predict(_SAMPLE_REQUEST)
    assert result.backend == "fake"
    assert result.mean_plddt == pytest.approx(87.5)


def test_structure_prediction_result_validates_plddt_range() -> None:
    """A pLDDT value outside [0, 100] is rejected at construction."""
    with pytest.raises(ValueError, match=r"\[0, 100\]"):
        StructurePredictionResult(pdb_path="x.pdb", backend="fake", mean_plddt=150.0)


def test_scaffolded_backend_raises_not_implemented(tmp_path: Path) -> None:
    """AlphaFold-Multimer is a clearly-marked scaffold until wired in.
    TCRDockPredictor and ESMFoldPredictor are covered separately (in
    test_tcrdock_adapter.py and test_esmfold_adapter.py) -- both now have
    real, input-dependent behavior (pairing / HLA-lookup requirements)
    rather than a single unconditional NotImplementedError."""
    predictor = AlphaFoldMultimerPredictor(output_dir=tmp_path)
    with pytest.raises(NotImplementedError):
        predictor.predict(_SAMPLE_REQUEST)


@pytest.mark.parametrize(
    "predictor_cls",
    [TCRDockPredictor, ESMFoldPredictor, AlphaFoldMultimerPredictor],
)
def test_adapter_constructors_create_output_dir(
    tmp_path: Path, predictor_cls: type[StructurePredictor]
) -> None:
    """Each adapter creates its output directory on construction."""
    output_dir = tmp_path / "nested" / "predictions"
    predictor_cls(output_dir=output_dir)  # type: ignore[call-arg]
    assert output_dir.exists()
