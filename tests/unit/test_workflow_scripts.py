"""Unit tests for the workflow/scripts/*.py Snakemake script bodies.

Each script's ``main()`` is plain Python taking explicit arguments -- see
each script's module docstring for why -- so these tests import and call
``main()`` directly. No Snakemake installation is required to run them.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import logging
import sys
from pathlib import Path
from types import ModuleType

import pytest

from mimicry_discovery.config import QcThresholds, ScoringConfig, StructurePredictionConfig
from mimicry_discovery.self_antigen.build_reference import SelfPeptide

_SCRIPTS_DIR = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
_FIXTURE_DIR = Path(__file__).resolve().parents[2] / "data" / "test"


def _load_script(name: str) -> ModuleType:
    """Import a workflow script as a plain module, by file path.

    Args:
        name: Script filename, e.g. ``"ingest_qc.py"``.

    Returns:
        The imported module (its top-level ``if "snakemake" in globals()``
        guard is never triggered here, since these tests don't inject a
        ``snakemake`` global).
    """
    spec = importlib.util.spec_from_file_location(name, _SCRIPTS_DIR / name)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def ingest_qc_script() -> ModuleType:
    """Load workflow/scripts/ingest_qc.py."""
    return _load_script("ingest_qc.py")


@pytest.fixture(scope="module")
def candidate_pairing_script() -> ModuleType:
    """Load workflow/scripts/candidate_pairing.py."""
    return _load_script("candidate_pairing.py")


@pytest.fixture(scope="module")
def structure_prediction_script() -> ModuleType:
    """Load workflow/scripts/structure_prediction.py."""
    return _load_script("structure_prediction.py")


@pytest.fixture(scope="module")
def risk_scoring_script() -> ModuleType:
    """Load workflow/scripts/risk_scoring.py."""
    return _load_script("risk_scoring.py")


@pytest.fixture(scope="module")
def wetlab_export_script() -> ModuleType:
    """Load workflow/scripts/wetlab_export.py."""
    return _load_script("wetlab_export.py")


def test_ingest_qc_main_parses_bundled_fixtures(
    tmp_path: Path, ingest_qc_script: ModuleType
) -> None:
    """ingest_qc.main() parses and QCs the bundled synthetic fixtures."""
    clonotypes_out = tmp_path / "clonotypes.json"
    neoantigens_out = tmp_path / "neoantigens.json"

    ingest_qc_script.main(
        tcr_vdj_path=str(_FIXTURE_DIR / "synthetic_tcr_clonotypes.csv"),
        neoantigen_report_path=str(_FIXTURE_DIR / "synthetic_neoantigens.tsv"),
        sample_id="synthetic-pt-001",
        qc=QcThresholds(),
        clonotypes_out=str(clonotypes_out),
        neoantigens_out=str(neoantigens_out),
    )

    clonotypes = json.loads(clonotypes_out.read_text())
    neoantigens = json.loads(neoantigens_out.read_text())
    assert len(clonotypes) == 2
    assert len(neoantigens) == 2
    assert clonotypes[0]["sample_id"] == "synthetic-pt-001"


def test_candidate_pairing_with_no_reference_yields_zero_candidates(
    tmp_path: Path, candidate_pairing_script: ModuleType, caplog: pytest.LogCaptureFixture
) -> None:
    """No self-peptide reference yet -> zero candidates, not a crash."""
    clonotypes_path = tmp_path / "c.json"
    neoantigens_path = tmp_path / "n.json"
    clonotypes_path.write_text(json.dumps([{
        "sample_id": "s1", "cdr3_aa": "AAA", "v_gene": "V1", "j_gene": "J1",
        "chain": "TRB", "umi_count": 5, "clonal_frequency": 1.0,
    }]))
    neoantigens_path.write_text(json.dumps([{
        "sample_id": "s1", "peptide": "SLLMWITQC", "hla_allele": "HLA-A*02:01",
        "wildtype_peptide": None, "gene_name": "TP53", "binding_affinity_nm": 40.0,
    }]))
    candidates_out = tmp_path / "candidates.json"

    with caplog.at_level(logging.WARNING):
        candidate_pairing_script.main(
            clonotypes_path=str(clonotypes_path),
            neoantigens_path=str(neoantigens_path),
            reference_dir=str(tmp_path / "nonexistent_reference"),
            candidates_out=str(candidates_out),
        )

    assert json.loads(candidates_out.read_text()) == []
    assert "still a stub" in caplog.text


def test_candidate_pairing_with_reference_generates_full_cross_product(
    tmp_path: Path, candidate_pairing_script: ModuleType
) -> None:
    """With a reference present, every matching (clonotype, neoantigen) pair
    that shares HLA allele and peptide length becomes a candidate."""
    clonotypes_path = tmp_path / "c.json"
    neoantigens_path = tmp_path / "n.json"
    clonotypes_path.write_text(json.dumps([
        {"sample_id": "s1", "cdr3_aa": "AAA", "v_gene": "V1", "j_gene": "J1",
         "chain": "TRB", "umi_count": 5, "clonal_frequency": 0.5},
        {"sample_id": "s1", "cdr3_aa": "BBB", "v_gene": "V2", "j_gene": "J2",
         "chain": "TRB", "umi_count": 5, "clonal_frequency": 0.5},
    ]))
    neoantigens_path.write_text(json.dumps([
        {"sample_id": "s1", "peptide": "SLLMWITQC", "hla_allele": "HLA-A*02:01",
         "wildtype_peptide": None, "gene_name": "TP53", "binding_affinity_nm": 40.0},
        {"sample_id": "s1", "peptide": "GLATEKSRW", "hla_allele": "HLA-A*02:01",
         "wildtype_peptide": None, "gene_name": "BRAF", "binding_affinity_nm": 45.0},
    ]))
    reference_dir = tmp_path / "reference"
    reference_dir.mkdir()
    (reference_dir / "neuronal_self_peptidome.json").write_text(json.dumps([
        {"peptide": "SLLMWITQC", "hla_allele": "HLA-A*02:01",
         "source_gene": "MBP", "evidence": "predicted"},
    ]))
    candidates_out = tmp_path / "candidates.json"

    candidate_pairing_script.main(
        clonotypes_path=str(clonotypes_path),
        neoantigens_path=str(neoantigens_path),
        reference_dir=str(reference_dir),
        candidates_out=str(candidates_out),
    )

    candidates = json.loads(candidates_out.read_text())
    # 2 clonotypes x 2 neoantigens, both neoantigens are 9-mers on
    # HLA-A*02:01 so both match the single reference entry -> 4 candidates.
    assert len(candidates) == 4
    assert all(c["self_peptide"] == "SLLMWITQC" for c in candidates)


def test_structure_prediction_smoke_test_mode(
    tmp_path: Path, structure_prediction_script: ModuleType
) -> None:
    """smoke_test=True substitutes a fixed placeholder, never the real backend."""
    candidates_path = tmp_path / "candidates.json"
    candidates_path.write_text(json.dumps([{
        "sample_id": "s1", "tcr_cdr3_beta": "AAA", "tcr_v_gene": "V1", "tcr_j_gene": "J1",
        "tumor_peptide": "SLLMWITQV", "self_peptide": "SLLMWITQC",
        "hla_allele": "HLA-A*02:01", "source_gene": "MBP",
    }]))
    results_out = tmp_path / "results.json"

    structure_prediction_script.main(
        candidates_path=str(candidates_path),
        structure_config=StructurePredictionConfig(),
        output_dir=str(tmp_path / "structures"),
        results_out=str(results_out),
        smoke_test=True,
    )

    results = json.loads(results_out.read_text())
    assert results[0]["structural_confidence"] == 0.5
    assert results[0]["backend"] == "smoke_test_placeholder"


def test_structure_prediction_production_mode_raises_not_implemented(
    tmp_path: Path, structure_prediction_script: ModuleType
) -> None:
    """Without smoke_test, the real (still-scaffolded) backend is called and
    correctly raises NotImplementedError -- this is expected, not a bug,
    until roadmap Sprint 2/3 lands."""
    candidates_path = tmp_path / "candidates.json"
    candidates_path.write_text(json.dumps([{
        "sample_id": "s1", "tcr_cdr3_beta": "AAA", "tcr_v_gene": "V1", "tcr_j_gene": "J1",
        "tumor_peptide": "SLLMWITQV", "self_peptide": "SLLMWITQC",
        "hla_allele": "HLA-A*02:01", "source_gene": "MBP",
    }]))

    with pytest.raises(NotImplementedError):
        structure_prediction_script.main(
            candidates_path=str(candidates_path),
            structure_config=StructurePredictionConfig(),
            output_dir=str(tmp_path / "structures"),
            results_out=str(tmp_path / "results.json"),
            smoke_test=False,
        )


def test_risk_scoring_ranks_candidates_descending(
    tmp_path: Path, risk_scoring_script: ModuleType
) -> None:
    """risk_scoring.main() scores every candidate and sorts descending."""
    structure_results_path = tmp_path / "structure_results.json"
    structure_results_path.write_text(json.dumps([
        {"sample_id": "s1", "tcr_cdr3_beta": "A", "tumor_peptide": "SLLMWITQV",
         "self_peptide": "SLLMWITQC", "hla_allele": "HLA-A*02:01", "source_gene": "MBP",
         "structural_confidence": 0.8, "backend": "smoke_test_placeholder"},
        {"sample_id": "s1", "tcr_cdr3_beta": "B", "tumor_peptide": "GPDEAKPRW",
         "self_peptide": "DDDDDDDDD", "hla_allele": "HLA-A*02:01", "source_gene": "MOG",
         "structural_confidence": 0.1, "backend": "smoke_test_placeholder"},
    ]))
    risk_scores_out = tmp_path / "risk_scores.json"

    risk_scoring_script.main(
        structure_results_path=str(structure_results_path),
        scoring_config=ScoringConfig(),
        risk_scores_out=str(risk_scores_out),
    )

    scored = json.loads(risk_scores_out.read_text())
    assert len(scored) == 2
    assert scored[0]["risk_score"] >= scored[1]["risk_score"]
    assert scored[0]["tcr_cdr3_beta"] == "A"


def test_wetlab_export_writes_top_n_as_csv(
    tmp_path: Path, wetlab_export_script: ModuleType
) -> None:
    """wetlab_export.main() writes exactly top_n rows as a tetramer-order CSV."""
    risk_scores_path = tmp_path / "risk_scores.json"
    risk_scores_path.write_text(json.dumps([
        {"sample_id": "s1", "tcr_cdr3_beta": "A", "tcr_v_gene": "V1", "tcr_j_gene": "J1",
         "self_peptide": "PEP1", "hla_allele": "HLA-A*02:01", "source_gene": "MBP",
         "risk_score": 0.9, "is_high_risk": True},
        {"sample_id": "s1", "tcr_cdr3_beta": "B", "tcr_v_gene": "V2", "tcr_j_gene": "J2",
         "self_peptide": "PEP2", "hla_allele": "HLA-A*02:01", "source_gene": "MOG",
         "risk_score": 0.3, "is_high_risk": False},
    ]))
    export_out = tmp_path / "export.csv"

    wetlab_export_script.main(
        risk_scores_path=str(risk_scores_path), export_out=str(export_out), top_n=1
    )

    with export_out.open() as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 1
    assert rows[0]["tcr_cdr3_beta"] == "A"
    assert "is_high_risk" not in rows[0]


def test_wetlab_export_rejects_non_positive_top_n(
    tmp_path: Path, wetlab_export_script: ModuleType
) -> None:
    """A non-positive top_n raises ValueError rather than silently writing nothing."""
    risk_scores_path = tmp_path / "risk_scores.json"
    risk_scores_path.write_text("[]")

    with pytest.raises(ValueError, match="positive"):
        wetlab_export_script.main(
            risk_scores_path=str(risk_scores_path),
            export_out=str(tmp_path / "export.csv"),
            top_n=0,
        )
