"""End-to-end smoke test exercising the currently-implemented pipeline
stages against the bundled synthetic fixtures in ``data/test/``.

TODO(pipeline): once ``workflow/Snakefile`` exists, add a sibling test
that invokes ``snakemake --configfile ... -n`` (dry run) or a real run
against these same fixtures, rather than calling package functions
directly the way this test does. Structure prediction is stubbed (see
``mimicry_discovery/structure/*_adapter.py``), so this test uses a fixed
placeholder ``structural_confidence`` rather than a real predicted
value -- clearly marked below.

This is also what CI's ``smoke-test`` job runs (see
``.github/workflows/ci.yml``): if this test passes, a fresh clone can
install the package and run every currently-implemented stage
end-to-end without touching real patient data or a GPU.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mimicry_discovery.config import QcThresholds, ScoringConfig
from mimicry_discovery.io.neoantigen_parsers import parse_pvacseq_report
from mimicry_discovery.io.tcr_parsers import parse_10x_vdj, repertoire_summary
from mimicry_discovery.lineage.provenance import write_provenance_manifest
from mimicry_discovery.scoring.features import anchor_conservation, peptide_similarity
from mimicry_discovery.scoring.risk_model import (
    MimicryRiskInput,
    MimicryRiskResult,
    compute_risk_score,
    rank_candidates,
)

FIXTURE_DIR = Path(__file__).resolve().parents[2] / "data" / "test"
_PLACEHOLDER_STRUCTURAL_CONFIDENCE = 0.5  # TODO(pipeline): replace with a real
# mimicry_discovery.structure.interface_metrics.mean_plddt(...) / 100 value
# once a structure-prediction backend is wired in.


@pytest.fixture
def tcr_fixture_path() -> Path:
    """Path to the bundled synthetic 10x VDJ fixture."""
    path = FIXTURE_DIR / "synthetic_tcr_clonotypes.csv"
    assert path.exists(), f"Missing fixture: {path}"
    return path


@pytest.fixture
def neoantigen_fixture_path() -> Path:
    """Path to the bundled synthetic pVACseq fixture."""
    path = FIXTURE_DIR / "synthetic_neoantigens.tsv"
    assert path.exists(), f"Missing fixture: {path}"
    return path


def test_ingest_stage_parses_and_qcs_the_bundled_fixtures(
    tcr_fixture_path: Path, neoantigen_fixture_path: Path
) -> None:
    """io.tcr_parsers and io.neoantigen_parsers run end-to-end on the fixtures.

    The fixture deliberately includes one non-productive and one
    low-confidence contig to confirm QC actually filters them, and one
    neoantigen weak enough (812nM) to be dropped by the default
    affinity cutoff.
    """
    clonotypes = parse_10x_vdj(
        tcr_fixture_path, sample_id="synthetic-pt-001", qc=QcThresholds()
    )
    neoantigens = parse_pvacseq_report(
        neoantigen_fixture_path, sample_id="synthetic-pt-001"
    )

    assert len(clonotypes) == 3  # the non-productive/low-confidence contigs are dropped
    assert len(neoantigens) == 2  # the 812nM KRAS call is dropped by the affinity cutoff
    assert {c.cdr3_beta for c in clonotypes} == {
        "CASSLGQGNTIYF", "CASSPGQGAYEQYF", "CASSIRSSYEQYF",
    }
    assert {n.gene_name for n in neoantigens} == {"TP53", "BRAF"}

    summary = repertoire_summary(clonotypes)
    assert summary["n_clonotypes"] == 3
    assert summary["n_paired"] == 1  # one clonotype recovered both chains
    assert summary["total_umis"] == 36


def test_scoring_stage_runs_end_to_end_and_ranks_candidates(
    tcr_fixture_path: Path, neoantigen_fixture_path: Path
) -> None:
    """scoring.features + scoring.risk_model run on parsed neoantigen pairs.

    Uses each neoantigen's mutant peptide vs. its own wild-type peptide
    as the similarity comparison -- a real, meaningful pair (not
    arbitrary dummy data), since point-mutation neoepitopes are
    guaranteed equal length to their wild-type counterpart.
    """
    neoantigens = parse_pvacseq_report(
        neoantigen_fixture_path, sample_id="synthetic-pt-001"
    )
    scoring_config = ScoringConfig()

    results: list[tuple[str, MimicryRiskResult]] = []
    for neoantigen in neoantigens:
        assert neoantigen.wildtype_peptide is not None
        similarity = peptide_similarity(neoantigen.peptide, neoantigen.wildtype_peptide)
        anchors = anchor_conservation(neoantigen.peptide, neoantigen.wildtype_peptide)
        features = MimicryRiskInput(
            peptide_similarity=similarity,
            anchor_conservation=anchors,
            structural_confidence=_PLACEHOLDER_STRUCTURAL_CONFIDENCE,
        )
        results.append((neoantigen.peptide, compute_risk_score(features, scoring_config)))

    ranked = rank_candidates(results)

    assert [candidate_id for candidate_id, _ in ranked] == ["GLATEKSRW", "SLLMWITQC"]
    scores = {candidate_id: result.risk_score for candidate_id, result in ranked}
    assert scores["SLLMWITQC"] == pytest.approx(0.6556, abs=1e-3)
    assert scores["GLATEKSRW"] == pytest.approx(0.6778, abs=1e-3)


def test_provenance_manifest_written_for_the_smoke_test_run(
    tmp_path: Path, tcr_fixture_path: Path
) -> None:
    """A provenance manifest can be written referencing the fixture inputs.

    Exercises the lineage-tracking hook every real pipeline run is
    expected to call (see the architecture rationale on data lineage).
    """
    manifest_path = write_provenance_manifest(
        output_path=tmp_path / "smoke_test_manifest.json",
        input_paths=[tcr_fixture_path],
        tool_versions={"mimicry_discovery": "0.1.0"},
        extra={"run_type": "ci_smoke_test"},
    )
    assert manifest_path.exists()
