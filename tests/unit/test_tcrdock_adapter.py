"""Unit tests for mimicry_discovery.structure.tcrdock_adapter.

Unlike the other two backends (still pure scaffolds), TCRDockPredictor
now does real work -- writing TCRdock's exact targets TSV and invoking
its real two-script subprocess pipeline -- so these tests cover more
than "does it raise NotImplementedError". Actually running TCRdock
itself is out of scope for a unit test (needs a real install, AlphaFold
weights, and a GPU); subprocess calls are mocked to verify the exact
arguments passed match TCRdock's documented CLI.
"""

from __future__ import annotations

import csv
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from mimicry_discovery.structure.base import StructurePredictionRequest
from mimicry_discovery.structure.tcrdock_adapter import TCRDockPredictor, _to_tcrdock_mhc_string

_PAIRED_REQUEST = StructurePredictionRequest(
    tcr_cdr3_beta="CASSIRSSYEQYF",
    tcr_v_gene="TRBV19",
    tcr_j_gene="TRBJ2-7",
    peptide="SLLMWITQV",
    hla_allele="HLA-A*02:01",
    tcr_cdr3_alpha="CAYRSAQGGSEKLVF",
    tcr_v_gene_alpha="TRAV38-2/DV8",
    tcr_j_gene_alpha="TRAJ52",
)
_UNPAIRED_REQUEST = StructurePredictionRequest(
    tcr_cdr3_beta="CASSIRSSYEQYF",
    tcr_v_gene="TRBV19",
    tcr_j_gene="TRBJ2-7",
    peptide="SLLMWITQV",
    hla_allele="HLA-A*02:01",
)


def test_to_tcrdock_mhc_string_strips_hla_prefix() -> None:
    """The HLA- prefix is stripped to match TCRdock's mhc column convention."""
    assert _to_tcrdock_mhc_string("HLA-A*02:01") == "A*02:01"


def test_write_targets_tsv_matches_tcrdock_schema_exactly(tmp_path: Path) -> None:
    """The written TSV has exactly TCRdock's 10 documented columns, correctly filled."""
    predictor = TCRDockPredictor(output_dir=tmp_path)
    targets_path = tmp_path / "targets.tsv"

    predictor._write_targets_tsv(_PAIRED_REQUEST, targets_path)

    with targets_path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert list(row.keys()) == [
        "organism",
        "mhc_class",
        "mhc",
        "peptide",
        "va",
        "ja",
        "cdr3a",
        "vb",
        "jb",
        "cdr3b",
    ]
    assert row["mhc"] == "A*02:01"
    assert row["cdr3a"] == "CAYRSAQGGSEKLVF"
    assert row["cdr3b"] == "CASSIRSSYEQYF"


def test_write_targets_tsv_rejects_unpaired_request(tmp_path: Path) -> None:
    """A beta-only request is rejected before ever touching TCRdock/AlphaFold."""
    predictor = TCRDockPredictor(output_dir=tmp_path)
    with pytest.raises(ValueError, match="paired"):
        predictor._write_targets_tsv(_UNPAIRED_REQUEST, tmp_path / "targets.tsv")


def test_predict_rejects_unpaired_request_without_calling_subprocess(tmp_path: Path) -> None:
    """predict() fails fast on unpaired input -- never shells out for it."""
    predictor = TCRDockPredictor(output_dir=tmp_path)
    with patch("subprocess.run") as mock_run:
        with pytest.raises(ValueError, match="paired"):
            predictor.predict(_UNPAIRED_REQUEST)
        mock_run.assert_not_called()


def test_predict_calls_the_real_two_script_pipeline_in_order(tmp_path: Path) -> None:
    """predict() calls setup_for_alphafold.py then run_prediction.py with the
    exact flags documented in TCRdock's README, then raises NotImplementedError
    (output parsing isn't wired in yet -- see the module docstring)."""
    predictor = TCRDockPredictor(
        output_dir=tmp_path,
        tcrdock_dir="/opt/TCRdock",
        alphafold_data_dir="/data/alphafold",
        model_name="model_2_ptm",
    )
    calls = []

    def fake_run(cmd, check, cwd):
        """Record the subprocess.run call instead of actually invoking anything."""
        calls.append(cmd)
        return MagicMock(returncode=0)

    with patch("subprocess.run", side_effect=fake_run):
        with pytest.raises(NotImplementedError):
            predictor.predict(_PAIRED_REQUEST)

    assert len(calls) == 2
    setup_call, run_call = calls
    assert setup_call[:2] == ["python", "/opt/TCRdock/setup_for_alphafold.py"]
    assert "--targets_tsvfile" in setup_call
    assert "--output_dir" in setup_call
    assert run_call[:2] == ["python", "/opt/TCRdock/run_prediction.py"]
    assert "--targets" in run_call
    assert "--outfile_prefix" in run_call
    assert "model_2_ptm" in run_call
    assert "/data/alphafold" in run_call


def test_predict_raises_file_not_found_if_tcrdock_not_installed(tmp_path: Path) -> None:
    """Without mocking, a non-existent tcrdock_dir surfaces a clear
    FileNotFoundError (not a cryptic subprocess failure) -- this is the
    real error a user sees today, since no environment has TCRdock
    actually installed yet (roadmap Sprint 2, in progress)."""
    predictor = TCRDockPredictor(
        output_dir=tmp_path, tcrdock_dir=tmp_path / "no_such_tcrdock_checkout"
    )
    with pytest.raises(FileNotFoundError):
        predictor.predict(_PAIRED_REQUEST)
