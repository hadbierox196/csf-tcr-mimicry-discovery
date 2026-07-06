"""Unit tests for the stubbed (non-priority) modules.

self_antigen/build_reference.py, cli.py's run_single_sample, and (until
wired in) the structure backends are intentionally incomplete -- tests
here assert the *documented* stub behavior (a clear NotImplementedError)
rather than real functionality. lineage/provenance.py and
self_antigen/query.py are fully implemented despite sitting in this
"stubs tier" module, since they're generic/simple enough to not need a
domain-specific TODO; they get real correctness tests below too.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from mimicry_discovery.cli import build_arg_parser, run_single_sample
from mimicry_discovery.config import PipelineConfig
from mimicry_discovery.lineage.provenance import (
    current_git_sha,
    hash_file,
    write_provenance_manifest,
)
from mimicry_discovery.self_antigen.build_reference import (
    SelfPeptide,
    build_neuronal_self_peptidome,
)
from mimicry_discovery.self_antigen.query import (
    find_self_peptides_by_length,
    find_self_peptides_for_allele,
)


def test_build_neuronal_self_peptidome_is_stubbed(tmp_path: Path) -> None:
    """The reference-building stub raises NotImplementedError, not a silent no-op."""
    with pytest.raises(NotImplementedError):
        build_neuronal_self_peptidome(
            neuronal_expression_reference=tmp_path / "expr.tsv",
            output_path=tmp_path / "out.json",
        )


def test_find_self_peptides_for_allele_filters_correctly() -> None:
    """The (already implemented) allele-query helper filters exactly."""
    reference = [
        SelfPeptide(
            peptide="AAAAAAAAA", hla_allele="HLA-A*02:01",
            source_gene="MBP", evidence="predicted",
        ),
        SelfPeptide(
            peptide="BBBBBBBBB", hla_allele="HLA-A*03:01",
            source_gene="MOG", evidence="predicted",
        ),
        SelfPeptide(
            peptide="CCCCCCC", hla_allele="HLA-A*02:01",
            source_gene="GFAP", evidence="ms_confirmed",
        ),
    ]
    matches = find_self_peptides_for_allele("HLA-A*02:01", reference)
    assert [m.peptide for m in matches] == ["AAAAAAAAA", "CCCCCCC"]


def test_find_self_peptides_by_length_filters_correctly() -> None:
    """The length-query helper returns only peptides of the exact length."""
    reference = [
        SelfPeptide(
            peptide="AAAAAAAAA", hla_allele="HLA-A*02:01",
            source_gene="MBP", evidence="predicted",
        ),
        SelfPeptide(
            peptide="BBBBBBB", hla_allele="HLA-A*02:01",
            source_gene="MOG", evidence="predicted",
        ),
    ]
    matches = find_self_peptides_by_length(9, reference)
    assert [m.peptide for m in matches] == ["AAAAAAAAA"]


def test_find_self_peptides_by_length_rejects_non_positive_length() -> None:
    """A non-positive requested length raises ValueError."""
    with pytest.raises(ValueError, match="must be positive"):
        find_self_peptides_by_length(0, [])


def test_hash_file_is_deterministic(tmp_path: Path) -> None:
    """Hashing the same file content twice yields the same digest."""
    path = tmp_path / "data.txt"
    path.write_text("hello world")
    assert hash_file(path) == hash_file(path)


def test_hash_file_missing_raises(tmp_path: Path) -> None:
    """Hashing a missing file raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        hash_file(tmp_path / "missing.txt")


def test_current_git_sha_outside_repo_returns_none(tmp_path: Path) -> None:
    """A directory with no .git ancestor returns None rather than raising."""
    assert current_git_sha(tmp_path) is None


def test_write_provenance_manifest_contains_expected_keys(tmp_path: Path) -> None:
    """The manifest is written to disk with the documented top-level keys."""
    input_file = tmp_path / "input.csv"
    input_file.write_text("a,b\n1,2\n")

    manifest_path = write_provenance_manifest(
        output_path=tmp_path / "manifest.json",
        input_paths=[input_file],
        tool_versions={"pandas": "2.2.0"},
    )

    assert manifest_path.exists()
    manifest = json.loads(manifest_path.read_text())
    assert set(manifest) == {
        "timestamp_utc", "git_sha", "input_hashes", "tool_versions", "extra",
    }
    assert manifest["tool_versions"] == {"pandas": "2.2.0"}


def test_cli_arg_parser_requires_config_and_sample_id() -> None:
    """The CLI parser rejects invocation missing required arguments."""
    parser = build_arg_parser()
    with pytest.raises(SystemExit):
        parser.parse_args([])


def test_cli_arg_parser_accepts_valid_args() -> None:
    """The CLI parser accepts a well-formed invocation."""
    parser = build_arg_parser()
    args = parser.parse_args(["--config", "config/config.yaml", "--sample-id", "pt-001"])
    assert args.sample_id == "pt-001"
    assert args.config == Path("config/config.yaml")


def test_run_single_sample_is_stubbed(tmp_config: PipelineConfig) -> None:
    """run_single_sample is a documented stub pending structure-backend TODOs."""
    with pytest.raises(NotImplementedError):
        run_single_sample(tmp_config, sample_id="pt-001")
