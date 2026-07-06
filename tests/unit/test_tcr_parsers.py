"""Unit tests for mimicry_discovery.io.tcr_parsers."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from mimicry_discovery.config import QcThresholds
from mimicry_discovery.io.tcr_parsers import (
    TCRClonotype,
    parse_10x_vdj,
    parse_adaptive_immunoseq,
    repertoire_summary,
)

_COLUMNS = [
    "chain", "v_gene", "j_gene", "cdr3", "productive",
    "high_confidence", "umis", "raw_clonotype_id",
]


def _contig_row(
    chain: str, v_gene: str, j_gene: str, cdr3: str,
    umis: int, clonotype_id: str, productive: str = "True", high_confidence: str = "True",
) -> dict:
    """Build one synthetic 10x contig row with sensible QC-passing defaults."""
    return {
        "chain": chain, "v_gene": v_gene, "j_gene": j_gene, "cdr3": cdr3,
        "productive": productive, "high_confidence": high_confidence,
        "umis": umis, "raw_clonotype_id": clonotype_id,
    }


def _write_10x_csv(path: Path, rows: list[dict]) -> None:
    """Write a minimal synthetic 10x contig annotation CSV for testing."""
    pd.DataFrame(rows, columns=_COLUMNS).to_csv(path, index=False)


def test_parse_10x_vdj_pairs_alpha_and_beta_by_clonotype_id(tmp_path: Path) -> None:
    """Contigs sharing a raw_clonotype_id are paired into one clonotype."""
    csv_path = tmp_path / "contigs.csv"
    _write_10x_csv(csv_path, [
        _contig_row("TRB", "TRBV19", "TRBJ2-7", "CASSIRSSYEQYF", 12, "clonotype4"),
        _contig_row("TRA", "TRAV38-2/DV8", "TRAJ52", "CAYRSAQGGSEKLVF", 10, "clonotype4"),
    ])

    clonotypes = parse_10x_vdj(csv_path, sample_id="pt-001", qc=QcThresholds())

    assert len(clonotypes) == 1
    clone = clonotypes[0]
    assert isinstance(clone, TCRClonotype)
    assert clone.is_paired is True
    assert clone.cdr3_beta == "CASSIRSSYEQYF"
    assert clone.cdr3_alpha == "CAYRSAQGGSEKLVF"
    assert clone.umi_count == 22  # summed across both chains
    assert clone.clonal_frequency == pytest.approx(1.0)


def test_parse_10x_vdj_drops_non_productive_chain_but_keeps_the_other(tmp_path: Path) -> None:
    """A clonotype with one non-productive chain survives as single-chain,
    not dropped entirely."""
    csv_path = tmp_path / "contigs.csv"
    _write_10x_csv(csv_path, [
        _contig_row("TRB", "TRBV7-2", "TRBJ1-1", "CASSPGQGAYEQYF", 3, "clonotype2"),
        _contig_row(
            "TRA", "TRAV12-2", "TRAJ42", "CAVNTGGFKTIF", 1, "clonotype2",
            productive="False",
        ),
    ])

    clonotypes = parse_10x_vdj(csv_path, sample_id="pt-001", qc=QcThresholds())

    assert len(clonotypes) == 1
    clone = clonotypes[0]
    assert clone.is_paired is False
    assert clone.cdr3_beta == "CASSPGQGAYEQYF"
    assert clone.cdr3_alpha is None
    assert clone.umi_count == 3  # only the surviving (productive) contig counted


def test_parse_10x_vdj_missing_file_raises(tmp_path: Path) -> None:
    """A missing input path raises FileNotFoundError, not a pandas error."""
    with pytest.raises(FileNotFoundError):
        parse_10x_vdj(tmp_path / "does_not_exist.csv", sample_id="pt-001")


def test_parse_10x_vdj_missing_columns_raises(tmp_path: Path) -> None:
    """A CSV missing required VDJ columns raises a clear ValueError."""
    bad_csv = tmp_path / "bad.csv"
    pd.DataFrame([{"foo": "bar"}]).to_csv(bad_csv, index=False)
    with pytest.raises(ValueError, match="missing expected 10x VDJ columns"):
        parse_10x_vdj(bad_csv, sample_id="pt-001")


def test_tcr_clonotype_rejects_invalid_cdr3() -> None:
    """Non-amino-acid characters in a CDR3 string are rejected at construction."""
    with pytest.raises(ValueError, match="non-standard amino acid"):
        TCRClonotype(
            sample_id="pt-001", cdr3_alpha=None, v_gene_alpha=None, j_gene_alpha=None,
            cdr3_beta="CASSL123", v_gene_beta="TRBV20-1", j_gene_beta="TRBJ2-1",
            umi_count=5, clonal_frequency=1.0,
        )


def test_tcr_clonotype_requires_at_least_one_chain() -> None:
    """A clonotype with neither chain resolved is rejected -- it means nothing."""
    with pytest.raises(ValueError, match="at least one resolved chain"):
        TCRClonotype(
            sample_id="pt-001", cdr3_alpha=None, v_gene_alpha=None, j_gene_alpha=None,
            cdr3_beta=None, v_gene_beta=None, j_gene_beta=None,
            umi_count=5, clonal_frequency=1.0,
        )


def test_parse_adaptive_immunoseq_is_stubbed() -> None:
    """The Adaptive parser is a documented stub, not a silent no-op."""
    with pytest.raises(NotImplementedError):
        parse_adaptive_immunoseq("dummy.tsv", sample_id="pt-001")


def test_repertoire_summary(tmp_path: Path) -> None:
    """Summary statistics are computed correctly, including the paired count."""
    csv_path = tmp_path / "contigs.csv"
    _write_10x_csv(csv_path, [
        _contig_row("TRB", "TRBV19", "TRBJ2-7", "CASSIRSSYEQYF", 12, "clonotype4"),
        _contig_row("TRA", "TRAV38-2/DV8", "TRAJ52", "CAYRSAQGGSEKLVF", 10, "clonotype4"),
        _contig_row("TRB", "TRBV7-2", "TRBJ1-1", "CASSPGQGAYEQYF", 3, "clonotype2"),
    ])
    clonotypes = parse_10x_vdj(csv_path, sample_id="pt-001")

    summary = repertoire_summary(clonotypes)

    assert summary["n_clonotypes"] == 2
    assert summary["n_paired"] == 1
    assert summary["total_umis"] == 25


def test_repertoire_summary_empty_input() -> None:
    """An empty clonotype list returns zeroed-out summary stats, not an error."""
    summary = repertoire_summary([])
    assert summary == {
        "n_clonotypes": 0, "n_paired": 0,
        "total_umis": 0, "top_clone_frequency": 0.0, "n_unique_v_genes_beta": 0,
    }


def test_parse_10x_vdj_ignores_unassigned_none_clonotype_id(tmp_path: Path) -> None:
    """Contigs with raw_clonotype_id == "None" must not be merged into a
    phantom clonotype, even if they came from otherwise-valid, high-UMI
    contigs of different chains."""
    csv_path = tmp_path / "contigs.csv"
    _write_10x_csv(
        csv_path,
        [
            _contig_row("TRA", "TRAV1-1", "TRAJ1", "CAVRDGGATNKLIF", 50, "None"),
            _contig_row("TRB", "TRBV2", "TRBJ1-2", "CASSEATNTGELFF", 45, "None"),
            _contig_row("TRB", "TRBV19", "TRBJ2-7", "CASSIRSSYEQYF", 12, "clonotype4"),
            _contig_row("TRA", "TRAV38-2/DV8", "TRAJ52", "CAYRSAQGGSEKLVF", 10, "clonotype4"),
        ],
    )

    clonotypes = parse_10x_vdj(csv_path, sample_id="pt-001", qc=QcThresholds())

    assert len(clonotypes) == 1
    clone = clonotypes[0]
    assert clone.cdr3_alpha == "CAYRSAQGGSEKLVF"
    assert clone.cdr3_beta == "CASSIRSSYEQYF"
    assert clone.umi_count == 22
