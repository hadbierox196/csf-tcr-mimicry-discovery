"""Unit tests for mimicry_discovery.io.neoantigen_parsers."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from mimicry_discovery.io.neoantigen_parsers import Neoantigen, parse_pvacseq_report


def _write_pvacseq_tsv(path: Path) -> None:
    """Write a minimal synthetic pVACseq-style TSV report for testing."""
    pd.DataFrame(
        [
            {
                "MT Epitope Seq": "SLLMWITQV",
                "WT Epitope Seq": "SLLMWITQC",
                "HLA Allele": "HLA-A*02:01",
                "Gene Name": "TP53",
                "Best MT Score": 45.2,
            },
            {
                "MT Epitope Seq": "GILGFVFTL",
                "WT Epitope Seq": "GILGFVATL",
                "HLA Allele": "HLA-A*02:01",
                "Gene Name": "KRAS",
                "Best MT Score": 800.0,  # weaker than the default 500nM cutoff
            },
        ]
    ).to_csv(path, sep="\t", index=False)


def test_parse_pvacseq_report_filters_by_affinity(tmp_path: Path) -> None:
    """The default 500nM cutoff drops the weak-binding second row."""
    tsv_path = tmp_path / "pvacseq.tsv"
    _write_pvacseq_tsv(tsv_path)

    neoantigens = parse_pvacseq_report(tsv_path, sample_id="pt-001")

    assert len(neoantigens) == 1
    neoantigen = neoantigens[0]
    assert isinstance(neoantigen, Neoantigen)
    assert neoantigen.peptide == "SLLMWITQV"
    assert neoantigen.hla_allele == "HLA-A*02:01"
    assert neoantigen.wildtype_peptide == "SLLMWITQC"
    assert neoantigen.binding_affinity_nm == pytest.approx(45.2)


def test_parse_pvacseq_report_affinity_filter_disabled(tmp_path: Path) -> None:
    """Passing max_binding_affinity_nm=None keeps both rows."""
    tsv_path = tmp_path / "pvacseq.tsv"
    _write_pvacseq_tsv(tsv_path)

    neoantigens = parse_pvacseq_report(tsv_path, sample_id="pt-001", max_binding_affinity_nm=None)

    assert len(neoantigens) == 2


def test_parse_pvacseq_report_missing_file_raises(tmp_path: Path) -> None:
    """A missing report path raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        parse_pvacseq_report(tmp_path / "missing.tsv", sample_id="pt-001")


def test_parse_pvacseq_report_missing_columns_raises(tmp_path: Path) -> None:
    """A report missing required columns raises ValueError, not a KeyError."""
    bad_tsv = tmp_path / "bad.tsv"
    pd.DataFrame([{"foo": "bar"}]).to_csv(bad_tsv, sep="\t", index=False)
    with pytest.raises(ValueError, match="missing required columns"):
        parse_pvacseq_report(bad_tsv, sample_id="pt-001")


def test_neoantigen_rejects_empty_peptide() -> None:
    """Constructing a Neoantigen with an empty peptide raises ValueError."""
    with pytest.raises(ValueError, match="must not be empty"):
        Neoantigen(sample_id="pt-001", peptide="", hla_allele="HLA-A*02:01")
