"""Snakemake script: parse + QC one sample's TCR-seq and neoantigen calls.

``main()`` is callable directly (e.g. from a test) with plain arguments and
does the real work; only the ``if "snakemake" in globals()`` guard at the
bottom knows about the object Snakemake injects when this file is run as a
rule's ``script:`` target. This split means the actual logic can be
exercised without Snakemake installed at all -- see
``tests/unit/test_workflow_scripts.py``.
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path

from mimicry_discovery.config import QcThresholds
from mimicry_discovery.io.neoantigen_parsers import parse_pvacseq_report
from mimicry_discovery.io.tcr_parsers import parse_10x_vdj


def main(
    tcr_vdj_path: str,
    neoantigen_report_path: str,
    sample_id: str,
    qc: QcThresholds,
    clonotypes_out: str,
    neoantigens_out: str,
) -> None:
    """Parse and QC one sample's raw TCR-seq and neoantigen inputs to JSON.

    Args:
        tcr_vdj_path: Path to that sample's 10x VDJ contig annotation CSV.
        neoantigen_report_path: Path to that sample's pVACseq TSV report.
        sample_id: De-identified sample ID.
        qc: QC thresholds to apply during TCR parsing.
        clonotypes_out: Where to write the parsed, QC'd clonotypes (JSON).
        neoantigens_out: Where to write the parsed neoantigens (JSON).
    """
    clonotypes = parse_10x_vdj(tcr_vdj_path, sample_id=sample_id, qc=qc)
    neoantigens = parse_pvacseq_report(neoantigen_report_path, sample_id=sample_id)

    Path(clonotypes_out).parent.mkdir(parents=True, exist_ok=True)
    Path(clonotypes_out).write_text(
        json.dumps([dataclasses.asdict(c) for c in clonotypes], indent=2)
    )
    Path(neoantigens_out).parent.mkdir(parents=True, exist_ok=True)
    Path(neoantigens_out).write_text(
        json.dumps([dataclasses.asdict(n) for n in neoantigens], indent=2)
    )


if "snakemake" in globals():
    main(
        tcr_vdj_path=snakemake.input.tcr_vdj,  # noqa: F821
        neoantigen_report_path=snakemake.input.neoantigens,  # noqa: F821
        sample_id=snakemake.params.sample_id,  # noqa: F821
        qc=snakemake.params.qc,  # noqa: F821
        clonotypes_out=snakemake.output.clonotypes,  # noqa: F821
        neoantigens_out=snakemake.output.neoantigens,  # noqa: F821
    )
