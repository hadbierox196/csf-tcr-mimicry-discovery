"""Parsers for tumor neoantigen calls (e.g. pVACseq-style output)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

DEFAULT_PVACSEQ_COLUMN_MAP: dict[str, str] = {
    "peptide": "MT Epitope Seq",
    "wildtype_peptide": "WT Epitope Seq",
    "hla_allele": "HLA Allele",
    "gene_name": "Gene Name",
    "binding_affinity_nm": "Best MT Score",
}
"""Default column mapping for pVACseq's aggregated/filtered TSV report.

TODO(io): pVACseq's exact column names have shifted across versions --
verify this mapping against the header of your actual pVACseq output
before relying on it, and override via the ``column_map`` argument if it
differs.
"""


@dataclass(frozen=True)
class Neoantigen:
    """A single tumor-derived candidate neoantigen-HLA pair.

    Attributes:
        sample_id: De-identified sample/patient identifier.
        peptide: Mutant (tumor-derived) peptide sequence.
        wildtype_peptide: Corresponding wild-type peptide sequence, when
            available, used later for self-similarity comparisons.
        hla_allele: Presenting HLA allele, IPD-IMGT/HLA nomenclature.
        gene_name: Source gene symbol, if available.
        binding_affinity_nm: Predicted peptide-HLA binding affinity in
            nM (lower is stronger); ``None`` if not provided.
    """

    sample_id: str
    peptide: str
    hla_allele: str
    wildtype_peptide: str | None = None
    gene_name: str | None = None
    binding_affinity_nm: float | None = None

    def __post_init__(self) -> None:
        """Validate the peptide sequence and binding affinity.

        Raises:
            ValueError: If ``peptide`` is empty, or
                ``binding_affinity_nm`` is negative.
        """
        if not self.peptide.strip():
            raise ValueError("peptide must not be empty.")
        if self.binding_affinity_nm is not None and self.binding_affinity_nm < 0:
            raise ValueError("binding_affinity_nm must be >= 0 if provided.")


def parse_pvacseq_report(
    report_path: Path | str,
    sample_id: str,
    column_map: dict[str, str] | None = None,
    max_binding_affinity_nm: float | None = 500.0,
) -> list[Neoantigen]:
    """Parse a pVACseq aggregated/filtered TSV report into Neoantigens.

    Args:
        report_path: Path to a pVACseq TSV report.
        sample_id: De-identified sample/patient identifier.
        column_map: Maps this function's field names to the actual
            column headers in ``report_path``. Defaults to
            :data:`DEFAULT_PVACSEQ_COLUMN_MAP` -- override this if your
            pVACseq version uses different headers (see its TODO note).
        max_binding_affinity_nm: Drop calls weaker than this predicted
            affinity (nM). ``None`` disables the filter.

    Returns:
        A list of parsed :class:`Neoantigen` objects.

    Raises:
        FileNotFoundError: If ``report_path`` does not exist.
        ValueError: If a required column from ``column_map`` is
            missing.
    """
    column_map = column_map or DEFAULT_PVACSEQ_COLUMN_MAP
    report_path = Path(report_path)
    if not report_path.exists():
        raise FileNotFoundError(f"pVACseq report not found: {report_path}")

    df = pd.read_csv(report_path, sep="\t")
    required = {column_map["peptide"], column_map["hla_allele"]}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Report is missing required columns {sorted(missing)}; "
            "pass a `column_map` matching your pVACseq version's headers."
        )

    affinity_col = column_map.get("binding_affinity_nm")
    wt_col = column_map.get("wildtype_peptide")
    gene_col = column_map.get("gene_name")

    neoantigens = []
    for _, row in df.iterrows():
        affinity_val: float | None = None
        if affinity_col and affinity_col in df.columns and pd.notna(row[affinity_col]):
            affinity_val = float(row[affinity_col])
        if (
            max_binding_affinity_nm is not None
            and affinity_val is not None
            and affinity_val > max_binding_affinity_nm
        ):
            continue
        neoantigens.append(
            Neoantigen(
                sample_id=sample_id,
                peptide=str(row[column_map["peptide"]]),
                hla_allele=str(row[column_map["hla_allele"]]),
                wildtype_peptide=(
                    str(row[wt_col])
                    if wt_col and wt_col in df.columns and pd.notna(row[wt_col])
                    else None
                ),
                gene_name=(
                    str(row[gene_col])
                    if gene_col and gene_col in df.columns and pd.notna(row[gene_col])
                    else None
                ),
                binding_affinity_nm=affinity_val,
            )
        )
    return neoantigens
