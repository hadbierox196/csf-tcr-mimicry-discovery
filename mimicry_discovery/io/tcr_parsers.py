"""Parsers for CSF-derived single-cell TCR sequencing data.

Currently supports 10x Genomics Cell Ranger V(D)J output
(``filtered_contig_annotations.csv``). Adaptive Biotechnologies
immunoSEQ support is stubbed -- see :func:`parse_adaptive_immunoseq`.

CORRECTION (roadmap Sprint 2): the original version of this module
treated each chain (TRA, TRB) as an independent "clonotype," which is
wrong for anything that needs the actual receptor -- TCRdock (and,
biologically, T-cell antigen recognition itself) needs the *paired*
alpha+beta combination from a single cell, not either chain in
isolation. This version groups by Cell Ranger's own ``raw_clonotype_id``
(which already encodes that pairing) instead of by chain+CDR3. Found
while wiring in the real TCRdock adapter -- see docs/roadmap.md.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from mimicry_discovery.config import QcThresholds

_VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def _validate_cdr3(value: str | None, field_name: str) -> str | None:
    """Validate an optional CDR3 amino acid string.

    Args:
        value: The CDR3 sequence, or None if that chain wasn't recovered.
        field_name: Field name, used in error messages.

    Returns:
        The validated, upper-cased sequence, or None unchanged.

    Raises:
        ValueError: If ``value`` is an empty string, or contains
            non-standard amino acid codes.
    """
    if value is None:
        return None
    value = value.strip().upper()
    if not value:
        raise ValueError(f"{field_name} must not be an empty string; use None if absent.")
    invalid = set(value) - _VALID_AA
    if invalid:
        raise ValueError(
            f"{field_name} '{value}' contains non-standard amino acid codes: {sorted(invalid)}"
        )
    return value


@dataclass(frozen=True)
class TCRClonotype:
    """A single, QC'd, paired-chain TCR clonotype observed in a CSF sample.

    Represents one clonotype as Cell Ranger defines it: the paired
    alpha+beta receptor shared by one or more cells with identical CDR3s.
    Either chain may be ``None`` if it wasn't recovered for every cell in
    the clonotype (common with single-cell chain dropout) -- but not
    both; a clonotype needs at least one resolved chain to mean anything.

    Attributes:
        sample_id: De-identified sample/patient identifier.
        cdr3_alpha: Alpha-chain CDR3 amino acid sequence, or None if not
            recovered.
        v_gene_alpha: Alpha-chain V-gene segment (IMGT nomenclature).
        j_gene_alpha: Alpha-chain J-gene segment (IMGT nomenclature).
        cdr3_beta: Beta-chain CDR3 amino acid sequence, or None if not
            recovered.
        v_gene_beta: Beta-chain V-gene segment (IMGT nomenclature).
        j_gene_beta: Beta-chain J-gene segment (IMGT nomenclature).
        umi_count: Total UMIs supporting this clonotype (summed across
            all contigs/chains/cells) in the sample.
        clonal_frequency: Proportion of the sample's total UMI-supported
            repertoire made up of this clonotype, in ``[0, 1]``.
    """

    sample_id: str
    cdr3_alpha: str | None
    v_gene_alpha: str | None
    j_gene_alpha: str | None
    cdr3_beta: str | None
    v_gene_beta: str | None
    j_gene_beta: str | None
    umi_count: int
    clonal_frequency: float

    def __post_init__(self) -> None:
        """Validate both CDR3s (if present) and numeric ranges.

        Raises:
            ValueError: If neither chain is present, if a present CDR3
                is invalid, if ``umi_count`` is negative, or if
                ``clonal_frequency`` is outside ``[0, 1]``.
        """
        if self.cdr3_alpha is None and self.cdr3_beta is None:
            raise ValueError("A clonotype needs at least one resolved chain (alpha or beta).")
        object.__setattr__(self, "cdr3_alpha", _validate_cdr3(self.cdr3_alpha, "cdr3_alpha"))
        object.__setattr__(self, "cdr3_beta", _validate_cdr3(self.cdr3_beta, "cdr3_beta"))
        if self.umi_count < 0:
            raise ValueError("umi_count must be >= 0.")
        if not 0.0 <= self.clonal_frequency <= 1.0:
            raise ValueError("clonal_frequency must be in [0, 1].")

    @property
    def is_paired(self) -> bool:
        """Whether both alpha and beta chains were recovered."""
        return self.cdr3_alpha is not None and self.cdr3_beta is not None


def parse_10x_vdj(
    contig_annotations_path: Path | str,
    sample_id: str,
    qc: QcThresholds | None = None,
) -> list[TCRClonotype]:
    """Parse a 10x Genomics Cell Ranger V(D)J contig annotation file.

    Groups contigs by Cell Ranger's own ``raw_clonotype_id`` -- which
    already encodes the paired alpha+beta grouping -- rather than by
    chain and CDR3 independently, so the alpha/beta pairing TCRdock (and
    real biology) needs is preserved. Within a clonotype, if a chain
    has multiple contigs (e.g. dual-alpha), the highest-UMI contig for
    that chain is used to represent it.

    Args:
        contig_annotations_path: Path to a Cell Ranger
            ``filtered_contig_annotations.csv`` file.
        sample_id: De-identified sample/patient identifier to attach to
            every resulting clonotype.
        qc: QC thresholds to apply. Defaults to
            :class:`~mimicry_discovery.config.QcThresholds` defaults if
            not provided.

    Returns:
        A list of :class:`TCRClonotype` objects that passed QC, sorted
        by descending clonal frequency.

    Raises:
        FileNotFoundError: If ``contig_annotations_path`` does not
            exist.
        ValueError: If required columns are missing from the input
            file.
    """
    qc = qc or QcThresholds()
    contig_annotations_path = Path(contig_annotations_path)
    if not contig_annotations_path.exists():
        raise FileNotFoundError(
            f"10x contig annotation file not found: {contig_annotations_path}"
        )

    df = pd.read_csv(contig_annotations_path)
    required_cols = {
        "chain", "v_gene", "j_gene", "cdr3", "productive",
        "high_confidence", "umis", "raw_clonotype_id",
    }
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(
            f"Input file is missing expected 10x VDJ columns: {sorted(missing)}"
        )

    if qc.require_productive:
        df = df[df["productive"].astype(str).str.lower() == "true"]
    if qc.require_high_confidence:
        df = df[df["high_confidence"].astype(str).str.lower() == "true"]
    df = df[df["cdr3"].notna() & (df["cdr3"] != "None")]
    df = df[df["chain"].isin(["TRA", "TRB"])]  # ignore TRG/TRD/multi/etc.
    # Cell Ranger marks contigs it couldn't confidently assign to a
    # clonotype with raw_clonotype_id == "None". Left unfiltered, groupby
    # would merge unrelated cells' contigs into one fictitious clonotype.
    df = df[df["raw_clonotype_id"] != "None"]

    clonotypes: list[TCRClonotype] = []
    for _clonotype_id, group in df.groupby("raw_clonotype_id"):
        chain_repr: dict[str, pd.Series] = {}
        for chain in ("TRA", "TRB"):
            chain_rows = group[group["chain"] == chain]
            if not chain_rows.empty:
                chain_repr[chain] = chain_rows.loc[chain_rows["umis"].idxmax()]

        if not chain_repr:
            continue  # neither chain survived QC for this clonotype

        total_umis = int(group["umis"].sum())
        alpha = chain_repr.get("TRA")
        beta = chain_repr.get("TRB")
        clonotypes.append(
            {
                "sample_id": sample_id,
                "cdr3_alpha": alpha["cdr3"] if alpha is not None else None,
                "v_gene_alpha": alpha["v_gene"] if alpha is not None else None,
                "j_gene_alpha": alpha["j_gene"] if alpha is not None else None,
                "cdr3_beta": beta["cdr3"] if beta is not None else None,
                "v_gene_beta": beta["v_gene"] if beta is not None else None,
                "j_gene_beta": beta["j_gene"] if beta is not None else None,
                "umi_count": total_umis,
            }
        )

    grand_total_umis = sum(c["umi_count"] for c in clonotypes) or 1
    result = []
    for c in clonotypes:
        if c["umi_count"] < qc.min_umi_count:
            continue
        clonal_frequency = c["umi_count"] / grand_total_umis
        if clonal_frequency < qc.min_clonal_frequency:
            continue
        result.append(TCRClonotype(clonal_frequency=clonal_frequency, **c))

    result.sort(key=lambda c: c.clonal_frequency, reverse=True)
    return result


def parse_adaptive_immunoseq(
    export_path: Path | str,
    sample_id: str,
    qc: QcThresholds | None = None,
) -> list[TCRClonotype]:
    """Parse an Adaptive Biotechnologies immunoSEQ export.

    Args:
        export_path: Path to an immunoSEQ TSV export.
        sample_id: De-identified sample/patient identifier.
        qc: QC thresholds to apply.

    Returns:
        A list of :class:`TCRClonotype` objects that passed QC.

    Raises:
        NotImplementedError: Always -- TODO(io): Adaptive's export
            column names vary by assay/export version (typically
            ``aminoAcid``, ``vGeneName``, ``jGeneName``, and a
            templates/reads count column), and standard bulk immunoSEQ
            exports are single-chain (beta-only) with no alpha pairing
            at all -- unlike 10x, there is no clonotype-level pairing
            to recover. Confirm the exact column mapping for your
            institution's current export before implementing, and treat
            the result as beta-only (cdr3_alpha=None for every
            clonotype).
    """
    raise NotImplementedError(
        "Adaptive immunoSEQ parsing is a stub -- confirm your export's "
        "exact column names before implementing (see docstring TODO)."
    )


def repertoire_summary(clonotypes: list[TCRClonotype]) -> dict[str, float | int]:
    """Compute basic repertoire diversity summary statistics.

    Args:
        clonotypes: QC'd clonotypes from a single sample, e.g. the
            output of :func:`parse_10x_vdj`.

    Returns:
        A dict with ``n_clonotypes``, ``n_paired`` (both chains
        recovered), ``total_umis``, ``top_clone_frequency``, and
        ``n_unique_v_genes_beta``.
    """
    if not clonotypes:
        return {
            "n_clonotypes": 0, "n_paired": 0, "total_umis": 0,
            "top_clone_frequency": 0.0, "n_unique_v_genes_beta": 0,
        }
    beta_v_genes = {c.v_gene_beta for c in clonotypes if c.v_gene_beta is not None}
    return {
        "n_clonotypes": len(clonotypes),
        "n_paired": sum(1 for c in clonotypes if c.is_paired),
        "total_umis": sum(c.umi_count for c in clonotypes),
        "top_clone_frequency": max(c.clonal_frequency for c in clonotypes),
        "n_unique_v_genes_beta": len(beta_v_genes),
    }

