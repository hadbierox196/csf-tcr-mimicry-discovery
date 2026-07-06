"""HLA allele nomenclature parsing and normalization.

Implements IPD-IMGT/HLA nomenclature (e.g. ``HLA-A*02:01:01:01``), the
standard maintained by the WHO Nomenclature Committee for Factors of
the HLA System.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

# HLA-<gene>*<field1>:<field2>[:<field3>[:<field4>]][<suffix>]
# e.g. HLA-A*02:01, HLA-DRB1*15:01:01, HLA-B*57:01N
_HLA_PATTERN = re.compile(
    r"^HLA-(?P<gene>[A-Z0-9]+)\*"
    r"(?P<field1>\d{2,3}):(?P<field2>\d{2,3})"
    r"(?::(?P<field3>\d{2,3}))?"
    r"(?::(?P<field4>\d{2,3}))?"
    r"(?P<suffix>[NLSCAQ])?$"
)


@dataclass(frozen=True)
class HLAAllele:
    """A single, normalized HLA allele call.

    Attributes:
        gene: HLA gene name, e.g. ``"A"``, ``"B"``, ``"DRB1"``.
        field1: Allele group (first field), e.g. ``"02"`` in ``A*02:01``.
        field2: Specific HLA protein (second field), e.g. ``"01"``.
        field3: Synonymous DNA substitution field, if resolved.
        field4: Non-coding region variation field, if resolved.
        suffix: Expression-status suffix (N/L/S/C/A/Q), if present.
    """

    gene: str
    field1: str
    field2: str
    field3: str | None = None
    field4: str | None = None
    suffix: str | None = None

    @property
    def two_field(self) -> str:
        """Return the canonical two-field allele string, e.g. ``"HLA-A*02:01"``.

        Two-field resolution is what most peptide-HLA binding predictors
        expect; higher-field information is kept on the object but not
        included here.
        """
        return f"HLA-{self.gene}*{self.field1}:{self.field2}"


def parse_hla_allele(raw: str) -> HLAAllele:
    """Parse and validate a raw HLA allele string into an HLAAllele.

    Accepts the standard ``HLA-``-prefixed form and also tolerates a
    missing ``HLA-`` prefix (e.g. ``"A*02:01"``) or embedded whitespace
    (as sometimes appears in clinical typing reports), normalizing both
    to the standard form before validation.

    Args:
        raw: A raw allele string, e.g. ``"HLA-A*02:01"``, ``"A*02:01"``,
            or ``"A* 02:01"``.

    Returns:
        A validated :class:`HLAAllele`.

    Raises:
        ValueError: If ``raw`` does not match IPD-IMGT/HLA nomenclature
            after normalization.
    """
    normalized = raw.strip().upper().replace(" ", "")
    if not normalized.startswith("HLA-"):
        normalized = f"HLA-{normalized}"

    match = _HLA_PATTERN.match(normalized)
    if match is None:
        raise ValueError(
            f"'{raw}' is not a valid IPD-IMGT/HLA allele string "
            f"(expected e.g. 'HLA-A*02:01')."
        )
    groups = match.groupdict()
    return HLAAllele(
        gene=groups["gene"],
        field1=groups["field1"],
        field2=groups["field2"],
        field3=groups.get("field3"),
        field4=groups.get("field4"),
        suffix=groups.get("suffix"),
    )


def shares_supertype(allele_a: HLAAllele, allele_b: HLAAllele) -> bool:
    """Check whether two HLA alleles share the same gene and first field.

    This is a coarse proxy for peptide-binding-groove similarity, used
    as an early, cheap filter before running expensive structure
    prediction on a candidate self-peptide presented on a different
    HLA than the tumor peptide it's being compared against (see
    :mod:`mimicry_discovery.self_antigen.query`).

    Args:
        allele_a: First allele.
        allele_b: Second allele.

    Returns:
        True if both alleles share the same gene and first-field group.
    """
    return allele_a.gene == allele_b.gene and allele_a.field1 == allele_b.field1
