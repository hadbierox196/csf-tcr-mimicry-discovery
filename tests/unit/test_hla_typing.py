"""Unit tests for mimicry_discovery.io.hla_typing."""

from __future__ import annotations

import pytest

from mimicry_discovery.io.hla_typing import parse_hla_allele, shares_supertype


def test_parse_hla_allele_normalizes_missing_prefix() -> None:
    """A prefix-less allele string parses to the same result as HLA--prefixed."""
    with_prefix = parse_hla_allele("HLA-A*02:01")
    without_prefix = parse_hla_allele("A*02:01")
    assert with_prefix == without_prefix
    assert with_prefix.two_field == "HLA-A*02:01"


def test_parse_hla_allele_tolerates_embedded_whitespace() -> None:
    """A space between the gene and field block (as seen in some clinical
    reports) is stripped before matching."""
    assert parse_hla_allele("A* 02:01") == parse_hla_allele("HLA-A*02:01")


def test_parse_hla_allele_handles_class_ii_and_extra_fields() -> None:
    """A four-part class II allele resolves field3 correctly."""
    allele = parse_hla_allele("HLA-DRB1*15:01:01")
    assert allele.gene == "DRB1"
    assert allele.field1 == "15"
    assert allele.field2 == "01"
    assert allele.field3 == "01"


def test_parse_hla_allele_handles_suffix() -> None:
    """A trailing expression-status suffix (e.g. N for null) is captured."""
    allele = parse_hla_allele("HLA-B*57:01N")
    assert allele.suffix == "N"


def test_parse_hla_allele_rejects_malformed_string() -> None:
    """A clearly malformed allele string raises ValueError, not a silent parse."""
    with pytest.raises(ValueError, match="not a valid"):
        parse_hla_allele("not-an-allele")


def test_shares_supertype() -> None:
    """Two alleles sharing gene+field1 are flagged; a different gene is not."""
    a02_01 = parse_hla_allele("HLA-A*02:01")
    a02_07 = parse_hla_allele("HLA-A*02:07")
    b07_02 = parse_hla_allele("HLA-B*07:02")
    assert shares_supertype(a02_01, a02_07) is True
    assert shares_supertype(a02_01, b07_02) is False
