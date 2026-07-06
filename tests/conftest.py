"""Shared pytest fixtures for the mimicry_discovery test suite."""

from __future__ import annotations

from pathlib import Path

import pytest

from mimicry_discovery.config import PathsConfig, PipelineConfig


@pytest.fixture
def tmp_config(tmp_path: Path) -> PipelineConfig:
    """Build a minimal, valid :class:`PipelineConfig` for tests.

    Args:
        tmp_path: Pytest's built-in temporary directory fixture.

    Returns:
        A validated config pointing at throwaway paths under
        ``tmp_path``, with all other sections left at their defaults.
    """
    return PipelineConfig(
        paths=PathsConfig(
            samples_tsv=tmp_path / "samples.tsv",
            reference_dir=tmp_path / "reference",
            models_dir=tmp_path / "models",
            output_dir=tmp_path / "output",
        )
    )
