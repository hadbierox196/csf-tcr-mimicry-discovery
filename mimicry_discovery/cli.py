"""Command-line entry point for the mimicry-discovery pipeline.

Thin wiring only -- real stage logic lives in
:mod:`mimicry_discovery.io`, :mod:`mimicry_discovery.structure`, and
:mod:`mimicry_discovery.scoring`. This CLI exists mainly for ad-hoc,
single-sample runs outside Snakemake (e.g. local debugging); the
Snakemake workflow is the source of truth for full-cohort runs.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from mimicry_discovery.config import PipelineConfig, load_config


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the top-level CLI argument parser.

    Returns:
        A configured :class:`argparse.ArgumentParser`.
    """
    parser = argparse.ArgumentParser(
        prog="mimicry-discovery",
        description="Structure-based TCR-pMHC molecular mimicry discovery.",
    )
    parser.add_argument(
        "--config", type=Path, required=True,
        help="Path to a validated config.yaml (see config/config.yaml).",
    )
    parser.add_argument(
        "--sample-id", type=str, required=True,
        help="De-identified sample/patient ID to process.",
    )
    return parser


def run_single_sample(config: PipelineConfig, sample_id: str) -> None:
    """Run the full pipeline for one sample outside of Snakemake.

    Args:
        config: A validated pipeline configuration.
        sample_id: De-identified sample/patient ID to process.

    Raises:
        NotImplementedError: Always -- TODO(cli): wire together
            io -> structure -> scoring for a single sample once the
            structure-prediction backends are implemented (see the
            TODOs in mimicry_discovery/structure/*_adapter.py).
    """
    raise NotImplementedError(
        f"Single-sample run for '{sample_id}' is a stub pending the "
        "structure-prediction backend TODOs."
    )


def main(argv: list[str] | None = None) -> int:
    """CLI entry point (wired to ``mimicry-discovery`` in pyproject.toml).

    Args:
        argv: Argument list to parse; defaults to ``sys.argv[1:]``.

    Returns:
        Process exit code (0 on success).
    """
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    config = load_config(args.config)
    run_single_sample(config, args.sample_id)
    return 0


if __name__ == "__main__":
    sys.exit(main())
