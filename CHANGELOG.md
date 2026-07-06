# Changelog

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
this project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.1.0] - 2026-07-05

Initial scaffold.

### Added
- `mimicry_discovery` package: `io` (10x VDJ + pVACseq parsing, HLA
  nomenclature and reference-sequence loading), `structure` (TCRdock,
  ESMFold, AlphaFold-Multimer adapters), `scoring` (composite mimicry
  risk score), `self_antigen`, `lineage`, and a CLI entry point.
- Snakemake workflow (`workflow/`) with SLURM and Kubernetes executor
  profiles, config-driven and container-aware.
- CI (lint, type-check, unit tests, smoke test), Docker (CPU + CUDA
  targets) with Singularity export, and full FAIR/publication
  scaffolding (`CITATION.cff`, `DATA_GOVERNANCE.md`, `docs/model_card.md`).
- `docs/roadmap.md` and `scripts/seed_github_project.sh` for sprint
  planning.

### Known limitations
See `docs/model_card.md` — structure-prediction backends are real but
incomplete (TCRdock's output parsing is the one remaining step; ESMFold
and AlphaFold-Multimer need model inference wired in), and the
self-antigen reference is not yet built. Tracked in `docs/roadmap.md`.
