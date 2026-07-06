#!/usr/bin/env bash
# Thin wrapper: `sbatch`-submits the Snakemake controller itself, which
# then submits each rule as its own SLURM job via the executor plugin.
#
# Usage: ./scripts/submit_slurm.sh [extra snakemake args]

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

sbatch \
  --job-name=mimicry-discovery \
  --partition=compute \
  --time=48:00:00 \
  --mem=4G \
  --wrap="snakemake --configfile config/config.yaml --profile config/profiles/slurm --rerun-incomplete $*"
