#!/usr/bin/env bash
# Thin wrapper around the cloud (Kubernetes) profile. Run from a machine
# with kubectl access to the target cluster and the storage plugin's
# credentials (e.g. AWS creds for the S3 storage provider) available.
#
# Usage: ./scripts/submit_cloud.sh [extra snakemake args]

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

snakemake \
  --configfile config/config.yaml \
  --profile config/profiles/cloud \
  --rerun-incomplete \
  "$@"
