#!/usr/bin/env bash
# Converts a published mimicry-discovery Docker image to a Singularity/
# Apptainer .sif image for HPC nodes that don't allow rootful Docker.
#
# Usage:
#   ./build_sif.sh <scoring|structure-prediction> [tag]
#
# Examples:
#   ./build_sif.sh scoring                  # pulls :scoring-latest
#   ./build_sif.sh structure-prediction v0.1.0
#
# TODO(containers): replace REGISTRY below with your actual GHCR/registry
# path once container publishing is set up (tracked separately from this
# repo's CI, which currently only lints/type-checks/tests/smoke-tests --
# see .github/workflows/ci.yml).

set -euo pipefail

REGISTRY="ghcr.io/YOUR-ORG/mimicry-discovery"

TARGET="${1:?Usage: build_sif.sh <scoring|structure-prediction> [tag]}"
TAG="${2:-latest}"

case "${TARGET}" in
  scoring|structure-prediction) ;;
  *)
    echo "ERROR: unknown target '${TARGET}' (expected 'scoring' or 'structure-prediction')" >&2
    exit 1
    ;;
esac

IMAGE="${REGISTRY}:${TARGET}-${TAG}"
OUT="mimicry-discovery-${TARGET}-${TAG}.sif"

if command -v apptainer >/dev/null 2>&1; then
  BUILDER="apptainer"
elif command -v singularity >/dev/null 2>&1; then
  BUILDER="singularity"
else
  echo "ERROR: neither 'apptainer' nor 'singularity' found on PATH." >&2
  echo "On most clusters: module load apptainer   (or: module load singularity)" >&2
  exit 1
fi

echo "Building ${OUT} from docker://${IMAGE} using ${BUILDER}..."
"${BUILDER}" build "${OUT}" "docker://${IMAGE}"

echo ""
echo "Done: ${OUT}"
if [ "${TARGET}" = "structure-prediction" ]; then
  echo "Run with GPU access:  ${BUILDER} exec --nv ${OUT} mimicry-discovery --help"
else
  echo "Run with:             ${BUILDER} exec ${OUT} mimicry-discovery --help"
fi
echo "Bind real paths at run time, e.g.:"
echo "  ${BUILDER} exec --bind /scratch/\$USER/data:/opt/mimicry-discovery/data ${OUT} ..."
