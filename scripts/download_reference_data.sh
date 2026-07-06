#!/usr/bin/env bash
# Downloads reference data into data/reference/. Currently covers the
# IPD-IMGT/HLA protein sequences (roadmap Sprint 2); IMGT germline
# V/J/C genes, PDB TCR-pMHC templates, and the neuronal self-peptidome
# are tracked separately (see docs/roadmap.md Sprint 3) and not yet
# wired into this script.
#
# Usage: ./scripts/download_reference_data.sh
#
# Source verified directly (fetched a real file to confirm the header
# format -- see mimicry_discovery/io/hla_sequences.py's module
# docstring) rather than assumed: https://github.com/ANHIG/IMGTHLA

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

DEST="data/reference/hla_alleles"
mkdir -p "${DEST}"

echo "Downloading IPD-IMGT/HLA class I protein sequences to ${DEST}/..."
for locus in A B C; do
  curl -fsSL \
    "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/fasta/${locus}_prot.fasta" \
    -o "${DEST}/${locus}_prot.fasta"
  echo "  ${locus}_prot.fasta: $(grep -c '^>' "${DEST}/${locus}_prot.fasta") alleles"
done

echo ""
echo "Done. Point ESMFoldPredictor.from_reference_fasta() at, e.g.,"
echo "  ${DEST}/A_prot.fasta"
echo "(the current default reference_dir/hla_alleles/A_prot.fasta path"
echo "used by workflow/rules/03_structure_prediction.smk)."
