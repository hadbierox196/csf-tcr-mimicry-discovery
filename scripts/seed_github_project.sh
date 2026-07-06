#!/usr/bin/env bash
# Seeds the GitHub issue backlog + Project (v2) board for docs/roadmap.md.
#
# Usage:
#   ./scripts/seed_github_project.sh --repo <owner>/<repo>            # dry run (default, safe)
#   ./scripts/seed_github_project.sh --repo <owner>/<repo> --apply    # actually create things
#
# Requires the GitHub CLI, authenticated with the `project` scope:
#   gh auth login
#   gh auth refresh -s project
#
# Re-running is mostly safe: `gh label create --force` overwrites, and
# milestone/issue creation will just add duplicates if run twice with
# --apply -- this is a one-time bootstrap script, not idempotent
# infrastructure-as-code, so review the dry-run output before --apply.

set -euo pipefail

REPO=""
APPLY=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo) REPO="$2"; shift 2 ;;
    --apply) APPLY=true; shift ;;
    -h|--help)
      echo "Usage: $0 --repo <owner>/<repo> [--apply]"
      exit 0
      ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

if [[ -z "$REPO" ]]; then
  echo "Usage: $0 --repo <owner>/<repo> [--apply]" >&2
  exit 1
fi
OWNER="${REPO%%/*}"

if ! command -v gh >/dev/null 2>&1; then
  echo "ERROR: GitHub CLI ('gh') not found. Install: https://cli.github.com" >&2
  exit 1
fi
if ! gh auth status >/dev/null 2>&1; then
  echo "ERROR: not authenticated. Run: gh auth login && gh auth refresh -s project" >&2
  exit 1
fi

if $APPLY; then
  echo "=== APPLY mode: this will create real issues/milestones/a project on ${REPO} ==="
else
  echo "=== DRY RUN (default) -- pass --apply to actually create anything ==="
fi
echo ""

run() {
  if $APPLY; then
    echo "+ $*"
    "$@"
  else
    echo "[dry run] $*"
  fi
}

# ---------------------------------------------------------------------------
# 1. Labels
# ---------------------------------------------------------------------------
run gh label create "roadmap" --repo "$REPO" --color "5319E7" \
  --description "Tracked in docs/roadmap.md" --force
for area in "area:orchestration" "area:structure-prediction" "area:self-antigen" \
            "area:wetlab" "area:io" "area:scoring" "area:containers" "area:docs"; do
  run gh label create "$area" --repo "$REPO" --color "1D76DB" --force
done

# ---------------------------------------------------------------------------
# 2. Milestones (one per sprint)
# ---------------------------------------------------------------------------
declare -a MILESTONES=(
  "Sprint 1: Orchestration foundation"
  "Sprint 2: Structure prediction backend #1"
  "Sprint 3: Structure prediction backends #2-3 + self-antigen"
  "Sprint 4: Wet-lab validation loop"
  "Sprint 5: Scale-out and publication readiness"
)
for m in "${MILESTONES[@]}"; do
  run gh api "repos/${REPO}/milestones" -f title="$m" --silent
done

# ---------------------------------------------------------------------------
# 3. Issues -- one explicit call per issue (verbose on purpose: easier to
#    read, diff, and edit individual issues than a clever data-driven loop).
# ---------------------------------------------------------------------------
create_issue() {
  local title="$1" milestone="$2" labels="$3" body="$4"
  run gh issue create --repo "$REPO" --title "$title" --milestone "$milestone" \
    --label "roadmap,${labels}" --body "$body"
}

create_issue \
  "Implement workflow/Snakefile + rules/*.smk" \
  "Sprint 1: Orchestration foundation" "area:orchestration" \
  "5-stage DAG (ingest_qc, candidate_pairing, structure_prediction, risk_scoring, wetlab_export) wired to the existing mimicry_discovery package. See docs/roadmap.md #1."

create_issue \
  "Set up DVC remote and convert data/reference + models to real .dvc pointers" \
  "Sprint 1: Orchestration foundation" "area:orchestration" \
  "Point .dvc/config at real institutional storage. See architecture rationale section 2.2 and docs/roadmap.md #2."

create_issue \
  "Implement real Snakemake profiles for SLURM and cloud" \
  "Sprint 1: Orchestration foundation" "area:orchestration" \
  "config/profiles/{slurm,cloud}/ currently only described, not implemented. See architecture rationale section 2.3 and docs/roadmap.md #3."

create_issue \
  "Implement TCRdock output parsing (predict() is otherwise real)" \
  "Sprint 2: Structure prediction backend #1" "area:structure-prediction" \
  "TCRDockPredictor.predict() now runs the real setup_for_alphafold.py + run_prediction.py subprocess pipeline with verified arguments -- only turning the resulting PDB+PAE output into a StructurePredictionResult remains (TCRdock's exact output file naming wasn't confirmable without a real install). See the TODO in tcrdock_adapter.py and docs/roadmap.md Sprint 2."

create_issue \
  "Implement ESMFoldPredictor.predict() (HLA lookup dependency now resolved)" \
  "Sprint 3: Structure prediction backends #2-3 + self-antigen" "area:structure-prediction" \
  "Only model loading/inference remain -- the HLA sequence lookup this needs is done (mimicry_discovery/io/hla_sequences.py, wired via ESMFoldPredictor.from_reference_fasta()). Remove the TODO(structure-prediction) in esmfold_adapter.py. See docs/roadmap.md #7."

create_issue \
  "Implement AlphaFoldMultimerPredictor.predict()" \
  "Sprint 3: Structure prediction backends #2-3 + self-antigen" "area:structure-prediction" \
  "Remove the TODO(structure-prediction) in af_multimer_adapter.py. See docs/roadmap.md #8."

create_issue \
  "Implement build_neuronal_self_peptidome()" \
  "Sprint 3: Structure prediction backends #2-3 + self-antigen" "area:self-antigen" \
  "Source in-house immunopeptidomics MS data or a public resource, filtered to CNS-expressed genes. Remove the TODO(self-antigen) in self_antigen/build_reference.py. See docs/roadmap.md #9."

create_issue \
  "Build out the wetlab/ module" \
  "Sprint 4: Wet-lab validation loop" "area:wetlab" \
  "protocols/tetramer_staining_sop.md, finalized assays/tetramer_results.schema.yaml against a real export, sample_tracking/, and validation_reports/build_concordance_report.py. See docs/roadmap.md #10."

create_issue \
  "Wire scoring/calibration.py into a real calibration workflow" \
  "Sprint 4: Wet-lab validation loop" "area:scoring" \
  "Once the first real tetramer results exist. See docs/model_card.md 'Known limitations' and docs/roadmap.md #11."

create_issue \
  "Implement cli.run_single_sample() end-to-end" \
  "Sprint 4: Wet-lab validation loop" "area:io" \
  "Depends on the structure-prediction backends being real. Remove the TODO(cli) in mimicry_discovery/cli.py. See docs/roadmap.md #12."

create_issue \
  "Implement the Adaptive immunoSEQ parser" \
  "Sprint 5: Scale-out and publication readiness" "area:io" \
  "Remove the TODO(io) in mimicry_discovery/io/tcr_parsers.py -- confirm exact column names against your institution's current export first. See docs/roadmap.md #13."

create_issue \
  "Upgrade peptide_similarity() to BLOSUM62 via Biopython" \
  "Sprint 5: Scale-out and publication readiness" "area:scoring" \
  "Replaces the property-group heuristic. Remove the TODO(scoring) in mimicry_discovery/scoring/features.py. See docs/roadmap.md #14."

create_issue \
  "Set up container registry + publish workflow" \
  "Sprint 5: Scale-out and publication readiness" "area:containers" \
  "build-containers.yml plus pointing containers/singularity/build_sif.sh's REGISTRY at the real path. See docs/roadmap.md #15."

create_issue \
  "MkDocs docs site + Zenodo DOI-per-release integration" \
  "Sprint 5: Scale-out and publication readiness" "area:docs" \
  "Publish docs/ as a site; wire CITATION.cff to Zenodo so tagged releases mint a DOI. See docs/roadmap.md #16."

# ---------------------------------------------------------------------------
# 4. Project (v2) board, with every issue just created added to it
# ---------------------------------------------------------------------------
PROJECT_TITLE="CSF TCR Mimicry Discovery Roadmap"

if $APPLY; then
  echo ""
  echo "+ gh project create --owner ${OWNER} --title \"${PROJECT_TITLE}\""
  PROJECT_URL=$(gh project create --owner "$OWNER" --title "$PROJECT_TITLE" --format json | python3 -c "import json,sys; print(json.load(sys.stdin)['url'])")
  echo "Created: ${PROJECT_URL}"

  echo "Adding every roadmap issue to the project..."
  gh issue list --repo "$REPO" --label roadmap --state open --limit 100 --json url \
    | python3 -c "import json,sys; [print(i['url']) for i in json.load(sys.stdin)]" \
    | while read -r issue_url; do
        gh project item-add --owner "$OWNER" --url "$issue_url" "$(echo "$PROJECT_URL" | grep -oE '[0-9]+$')" >/dev/null
      done
  echo "Done."
else
  echo ""
  echo "[dry run] gh project create --owner ${OWNER} --title \"${PROJECT_TITLE}\""
  echo "[dry run] (then) add every issue labeled 'roadmap' to that project"
fi

echo ""
if $APPLY; then
  echo "=== Done. 5 milestones, 16 issues, 1 project board created on ${REPO}. ==="
else
  echo "=== Dry run complete. Re-run with --apply to actually create these. ==="
fi
