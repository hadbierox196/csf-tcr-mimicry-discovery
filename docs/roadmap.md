# Roadmap

Every item below traces to an actual `TODO(...)` or an explicitly-scoped-but-
not-yet-built piece of the architecture from earlier phases — this isn't a
generic template backlog, it's what's actually left. `scripts/seed_github_project.sh`
turns this document into real GitHub milestones, issues, and a Project board.

Five two-week sprints, ordered by dependency (structure prediction before
self-antigen search makes sense of it; wet-lab loop needs real candidates
to validate first; polish/scale is last on purpose).

## Sprint 1 — Orchestration foundation

| # | Issue | Traces to |
|---|---|---|
| 1 | Write `workflow/Snakefile` + `workflow/rules/*.smk` (5-stage DAG: ingest_qc, candidate_pairing, structure_prediction, risk_scoring, wetlab_export), wired to the existing `mimicry_discovery` package | Scoped in the original architecture; never implemented |
| 2 | Set up the DVC remote (`.dvc/config`) against real institutional storage; convert `data/reference/*` and `models/` from plain directories to real `.dvc` pointer files | Architecture rationale, §2.2 (data lineage) |
| 3 | Implement `config/profiles/slurm/` and `config/profiles/cloud/` as real Snakemake executor profiles | Architecture rationale, §2.3 (HPC/cloud divergence) |

## Sprint 2 — Structure prediction backend #1 ✅ mostly done

| # | Issue | Status |
|---|---|---|
| 4 | Implement `TCRDockPredictor.predict()` | ✅ Real subprocess pipeline (verified against TCRdock's actual README) writes the exact targets TSV and calls `setup_for_alphafold.py` + `run_prediction.py`. **Still open:** output parsing (turning TCRdock's PDB+PAE output into a `StructurePredictionResult`) — TCRdock's output-file naming wasn't fetchable to confirm, see the TODO in `tcrdock_adapter.py`. |
| 5 | Build a real IPD-IMGT/HLA protein sequence lookup | ✅ `mimicry_discovery/io/hla_sequences.py` parses the real `ANHIG/IMGTHLA` FASTA format (verified against an actually-fetched file), wired into `ESMFoldPredictor.from_reference_fasta()`. `scripts/download_reference_data.sh` fetches the real files. |
| 6 | Uncomment/pin TCRdock in the structure-prediction Docker target | ✅ Clones the real TCRdock repo. AlphaFold's own setup (genetic DBs, params) is intentionally left to TCRdock's own README rather than duplicated here. |
| — | *(unplanned, found during #4)* Fix TCR alpha/beta pairing | ✅ `TCRClonotype` and `parse_10x_vdj` originally treated each chain as an independent clonotype — wrong for TCRdock, which needs the paired receptor. Reworked to group by Cell Ranger's own `raw_clonotype_id`. See the correction note in `tcr_parsers.py`'s module docstring. |

**What's NOT done and shouldn't be assumed working:** TCRdock's output parsing (above), and therefore an actual end-to-end real structure prediction. Running `TCRDockPredictor.predict()` today gets you two genuine subprocess calls with correct arguments, then either `FileNotFoundError` (TCRdock/AlphaFold not installed) or, once they are, a `NotImplementedError` at the parsing step.

## Sprint 3 — Structure prediction backends #2/#3 + self-antigen

| # | Issue | Traces to |
|---|---|---|
| 7 | Implement `ESMFoldPredictor.predict()` (HLA lookup dependency from Sprint 2 is resolved) | `TODO(structure-prediction)` in `esmfold_adapter.py` — only model loading/inference remain |
| 8 | Implement `AlphaFoldMultimerPredictor.predict()` | `TODO(structure-prediction)` in `af_multimer_adapter.py` |
| 9 | Implement `build_neuronal_self_peptidome()` — source in-house immunopeptidomics MS data or a public resource, filtered to CNS-expressed genes | `TODO(self-antigen)` in `self_antigen/build_reference.py` |

## Sprint 4 — Wet-lab validation loop

| # | Issue | Traces to |
|---|---|---|
| 10 | Build out `wetlab/`: `protocols/tetramer_staining_sop.md`, finalize `assays/tetramer_results.schema.yaml` against a real export, `sample_tracking/`, and `validation_reports/build_concordance_report.py` | Scoped in the original architecture; never implemented |
| 11 | Wire `scoring/calibration.py` into an actual calibration workflow once the first real tetramer results exist | `docs/model_card.md`, "Known limitations" |
| 12 | Implement `cli.run_single_sample()` end-to-end now that structure prediction is real (depends on #4/#7/#8) | `TODO(cli)` in `mimicry_discovery/cli.py` |

## Sprint 5 — Scale-out and publication readiness

| # | Issue | Traces to |
|---|---|---|
| 13 | Implement the Adaptive immunoSEQ parser | `TODO(io)` in `mimicry_discovery/io/tcr_parsers.py` |
| 14 | Upgrade `peptide_similarity()` from the property-group heuristic to BLOSUM62 via Biopython | `TODO(scoring)` in `mimicry_discovery/scoring/features.py` |
| 15 | Set up a container registry + `build-containers.yml` publish workflow; point `containers/singularity/build_sif.sh`'s `REGISTRY` at the real path | `TODO(containers)` in `build_sif.sh` |
| 16 | MkDocs docs site (publish `docs/`) + Zenodo DOI-per-release integration | README "Citation" section; original architecture §3 (GitHub settings) |

## Using this with the seed script

```bash
gh auth login
gh auth refresh -s project     # milestones/issues use default scopes; the
                                 # Project board needs this one specifically
./scripts/seed_github_project.sh --repo YOUR-ORG/csf-tcr-mimicry-discovery
# review the dry-run output, then:
./scripts/seed_github_project.sh --repo YOUR-ORG/csf-tcr-mimicry-discovery --apply
```

This creates 5 milestones, 16 labeled issues, a Project (v2) board titled
"CSF TCR Mimicry Discovery Roadmap," and adds every issue to it. Status
columns (Todo/In Progress/Done) come from the board's default view — no
extra setup needed unless you want to customize them.
