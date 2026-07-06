# Model Card: Molecular Mimicry Risk Score

## Intended use

**Research candidate-prioritization only.** The risk score produced by
`mimicry_discovery.scoring` ranks candidate TCR / tumor-peptide /
neuronal-self-peptide triples so a wet-lab team can prioritize which ones
are worth an actual tetramer-binding assay. It is a hypothesis-generation
tool for an active research pipeline.

**It is not a validated clinical diagnostic or screening test.** It has
not been through clinical validation, is not FDA-cleared or equivalent,
and should not inform a clinical decision about any individual patient.
"Pre-symptomatic screening" in the project's stated long-term goal refers
to a future validated version of this work, not the current state of this
repository.

## How the score is currently computed

A transparent, config-weighted linear combination of three features (see
`mimicry_discovery/scoring/risk_model.py` and
`mimicry_discovery/scoring/features.py`):

1. **Peptide similarity** between the tumor-derived peptide and a
   candidate neuronal self-peptide (currently a physicochemical
   property-group heuristic — see the `TODO` in `scoring/features.py` for
   the planned BLOSUM62 upgrade).
2. **HLA anchor-residue conservation** between the same pair.
3. **Structural interface confidence** from a TCR:self-peptide-HLA
   complex prediction (TCRdock / ESMFold / AlphaFold-Multimer — these
   backends are currently scaffolded; see
   `mimicry_discovery/structure/*_adapter.py`).

There is no learned/trained component yet — weights in
`config/config.yaml` are analyst-set defaults, not fit to labeled data.
`mimicry_discovery/scoring/calibration.py` exists specifically so that,
once tetramer-confirmed positive/negative pairs accumulate via
`wetlab/validation_reports/`, the score can be calibrated (and eventually
replaced with a learned model) against real ground truth.

## Known limitations

- **No structural signal is currently real** — until the structure
  backends are implemented, any score computed today substitutes a
  placeholder in place of true structural confidence (see the smoke test
  in `tests/integration/` for exactly how that placeholder is used).
- **No calibration data yet.** Precision/recall/AUROC against real
  tetramer results (`scoring/calibration.py`) are not yet computed
  because no such results exist yet in this project's timeline.
- **Small expected sample sizes.** Paraneoplastic neurodegeneration is
  rare; even a successful version of this pipeline should be expected to
  operate on cohorts too small for the usual large-N ML validation
  standards, and results should be interpreted with that in mind.
- **Self-peptide reference is not yet built** — the neuronal
  self-peptidome reference this comparison depends on
  (`mimicry_discovery/self_antigen/build_reference.py`) is itself a
  documented stub.

## Who to contact

Questions about the current validation status of this score: the PI or
computational lead listed in [`CODEOWNERS`](.github/CODEOWNERS).
