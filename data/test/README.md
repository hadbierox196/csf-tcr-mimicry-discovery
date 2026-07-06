# Test fixtures

Small, fully synthetic files used by `tests/integration/` and CI's
smoke-test job. **None of this is real patient data** -- sequences and
IDs are hand-constructed or drawn from public reference examples
(e.g. TP53/KRAS/BRAF hotspot-adjacent peptides) purely to exercise the
parsers and scoring code with realistic shapes.

- `synthetic_tcr_clonotypes.csv` -- 10x Cell Ranger V(D)J-style contig
  annotations (3 barcodes, one non-productive and one low-confidence
  contig included on purpose to exercise QC filtering).
- `synthetic_neoantigens.tsv` -- pVACseq-style aggregated report (3
  candidate neoantigen-HLA pairs; mutant/wild-type peptides are the
  same length by construction, as real point-mutation neoepitopes are).
- `expected_outputs/` -- golden values the integration test asserts
  against. Regenerate deliberately (not by copy-pasting a failing
  test's actual output) if you intentionally change parsing or scoring
  behavior.

Do not add real sequencing output, real HLA typing, or any
patient-linkable data here -- see `DATA_GOVERNANCE.md`.
