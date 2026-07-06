# Data Governance

This project processes patient-derived CSF T-cell receptor sequencing and
tumor sequencing/HLA typing data. This document is the short version of
what that means for how data moves through this repository. It is not a
substitute for your institution's IRB protocol, data use agreement (DUA),
or HIPAA/GDPR compliance program — it's the practical summary contributors
need before touching anything in `data/`.

## The one rule

**Real, patient-derived data never enters this git repository — not in a
commit, not in an issue, not in a PR description, not in a screenshot.**
This applies regardless of whether the repository is public or private;
private-but-committed is not the same as properly access-controlled.

## Where real data actually lives

- **Raw sequencing, clinical metadata (`data/raw/`):** stays entirely
  outside git, in your institution's secure/HIPAA-compliant storage.
  `config/samples.tsv` references it only by de-identified sample ID; the
  actual filesystem/bucket path is resolved via an institution-specific
  convention documented separately (ask your data manager), never
  hardcoded into code or config that gets committed.
- **Reference and intermediate data that IS DVC-tracked
  (`data/reference/`, `models/`, some of `data/processed/`):** the `.dvc`
  pointer files (small, hash-only) are committed to git and are safe even
  in a public repo — they reveal that a versioned dataset exists and its
  hash, not its contents. The actual bytes live in a DVC remote with its
  own access control (e.g., an S3 bucket restricted by IAM role to
  approved lab members). `dvc pull` will simply fail for anyone without
  those credentials, including a random person who finds the public repo.
- **Wet-lab results (`wetlab/sample_tracking/`, `wetlab/assays/`):** IDs
  only, matching the same de-identified scheme as `config/samples.tsv`.
  No raw instrument files (e.g. `.fcs` flow cytometry files) belong here —
  those follow the same secure-storage path as raw sequencing.

## Access requirements

Before you can access real data (not the synthetic fixtures in
`data/test/`, which anyone can use freely):

1. Current IRB approval covering this specific analysis use.
2. A signed data use agreement (DUA) if the data originated at a different
   institution than yours.
3. Addition to the approved-personnel list for the relevant DVC remote /
   secure storage bucket.

Contact the PI or the lab manager listed in
[`CODEOWNERS`](.github/CODEOWNERS) to start that process — it is
independent of, and slower than, getting `git` access to this repository.

## De-identification

Every sample is referenced everywhere in this codebase — code, config,
filenames, test fixtures, issue reports — by a de-identified sample ID
only. If you're generating that ID for a new sample, follow your lab's
existing de-identification SOP; don't invent a new scheme locally.

## Before making the repository public

If/when this project is prepared for publication, get explicit sign-off
that:

- No `data/raw/`, real `models/` weights, or real wet-lab results have
  ever been committed (check history, not just the current tree —
  `git log --all -- data/raw` and similar).
- `data/test/` fixtures are confirmed fully synthetic (see
  `data/test/README.md`).
- Your DVC remote's access control is configured correctly for a public
  repo pointing at it (pointer files public, contents still restricted).

This is a human sign-off step, not an automated CI check — CI (see
`.github/workflows/ci.yml`) verifies the code works, not that governance
requirements have been met.
