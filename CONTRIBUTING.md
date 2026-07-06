# Contributing

This is a mixed computational/wet-lab project. "Contributing" means different
things depending on which side of the lab you sit on — this doc has a path
for both, plus the parts that apply to everyone.

## Everyone, regardless of role

- **Never commit real patient data.** No raw CSF TCR-seq, no tumor sequencing,
  no HLA typing, no anything patient-linkable — not even in a "temporary" test
  file. Use the fixtures in `data/test/` or extend them with more synthetic
  examples. See [`DATA_GOVERNANCE.md`](DATA_GOVERNANCE.md) if you're not sure
  whether something counts.
- **Open an issue before a big change.** Use the right [issue template](.github/ISSUE_TEMPLATE/)
  — `data_issue` if something about the data itself looks wrong (a mismatched
  sample ID, a suspicious clonotype count, an SOP that doesn't match what the
  pipeline expects), `bug_report` / `feature_request` otherwise.
- **Authorship:** if your contribution is substantial enough that you'd expect
  to be a paper co-author, raise that with the PI early — before the PR, not
  after the manuscript is drafted. This repo's commit history and CODEOWNERS
  reviews are not a substitute for that conversation.
- **PRs need the right reviewer.** [`CODEOWNERS`](.github/CODEOWNERS) routes
  `wetlab/**` to wet-lab leads and `mimicry_discovery/**` + `workflow/**` to
  computational leads — this is so a computational reviewer isn't the sole
  approver of an assay protocol change, and vice versa. If your PR touches
  both, expect two reviewers.

## If you're primarily wet-lab

You don't need to be comfortable with Python to contribute meaningfully here.

- **Protocols and SOPs** live in `wetlab/protocols/` as plain Markdown — edit
  them like any shared lab document, then open a PR (ask a computational
  teammate to walk you through `git` once; after that it's copy-paste).
- **Assay results** go through `wetlab/assays/tetramer_results.schema.yaml` —
  that file defines exactly which columns a results file needs. If your
  export doesn't match it, that's useful signal either way: either the schema
  needs updating (open an issue) or the export needs a column renamed.
- **Sample tracking** (`wetlab/sample_tracking/`) should only ever contain
  de-identified IDs — the same IDs used in `config/samples.tsv` on the
  computational side, so results join up correctly without either side
  needing access to the other's identifiers.
- **You don't need to run tests or CI** for a docs/protocol-only change —
  those checks only trigger on code changes.

## If you're primarily dry-lab

Standard software workflow, with a couple of project-specific expectations:

```bash
git clone https://github.com/YOUR-ORG/csf-tcr-mimicry-discovery.git
cd csf-tcr-mimicry-discovery
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pre-commit install   # runs ruff/mypy locally before each commit, matching CI
```

- **Branch naming:** `feat/<short-description>`, `fix/<short-description>`,
  or `data/<short-description>` for anything touching fixtures/reference data.
- **Before opening a PR:** `ruff check . && ruff format --check . && mypy mimicry_discovery && pytest tests/unit tests/integration`
  — this is exactly what CI runs, so a clean local run means a clean CI run.
- **New functions need a docstring (Google style), full type hints, and at
  least one test** — this isn't a style preference, it's load-bearing: this
  codebase gets picked up by people without your context in the room, months
  apart, and the docstrings/tests are how they reconstruct that context.
- **If you're implementing one of the `TODO(...)`-marked stubs** (e.g. an
  actual TCRdock/ESMFold/AlphaFold-Multimer call in
  `mimicry_discovery/structure/`), keep the existing `StructurePredictor`
  interface — the scoring stage and tests depend on that contract, not on
  any particular backend's internals.
- **Changing `data/reference/` or `models/`?** These are DVC-tracked — commit
  the updated `.dvc` pointer file, not the raw asset, and note the source/
  version change in the PR description.

## Commit messages

Conventional-commit-style prefixes (`feat:`, `fix:`, `docs:`, `test:`,
`data:`, `chore:`) — `CHANGELOG.md` generation depends on these being
reasonably consistent.

## Questions

If it's not covered here or in `README.md` / `DATA_GOVERNANCE.md`, open a
discussion or ask in whichever channel your lab already uses — don't guess
on anything involving patient data or authorship.
