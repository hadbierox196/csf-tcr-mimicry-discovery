## What does this PR do?

<!-- Brief description. Link any related issue. -->

## Checklist

- [ ] Tests added/updated; `pytest tests/unit tests/integration` passes locally
- [ ] `ruff check . && ruff format --check . && mypy mimicry_discovery` passes locally
- [ ] New/changed functions have docstrings (Google style) and full type hints
- [ ] `CHANGELOG.md` updated, if this is user-facing
- [ ] No real patient data added anywhere in this diff (see `DATA_GOVERNANCE.md`)
- [ ] If touching `data/reference/` or `models/`: the `.dvc` pointer file is
      updated, not the raw asset committed directly
- [ ] If touching `wetlab/`: reviewed by a wet-lab-side CODEOWNER, not only a
      computational one
