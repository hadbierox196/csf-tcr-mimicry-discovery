"""Snakemake script: predict TCR:self-peptide-HLA complex structures.

Real backend calls (TCRdock/ESMFold/AlphaFold-Multimer) are still
scaffolded -- see the TODOs in ``mimicry_discovery/structure/*_adapter.py``
(roadmap Sprints 2-3). Invoking this rule normally surfaces that
``NotImplementedError``, which is correct: this stage genuinely can't
produce real structures yet.

For smoke-testing the rest of the DAG before those backends land, pass
``--config smoke_test=true`` on the ``snakemake`` command line -- this
substitutes a fixed placeholder confidence instead of calling a real
backend, exactly mirroring what
``tests/integration/test_pipeline_end_to_end.py`` already does at the
package level. It is never the default.

Note: if ``structure_prediction.backend`` is set to ``tcrdock`` and a
candidate's clonotype has no resolved alpha chain, TCRDockPredictor
raises ValueError rather than silently modeling beta alone (TCRdock's
input schema has no way to represent a single-chain target) -- route
beta-only clonotypes through the ``esmfold`` backend instead, or filter
them out upstream if TCRdock-only analysis is intended.
"""

from __future__ import annotations

import json
from pathlib import Path

from mimicry_discovery.config import StructureBackend, StructurePredictionConfig
from mimicry_discovery.structure.af_multimer_adapter import AlphaFoldMultimerPredictor
from mimicry_discovery.structure.base import StructurePredictionRequest, StructurePredictor
from mimicry_discovery.structure.esmfold_adapter import ESMFoldPredictor
from mimicry_discovery.structure.tcrdock_adapter import TCRDockPredictor

_PLACEHOLDER_CONFIDENCE = 0.5  # matches tests/integration's placeholder exactly

_BACKENDS: dict[StructureBackend, type[StructurePredictor]] = {
    StructureBackend.TCRDOCK: TCRDockPredictor,
    StructureBackend.ESMFOLD: ESMFoldPredictor,
    StructureBackend.AF_MULTIMER: AlphaFoldMultimerPredictor,
}


def _build_predictor(
    structure_config: StructurePredictionConfig,
    output_dir: str,
    hla_reference_path: str | None,
) -> StructurePredictor:
    """Construct the configured backend, wiring in the HLA reference if usable.

    Args:
        structure_config: Which backend to use and its settings.
        output_dir: Directory predicted structures are written to.
        hla_reference_path: Path to an IMGT/HLA protein FASTA file, or
            None/nonexistent if not built yet (roadmap Sprint 2). Only
            consulted for the ``esmfold`` backend, which needs it to
            build its complex sequence -- ``tcrdock`` and
            ``af_multimer`` don't use this lookup.

    Returns:
        A configured predictor instance for ``structure_config.backend``.
    """
    predictor_cls = _BACKENDS[structure_config.backend]
    if (
        structure_config.backend == StructureBackend.ESMFOLD
        and hla_reference_path
        and Path(hla_reference_path).exists()
    ):
        return ESMFoldPredictor.from_reference_fasta(
            output_dir=output_dir, hla_fasta_path=hla_reference_path
        )
    return predictor_cls(output_dir=output_dir)


def main(
    candidates_path: str,
    structure_config: StructurePredictionConfig,
    output_dir: str,
    results_out: str,
    smoke_test: bool = False,
    hla_reference_path: str | None = None,
) -> None:
    """Predict, or in smoke-test mode placeholder, structural confidence.

    Args:
        candidates_path: Path to the ``candidates.json`` written by the
            ``candidate_pairing`` rule.
        structure_config: Which backend to use and its settings.
        output_dir: Directory predicted structures are written to.
        results_out: Where to write per-candidate structural results
            (JSON).
        smoke_test: If True, skip the real (currently-scaffolded) backend
            entirely and use a fixed placeholder confidence. Never the
            default -- only for exercising the rest of the DAG before
            roadmap Sprints 2-3 land.
        hla_reference_path: Path to an IMGT/HLA protein FASTA (see
            ``scripts/download_reference_data.sh``). Only used by the
            ``esmfold`` backend; ignored (and safe to omit) otherwise.
            If not provided or the file doesn't exist yet, ``esmfold``
            falls back to an empty lookup and raises ``KeyError`` per
            candidate -- an honest failure, not a silent wrong answer.
    """
    candidates = json.loads(Path(candidates_path).read_text())
    results = []

    if smoke_test:
        for candidate in candidates:
            results.append(
                {
                    **candidate,
                    "structural_confidence": _PLACEHOLDER_CONFIDENCE,
                    "backend": "smoke_test_placeholder",
                }
            )
    else:
        predictor = _build_predictor(structure_config, output_dir, hla_reference_path)
        for candidate in candidates:
            request = StructurePredictionRequest(
                tcr_cdr3_beta=candidate["tcr_cdr3_beta"],
                tcr_v_gene=candidate["tcr_v_gene"],
                tcr_j_gene=candidate["tcr_j_gene"],
                peptide=candidate["self_peptide"],
                hla_allele=candidate["hla_allele"],
                tcr_cdr3_alpha=candidate.get("tcr_cdr3_alpha"),
                tcr_v_gene_alpha=candidate.get("tcr_v_gene_alpha"),
                tcr_j_gene_alpha=candidate.get("tcr_j_gene_alpha"),
            )
            # NOTE: raises NotImplementedError (or, for tcrdock with
            # unpaired input, ValueError; for esmfold with no reference
            # loaded, KeyError) until roadmap Sprint 2/3 finishes wiring
            # in a real backend -- expected behavior today, not a bug in
            # this script. See the module docstring above.
            result = predictor.predict(request)
            confidence = (result.interface_plddt or result.mean_plddt or 0.0) / 100.0
            results.append(
                {
                    **candidate,
                    "structural_confidence": confidence,
                    "backend": structure_config.backend.value,
                }
            )

    Path(results_out).parent.mkdir(parents=True, exist_ok=True)
    Path(results_out).write_text(json.dumps(results, indent=2))


if "snakemake" in globals():
    main(
        candidates_path=snakemake.input.candidates,  # noqa: F821
        structure_config=snakemake.params.structure_config,  # noqa: F821
        output_dir=snakemake.params.output_dir,  # noqa: F821
        results_out=snakemake.output.results,  # noqa: F821
        smoke_test=snakemake.params.smoke_test,  # noqa: F821
        hla_reference_path=snakemake.params.hla_reference_path,  # noqa: F821
    )
