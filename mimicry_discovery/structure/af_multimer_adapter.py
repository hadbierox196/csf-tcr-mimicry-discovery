"""AlphaFold-Multimer backend for TCR-pMHC structure prediction.

AlphaFold-Multimer is the highest-confidence, highest-cost backend of
the three supported here -- typically reserved for a shortlist of
top-ranked candidates from a faster backend (TCRdock/ESMFold) given its
runtime and the genetic-database search it depends on.
"""

from __future__ import annotations

from pathlib import Path

from mimicry_discovery.structure.base import (
    StructurePredictionRequest,
    StructurePredictionResult,
    StructurePredictor,
)


class AlphaFoldMultimerPredictor(StructurePredictor):
    """Structure predictor backed by AlphaFold-Multimer.

    Attributes:
        output_dir: Directory predicted structures, PAE, and confidence
            JSON files are written to.
        model_preset: Which AlphaFold-Multimer model preset to run.
    """

    def __init__(self, output_dir: Path | str, model_preset: str = "multimer") -> None:
        """Initialize the AlphaFold-Multimer predictor.

        Args:
            output_dir: Directory to write predicted structures and
                confidence outputs to.
            model_preset: AlphaFold model preset name.
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.model_preset = model_preset

    def predict(self, request: StructurePredictionRequest) -> StructurePredictionResult:
        """Predict a TCR-pMHC complex structure using AlphaFold-Multimer.

        Args:
            request: The TCR/peptide/HLA combination to model.

        Returns:
            A :class:`StructurePredictionResult` with the modeled
            complex's PDB path, mean pLDDT, and interface pLDDT.

        Raises:
            RuntimeError: If AlphaFold-Multimer fails to produce a
                structure (e.g. an MSA/template search failure).
            NotImplementedError: Always, until the TODO below is
                implemented.
        """
        # TODO(structure-prediction): shell out to `run_alphafold.py`
        # (or the containerized equivalent in
        # containers/Dockerfile.structure_prediction) with a FASTA built
        # from `request`, using --model_preset=multimer. Parse the
        # resulting ranking_debug.json + relaxed PDB into the result
        # below. This is the most expensive backend -- gate its use to a
        # shortlist from a faster backend rather than every candidate.
        raise NotImplementedError(
            "AlphaFoldMultimerPredictor.predict is a scaffold -- wire in "
            "the real subprocess/API call (see TODO above) before use."
        )
