"""ESMFold backend for TCR-pMHC structure prediction.

ESMFold predicts structure from a single sequence without an MSA, which
makes it attractive for the throughput this pipeline needs when
scanning many candidate TCR-peptide-HLA triples. We model the complex
as a single chain-break-separated sequence (TCR beta chain + peptide +
HLA heavy chain), following the convention AlphaFold-Multimer-style
complex folding uses.
"""

from __future__ import annotations

from pathlib import Path

from mimicry_discovery.io.hla_sequences import (
    load_hla_protein_sequences,
    to_two_field_lookup,
)
from mimicry_discovery.structure.base import (
    StructurePredictionRequest,
    StructurePredictionResult,
    StructurePredictor,
)


class ESMFoldPredictor(StructurePredictor):
    """Structure predictor backed by Meta AI's ESMFold.

    Attributes:
        output_dir: Directory predicted structures are written to.
        device: Torch device string, e.g. ``"cuda"`` or ``"cpu"``.
        hla_sequence_lookup: Two-field allele -> protein sequence map,
            e.g. from ``data/reference/hla_alleles/A_prot.fasta`` via
            :meth:`from_reference_fasta`. Empty by default.
    """

    def __init__(
        self,
        output_dir: Path | str,
        device: str = "cuda",
        hla_sequence_lookup: dict[str, str] | None = None,
    ) -> None:
        """Initialize the ESMFold predictor.

        Args:
            output_dir: Directory to write predicted PDB structures to.
            device: Torch device to run inference on.
            hla_sequence_lookup: Two-field allele -> protein sequence
                map. Defaults to empty; requests for an allele not in
                this map raise ``KeyError`` at predict time (see
                :meth:`_build_complex_sequence`).
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.device = device
        self.hla_sequence_lookup = hla_sequence_lookup or {}
        self._model: object | None = None  # lazily loaded, see _load_model

    @classmethod
    def from_reference_fasta(
        cls, output_dir: Path | str, hla_fasta_path: Path | str, device: str = "cuda"
    ) -> "ESMFoldPredictor":
        """Build a predictor with its HLA lookup loaded from a reference FASTA.

        Args:
            output_dir: Directory to write predicted PDB structures to.
            hla_fasta_path: Path to an IMGT/HLA protein FASTA file (see
                ``scripts/download_reference_data.sh``).
            device: Torch device to run inference on.

        Returns:
            A configured :class:`ESMFoldPredictor`.
        """
        full_res = load_hla_protein_sequences(hla_fasta_path)
        return cls(
            output_dir=output_dir,
            device=device,
            hla_sequence_lookup=to_two_field_lookup(full_res),
        )

    def _load_model(self) -> None:
        """Lazily load the ESMFold model onto ``self.device``.

        Raises:
            NotImplementedError: Always, until the TODO below is
                implemented.
        """
        # TODO(structure-prediction): load the real model, e.g.:
        #     from transformers import AutoTokenizer, EsmForProteinFolding
        #     self._tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
        #     self._model = EsmForProteinFolding.from_pretrained(
        #         "facebook/esmfold_v1"
        #     ).to(self.device)
        # Confirm the exact class/weights path against the `transformers`
        # version pinned in environment.yml before wiring this in.
        raise NotImplementedError("ESMFold model loading is a scaffold -- see TODO above.")

    def _build_complex_sequence(self, request: StructurePredictionRequest) -> str:
        """Build the single-chain, chain-break-joined complex sequence.

        Args:
            request: The TCR/peptide/HLA combination to model.

        Returns:
            A single amino-acid string representing the TCR beta chain,
            peptide, and HLA heavy chain joined with ``":"`` chain
            breaks, in the convention ESMFold's multimer mode expects.

        Raises:
            KeyError: If ``request.hla_allele``'s two-field form has no
                entry in ``self.hla_sequence_lookup`` -- populate it via
                :meth:`from_reference_fasta` or pass a lookup covering
                the alleles you need.
        """
        two_field = request.hla_allele.removeprefix("HLA-")
        if two_field not in self.hla_sequence_lookup:
            raise KeyError(
                f"No reference sequence for HLA allele '{request.hla_allele}' "
                f"(looked up as '{two_field}') -- build hla_sequence_lookup via "
                "ESMFoldPredictor.from_reference_fasta() first."
            )
        hla_sequence = self.hla_sequence_lookup[two_field]
        return f"{request.tcr_cdr3_beta}:{request.peptide}:{hla_sequence}"

    def predict(self, request: StructurePredictionRequest) -> StructurePredictionResult:
        """Predict a TCR-pMHC complex structure using ESMFold.

        Args:
            request: The TCR/peptide/HLA combination to model.

        Returns:
            A :class:`StructurePredictionResult` with the modeled
            complex's PDB path and mean pLDDT.

        Raises:
            KeyError: If the request's HLA allele has no reference
                sequence loaded (see :meth:`_build_complex_sequence`).
            NotImplementedError: Always, until model loading (see
                :meth:`_load_model`) is implemented.
        """
        sequence = self._build_complex_sequence(request)
        if self._model is None:
            self._load_model()
        # TODO(structure-prediction): run real inference, e.g.:
        #     with torch.no_grad():
        #         output = self._model.infer_pdb(sequence)
        #     pdb_path = self.output_dir / f"{request.hla_allele}_{request.peptide}.pdb"
        #     pdb_path.write_text(output)
        #     return StructurePredictionResult(
        #         pdb_path=str(pdb_path), backend="esmfold",
        #     )
        raise NotImplementedError(
            "ESMFoldPredictor.predict inference is a scaffold -- wire in "
            "real inference (see TODO above) before use."
        )
