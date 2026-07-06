"""TCRdock backend for TCR-pMHC structure prediction.

TCRdock (Bradley lab, https://github.com/phbradley/TCRdock) is a
CLI/subprocess-driven tool, not an importable Python API: it works by
writing a TSV of TCR/peptide/MHC targets, running ``setup_for_alphafold.py``
to build AlphaFold inputs from them, then ``run_prediction.py`` to
actually fold them with a TCR-specialized AlphaFold pipeline (three
``model_2_ptm`` simulations per target, keeping the one with the lowest
TCR:pMHC interface PAE). This adapter shells out to both scripts in
sequence, matching that real interface, rather than the earlier
placeholder that guessed at a Python function call.

Verified against TCRdock's current README (fetched 2026-07) rather than
assumed from memory -- see ``_to_tcrdock_mhc_string``'s docstring for the
one specific detail (exact MHC allele string formatting) that couldn't
be confirmed from the README alone and still needs checking against a
real TCRdock install before first use. Output parsing (the final step,
turning TCRdock's PDB + PAE output into a
:class:`~mimicry_discovery.structure.base.StructurePredictionResult`) is
also not yet wired in for the same reason -- see the TODO in
:meth:`TCRDockPredictor.predict`.
"""

from __future__ import annotations

import csv
import subprocess
from pathlib import Path

from mimicry_discovery.structure.base import (
    StructurePredictionRequest,
    StructurePredictionResult,
    StructurePredictor,
)

_TARGETS_TSV_COLUMNS = [
    "organism", "mhc_class", "mhc", "peptide",
    "va", "ja", "cdr3a", "vb", "jb", "cdr3b",
]


def _to_tcrdock_mhc_string(hla_allele: str) -> str:
    """Convert a two-field IPD-IMGT/HLA string to TCRdock's ``mhc`` column format.

    TCRdock's README documents the 10-column targets TSV schema but its
    worked example file (``examples/benchmark/single_target.tsv``) wasn't
    fetchable from this environment, so the exact string format TCRdock
    expects here isn't independently confirmed. This strips the ``"HLA-"``
    prefix as the most likely convention (IMGT/HLA nomenclature without
    the database prefix, e.g. ``"A*02:01"``), matching how most non-IPD
    tools in this space format class I alleles -- but treat this function
    as the thing to verify first if TCRdock rejects a target.

    Args:
        hla_allele: Two-field IPD-IMGT/HLA string, e.g. ``"HLA-A*02:01"``.

    Returns:
        The allele string with the ``"HLA-"`` prefix removed.
    """
    return hla_allele.removeprefix("HLA-")


class TCRDockPredictor(StructurePredictor):
    """Structure predictor backed by TCRdock.

    Requires a TCRdock repo checkout (for ``setup_for_alphafold.py`` and
    ``run_prediction.py``) and a working AlphaFold installation with
    parameters available locally -- TCRdock's own README covers that
    AlphaFold-side setup in detail; this class doesn't attempt to
    replicate it, only to drive the two scripts once they're available.

    Attributes:
        output_dir: Directory targets TSVs, setup files, and predicted
            structures are written to.
        tcrdock_dir: Path to a TCRdock repo checkout.
        alphafold_data_dir: Path to the AlphaFold ``params/`` parent
            directory.
        model_name: AlphaFold parameter set. TCRdock's own benchmark
            used ``"model_2_ptm"``.
    """

    def __init__(
        self,
        output_dir: Path | str,
        tcrdock_dir: Path | str = "/opt/TCRdock",
        alphafold_data_dir: Path | str = "/data/alphafold",
        model_name: str = "model_2_ptm",
    ) -> None:
        """Initialize the TCRdock predictor.

        Args:
            output_dir: Directory to write targets TSVs, setup files, and
                predicted structures to. Created if it doesn't exist.
            tcrdock_dir: Path to a TCRdock repo checkout.
            alphafold_data_dir: Path to the AlphaFold ``params/`` parent
                directory (passed as AlphaFold's ``--data_dir``).
            model_name: AlphaFold parameter set to request.
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.tcrdock_dir = Path(tcrdock_dir)
        self.alphafold_data_dir = Path(alphafold_data_dir)
        self.model_name = model_name

    def _write_targets_tsv(
        self, request: StructurePredictionRequest, targets_path: Path
    ) -> None:
        """Write TCRdock's required 10-column targets TSV for one request.

        Args:
            request: The TCR/peptide/HLA combination to model. Must be
                paired (both chains present) -- TCRdock's schema has no
                way to represent a single-chain target.
            targets_path: Where to write the one-row TSV.

        Raises:
            ValueError: If ``request.is_paired`` is False.
        """
        if not request.is_paired:
            raise ValueError(
                "TCRDockPredictor requires a paired (alpha+beta) TCR -- "
                f"got beta-only input for CDR3-beta '{request.tcr_cdr3_beta}'. "
                "Route beta-only clonotypes to ESMFoldPredictor instead, "
                "or resolve the alpha chain before calling TCRdock."
            )
        with targets_path.open("w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(_TARGETS_TSV_COLUMNS)
            writer.writerow([
                request.organism,
                request.mhc_class,
                _to_tcrdock_mhc_string(request.hla_allele),
                request.peptide,
                request.tcr_v_gene_alpha,
                request.tcr_j_gene_alpha,
                request.tcr_cdr3_alpha,
                request.tcr_v_gene,
                request.tcr_j_gene,
                request.tcr_cdr3_beta,
            ])

    def predict(self, request: StructurePredictionRequest) -> StructurePredictionResult:
        """Predict a TCR-pMHC complex structure using TCRdock.

        Runs the real two-step TCRdock pipeline: writes a one-row
        targets TSV, calls ``setup_for_alphafold.py`` to build AlphaFold
        inputs, then ``run_prediction.py`` to fold them.

        Args:
            request: The TCR/peptide/HLA combination to model. Must be
                paired -- see :meth:`_write_targets_tsv`.

        Returns:
            A :class:`StructurePredictionResult` with the modeled
            complex's PDB path and interface confidence.

        Raises:
            ValueError: If ``request`` isn't paired.
            subprocess.CalledProcessError: If either TCRdock script
                exits non-zero (e.g. TCRdock/AlphaFold isn't actually
                installed at the configured paths).
            NotImplementedError: Always, after both subprocess calls
                succeed -- output parsing isn't wired in yet (see the
                module docstring).
        """
        targets_path = self.output_dir / "targets.tsv"
        self._write_targets_tsv(request, targets_path)

        setup_dir = self.output_dir / "setup"
        subprocess.run(
            [
                "python", str(self.tcrdock_dir / "setup_for_alphafold.py"),
                "--targets_tsvfile", str(targets_path),
                "--output_dir", str(setup_dir),
            ],
            check=True,
            cwd=self.tcrdock_dir,
        )

        outfile_prefix = str(self.output_dir / "prediction")
        subprocess.run(
            [
                "python", str(self.tcrdock_dir / "run_prediction.py"),
                "--targets", str(setup_dir / "targets.tsv"),
                "--outfile_prefix", outfile_prefix,
                "--model_names", self.model_name,
                "--data_dir", str(self.alphafold_data_dir),
            ],
            check=True,
            cwd=self.tcrdock_dir,
        )

        # TODO(structure-prediction): parse the real output here. TCRdock's
        # repo includes add_pmhc_tcr_pae_to_tsvfile.py, which strongly
        # suggests PAE/confidence lands in a TSV alongside the PDB, but
        # this environment couldn't fetch that script to confirm its
        # exact output column names -- run the pipeline once against a
        # real install, inspect `{outfile_prefix}*`, and finish this.
        raise NotImplementedError(
            "TCRDockPredictor.predict ran the real TCRdock subprocess "
            "pipeline above, but output parsing is not yet wired in -- "
            "see the TODO above."
        )
