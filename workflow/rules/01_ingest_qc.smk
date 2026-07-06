"""Rule: parse + QC each sample's CSF TCR-seq and tumor neoantigen calls.

Reads raw paths from the ``samples`` DataFrame (parsed once at the top of
the Snakefile from ``config/samples.tsv``) and writes standardized JSON
per sample. This is the only rule that touches raw input files directly --
every downstream rule reads the JSON this one produces.
"""


def _sample_path(wildcards, column: str) -> str:
    """Look up a per-sample raw file path from the samples manifest.

    Args:
        wildcards: Snakemake wildcards object; only ``.sample`` is used.
        column: Which samples.tsv column to read.

    Returns:
        The raw file path for that sample and column.
    """
    return samples.set_index("sample_id").loc[wildcards.sample, column]


rule ingest_qc:
    input:
        tcr_vdj=lambda wc: _sample_path(wc, "tcr_vdj_path"),
        neoantigens=lambda wc: _sample_path(wc, "neoantigen_report_path"),
    output:
        clonotypes=str(pipeline_config.paths.output_dir / "{sample}" / "clonotypes.json"),
        neoantigens=str(pipeline_config.paths.output_dir / "{sample}" / "neoantigens.json"),
    params:
        sample_id=lambda wc: wc.sample,
        qc=pipeline_config.qc,
    log:
        str(pipeline_config.paths.output_dir / "{sample}" / "logs" / "ingest_qc.log"),
    conda:
        "envs/tcr_qc.yaml"
    container:
        "docker://ghcr.io/YOUR-ORG/mimicry-discovery:scoring-latest"
    script:
        "../scripts/ingest_qc.py"
