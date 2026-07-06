"""Rule: generate TCR x tumor-neoantigen x self-peptide candidates.

See workflow/scripts/candidate_pairing.py's docstring for why this rule
runs (and produces a valid, if currently-empty, candidate list) even
before the self-antigen reference exists.
"""


rule candidate_pairing:
    input:
        clonotypes=rules.ingest_qc.output.clonotypes,
        neoantigens=rules.ingest_qc.output.neoantigens,
    output:
        candidates=str(pipeline_config.paths.output_dir / "{sample}" / "candidates.json"),
    params:
        reference_dir=str(pipeline_config.paths.reference_dir),
    log:
        str(pipeline_config.paths.output_dir / "{sample}" / "logs" / "candidate_pairing.log"),
    conda:
        "envs/tcr_qc.yaml"
    container:
        "docker://ghcr.io/YOUR-ORG/mimicry-discovery:scoring-latest"
    script:
        "../scripts/candidate_pairing.py"
