"""Rule: predict TCR:self-peptide-HLA complex structures per candidate.

See workflow/scripts/structure_prediction.py's docstring -- this rule
raises NotImplementedError today unless invoked with
`--config smoke_test=true`, which is correct (the real backends are
still scaffolded, roadmap Sprints 2-3).
"""


rule structure_prediction:
    input:
        candidates=rules.candidate_pairing.output.candidates,
    output:
        results=str(pipeline_config.paths.output_dir / "{sample}" / "structure_results.json"),
    params:
        structure_config=pipeline_config.structure_prediction,
        output_dir=str(pipeline_config.paths.output_dir / "{sample}" / "structures"),
        smoke_test=config.get("smoke_test", False),
        hla_reference_path=str(pipeline_config.paths.reference_dir / "hla_alleles" / "A_prot.fasta"),
    log:
        str(pipeline_config.paths.output_dir / "{sample}" / "logs" / "structure_prediction.log"),
    conda:
        "envs/structure_prediction.yaml"
    container:
        "docker://ghcr.io/YOUR-ORG/mimicry-discovery:structure-prediction-latest"
    resources:
        gpu=1 if not config.get("smoke_test", False) else 0,
        mem_mb=32000,
    script:
        "../scripts/structure_prediction.py"
