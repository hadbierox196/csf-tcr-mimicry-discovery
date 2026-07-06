"""Rule: export the top-ranked candidates for tetramer ordering."""


rule wetlab_export:
    input:
        risk_scores=rules.risk_scoring.output.risk_scores,
    output:
        export=str(pipeline_config.paths.output_dir / "{sample}" / "tetramer_order_export.csv"),
    params:
        top_n=config.get("wetlab_export_top_n", 10),
    log:
        str(pipeline_config.paths.output_dir / "{sample}" / "logs" / "wetlab_export.log"),
    conda:
        "envs/scoring.yaml"
    container:
        "docker://ghcr.io/YOUR-ORG/mimicry-discovery:scoring-latest"
    script:
        "../scripts/wetlab_export.py"
