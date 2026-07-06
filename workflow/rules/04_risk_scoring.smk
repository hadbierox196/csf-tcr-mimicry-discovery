"""Rule: compute the composite mimicry risk score per candidate.

The only rule in the DAG with no NotImplementedError caveat -- scoring is
fully implemented (see mimicry_discovery/scoring/).
"""


rule risk_scoring:
    input:
        results=rules.structure_prediction.output.results,
    output:
        risk_scores=str(pipeline_config.paths.output_dir / "{sample}" / "risk_scores.json"),
    params:
        scoring_config=pipeline_config.scoring,
    log:
        str(pipeline_config.paths.output_dir / "{sample}" / "logs" / "risk_scoring.log"),
    conda:
        "envs/scoring.yaml"
    container:
        "docker://ghcr.io/YOUR-ORG/mimicry-discovery:scoring-latest"
    script:
        "../scripts/risk_scoring.py"
