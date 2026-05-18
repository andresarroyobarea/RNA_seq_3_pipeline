
# TODO: Feature counts per sample.
rule feature_counts:
    input:
        bam = expand(aligned_reads, sample = config["sample"])
    output: 
        feature_table = "results/feature_counts/counts.tsv" 
    conda: 
        config["conda_envs"]["quantification"]
    threads: 
        get_resource(config, "feature_counts", "threads")
    resources:
        mem_mb = get_resource(config, "feature_counts", "mem_mb"),
        runtime = get_resource(config, "feature_counts", "runtime")
    params:
        annotations = config["annotation"],
        extra = config["parameters"]["feature_counts"]["extra"]
    log:
        "log/featureCounts/featureCounts.log"
    benchmark:
        "benchmarks/featureCounts/subread_featureCounts.bmk"
    shell:"""
        featureCounts -a {params.annotations} \ 
            -T {threads} \
            -o {output.feature_table} \
            {input.bam} \
            {params.extra} > {log}
    """