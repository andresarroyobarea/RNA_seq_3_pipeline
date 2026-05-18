
rule feature_counts:
    input:
        bam = expand(aligned_reads, sample = config["sample"])
    output: 
        feature_table = "results/feature_counts/counts.tsv" 
    conda: 
        config["conda_envs"]["quantification"]
    threads: get_resource(config, "feature_counts", "threads")
    resources:
        mem_mb = get_resource(config, "feature_counts", "mem_mb"),
        runtime = get_resource(config, "feature_counts", "runtime")
    params:
        annotations = config["annotation"]
    log:
        "log/featureCounts/featureCounts.log"
    benchmark:
        "benchmarks/subread_featureCounts.bmk"
    shell:"""
        featureCounts -a {params.annotations} -O -F GTF -t gene -g gene_id \
            --extraAttributes gene_name,transcript_name -s 1 -T {threads} \
            -o {output.feature_table} {input.bam} 2> {log}
    """