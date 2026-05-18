rule merged_fastq:
    input: 
        lambda wildcards: expand(
            "results/trimmed/{sample}/{sample}_{seq_lane}_trimmed.fastq.gz",
            sample = wildcards.sample,
            seq_lane = config["seq_lane"]
        )
    output: 
        fastq_merged = "results/merged/{sample}.fastq.gz"
    conda: 
        config["conda_envs"]["qc"]
    threads: get_resource(config, "default", "threads")
    resources:
        mem_mb = get_resource(config, "default", "mem_mb"),
        runtime = get_resource(config, "default", "runtime")
    log:
        "log/merged/{sample}_merging.log"
    shell:
        "cat {input} > {output.fastq_merged} 2> {log}"
