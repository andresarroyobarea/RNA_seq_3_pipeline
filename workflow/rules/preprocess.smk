if UMIs:

    rule umi_extract:
        input:
            "data/{sample}_{seq_lane}.fastq.gz"
        output:
            "results/umi_extract/{sample}_{seq_lane}.fastq.gz"
        conda:
            config["conda_envs"]["umi_tools"]
        threads: get_resource(config, "umi_extract", "threads")
        resources:
            mem_mb = get_resource(config, "umi_extract", "mem_mb"),
            runtime = get_resource(config, "umi_extract", "walltime")
        params:
            pattern = lambda wildcards: config["umi_processing"]["pattern"]
        log:
            "log/umitools/extract/{sample}_{seq_lane}.log"
        benchmark:
            "benchmarks/{sample}_{seq_lane}_umi_tools_extract.bmk"
        shell: """
            umi_tools extract --stdin={input} \
                --extract-method regex --bc-pattern="{params.pattern}" \
                --log={log} --stdout={output}
        """

    rule umi_dedup:
        input:
            bam = "results/alignment/{sample}_Aligned.sortedByCoord.out.bam",
            bam_bai = "results/alignment/{sample}_Aligned.sortedByCoord.out.bam.bai"
        output:
            dedup = "results/alignment/dedup/{sample}_Aligned.sortedByCoord.out.bam",
            stat = "results/alignment/dedup/{sample}_per_umi.tsv"
        conda:
            config["conda_envs"]["umi_tools"]
        threads: get_resource(config, "umi_dedup", "threads")
        resources:
            mem_mb = get_resource(config, "umi_dedup", "mem_mb"),
            runtime = get_resource(config, "umi_dedup", "runtime")
        params:
            stats_dir = lambda wildcards: f"results/alignment/dedup/{wildcards.sample}"
        log:
            "log/umitools/dedup/{sample}.log"
        benchmark:
            "benchmarks/{sample}_umi_dedup.bmk"
        shell:"""
            umi_tools dedup -I {input.bam} --log={log} -S {output.dedup} --output-stats={params.stats_dir} \
                --method=unique --multimapping-detection-method=NH
        """