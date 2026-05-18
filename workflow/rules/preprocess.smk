if UMIs:

    rule umi_extract:
        input:
            "data/{sample}_{seq_lane}.fastq.gz"
        output:
            "results/umi_extract/{sample}_{seq_lane}.fastq.gz"
        conda:
            config["conda_envs"]["umi_tools"]
        threads: 
            get_resource(config, "umi_extract", "threads")
        resources:
            mem_mb = get_resource(config, "umi_extract", "mem_mb"),
            runtime = get_resource(config, "umi_extract", "walltime")
        params:
            extract_method = config["umi_extract"]["extract_method"],
            extra = config["umi_extract"]["extra"]
        log:
            "log/umi_extract/{sample}_{seq_lane}.log"
        benchmark:
            "benchmarks/umi_extract/{sample}_{seq_lane}.bmk"
        shell: """
            umi_tools extract \ 
                --stdin={input} \
                --extract-method {params.extract_method} \
                {params.extra} \
                --log={log} \ 
                --stdout={output}
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
        threads: 
            get_resource(config, "umi_dedup", "threads")
        resources:
            mem_mb = get_resource(config, "umi_dedup", "mem_mb"),
            runtime = get_resource(config, "umi_dedup", "runtime")
        params:
            method = config["umi_dedup"]["method"],
            stats_dir = lambda wildcards: f"results/alignment/dedup/{wildcards.sample}"
            extra = config["umi_dedup"]["extra"]
        log:
            "log/umi_dedup/{sample}.log"
        benchmark:
            "benchmarks/umi_dedup/{sample}.bmk"
        shell:"""
            umi_tools dedup \
                -I {input.bam} \
                --output-stats={params.stats_dir} \
                --method={params.method} \
                {params.extra} \
                --log={log} \
                -S {output.dedup} 
        """