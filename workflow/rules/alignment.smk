rule alignment:
    input: 
        fastq_merged = "results/merged/{sample}.fastq.gz"
    output: 
        bam = "results/alignment/{sample}_Aligned.sortedByCoord.out.bam"
    conda:
        config["conda_envs"]["aligners"]
    threads: get_resource("alignment", "threads")
    resources:
        mem_mb = get_resource(config, "alignment", "mem_mb"),
        runtime = get_resource(config, "alignment", "runtime")
    params:
        genome_index = config["genome"]["index"],
        outdir = lambda wildcards: os.path.join("results", "alignment", f"{wildcards.sample}_")
        extra = config["parameters"]["alignment"]["extra"]
    log:
        "log/alignment/{sample}.log"
    benchmark:
        "benchmarks/alignment/{sample}.bmk"
    shell: """
        STAR --runThreadN {threads} \
            --genomeDir {params.genome_index} 
            --readFilesIn {input.fastq_merged} \
            --readFilesCommand gunzip -c \
            --outFileNamePrefix {params.outdir} \
            --outSAMtype BAM SortedByCoordinate
            {params.extra} 2> {log}
    """

rule bam_indexing:
    input:
        bam = "results/alignment/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        bam_bai = "results/alignment/{sample}_Aligned.sortedByCoord.out.bam.bai"
    conda:
        config["conda_envs"]["aligners"]
    threads: get_resource(config, "bam_indexing", "threads")
    resources:
        mem_mb = get_resource(config, "bam_indexing", "mem_mb"),
        runtime = get_resource(config, "bam_indexing", "runtime"),
        extra = config["parameters"]["bam_indexing"]["extra"]
    log:
        "log/bam_indexing/{sample}.log"
    benchmark:
        "benchmarks/bam_indexing/{sample}.bmk"
    shell:"""
        samtools index \
            -@ {threads} \
            {input.bam} \
            {params.extra} 2> {log}
    """