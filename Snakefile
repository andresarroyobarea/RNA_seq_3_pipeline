configfile: "config/config.yaml"

import os

UMIs = config["umi_processing"]["enabled"]

## Let stablish where the trimming input will be saved.
if UMIs:
    dir_in_trim = "results/umi_extract"
else:
    dir_in_trim = "data"

## Let stablish the path for BAM files in each case
if UMIs:
    aligned_reads = "results/alignment/dedup/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
else:
    aligned_reads = "results/alignment/{sample}/{sample}_Aligned.sortedByCoord.out.bam"

rule all:
    input:
        "results/feature_counts/counts.tsv",
        "results/QC/alignment/qualimap/multi_bamqc/multisampleBamQcReport.html",
        "results/MultiQC/multiqc_report.html"

rule fastqc_raw:
    input: 
        fastq = "data/{sample}_{seq_lane}.fastq.gz"
    output:
        html = "results/QC/raw/qc_per_lane/{sample}/{sample}_{seq_lane}_fastqc.html",
        zip = "results/QC/raw/qc_per_lane/{sample}/{sample}_{seq_lane}_fastqc.zip"
    conda:
        config["conda_envs"]["qc"]
    threads: 2
    params: 
        outdir = lambda wildcards, output: os.path.dirname(output.html)
    log:
        fastqc = "log/QC/raw/qc_per_lane/{sample}_{seq_lane}_fastqc.log"
    shell:
        "mkdir -p {params.outdir} &&"
        "fastqc --outdir {params.outdir} --threads {threads} {input.fastq} 2> {log.fastqc} "

rule bbduk_se:
    input:
        sample = [dir_in_trim + "/{sample}_{seq_lane}.fastq.gz"]
    output:
        trimmed = "results/trimmed/{sample}/{sample}_{seq_lane}_trimmed.fastq.gz",
        singleton = "results/trimmed/{sample}/{sample}_{seq_lane}_single.fastq.gz",
        discarded = "results/trimmed/{sample}/{sample}_{seq_lane}_discarded.fastq.gz",
        stats = "results/trimmed/{sample}/{sample}_{seq_lane}_stats.txt"
    conda: 
        config["conda_envs"]["preprocessing"]
    threads: 2
    params:
        adapters = config["adapters"],
        polyA = config["polyA"]
    log:
        "log/bbduk/{sample}_{seq_lane}_bbduk.log"
    benchmark:
        "benchmarks/{sample}_{seq_lane}_bbduk.bmk"
    shell:
        """
        bbduk.sh in={input.sample} \
            out={output.trimmed} \
            outs={output.singleton} \
            outm={output.discarded} \
            stats={output.stats} \
            ref={params.adapters},{params.polyA} \
            k=13 ktrim=r mink=5 qtrim=r trimq=20 useshortkmers=t minlength=20 > {log} 2>&1
        """

rule fastqc_trim:
    input: 
        fastq = "results/trimmed/{sample}/{sample}_{seq_lane}_trimmed.fastq.gz"
    output:
        html = "results/QC/trimmed/qc_per_lane/{sample}/{sample}_{seq_lane}_trimmed_fastqc.html",
        zip = "results/QC/trimmed/qc_per_lane/{sample}/{sample}_{seq_lane}_trimmed_fastqc.zip",
    conda:
        config["conda_envs"]["qc"]
    threads: 2
    params: 
        outdir = lambda wildcards, output: os.path.dirname(output.html)
    log:
        fastqc = "log/QC/trimmed/qc_per_lane/{sample}_{seq_lane}_trimmed_fastqc.log"
    shell:
        "mkdir -p {params.outdir} &&"
        "fastqc --outdir {params.outdir} --threads {threads} {input.fastq} 2> {log.fastqc} "

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
    log:
        "log/merged/{sample}_merging.log"
    shell:
        "cat {input} > {output.fastq_merged} 2> {log}"

rule fastq_screen:
    input: 
        fastq_merged = expand("results/merged/{sample}.fastq.gz", sample = config["sample"])
    output: 
        fastq_screen_txt = "results/QC/fastq_screen/{sample}_fastq_screen.txt",
        fastq_screen_png = "results/QC/fastq_screen/{sample}_fastq_screen.png"
    conda:
        config["conda_envs"]["fastq_screen"]
    threads: 1
    resources:
        mem_mb=28728
    params:
        fastq_screen_config = config["fastq_screen_conf"],
        aligner = config["fastq_screen_aling"],
        outdir = "results/QC/fastq_screen/"
    log:
        log = "log/QC/fastq_screen/{sample}_fastq_screen.log"
    benchmark:
        "benchmarks/{sample}_fastq_screen.bmk"
    shell:"""
        fastq_screen {input.fastq_merged} --aligner {params.aligner} \
            --conf {params.fastq_screen_config} --outdir {params.outdir} \
            -threads {threads} 2> {log}
    """

rule alignment:
    input: 
        fastq_merged = "results/merged/{sample}.fastq.gz"
    output: 
        bam = "results/alignment/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    conda:
        config["conda_envs"]["aligners"]
    threads: 5
    resources:
        mem_mb=18432
    params:
        genome_index = config["genome_index"],
        outdir = lambda wildcards, output: os.path.dirname(output.bam)
    log:
        "log/alignment/{sample}_alignment.log"
    benchmark:
        "benchmarks/{sample}_alignment.bmk"
    shell: """
        STAR --runThreadN {threads} --genomeDir {params.genome_index} --genomeLoad LoadAndKeep --readFilesIn {input.fastq_merged} \
            --readFilesCommand gunzip -c --outFilterType BySJout --outFilterMultimapNmax 25 --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.3 \
            --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --limitBAMsortRAM 12000000000 \
            --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.outdir} 2> {log}
    """

rule bam_indexing:
    input:
        bam = "results/alignment/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        bam_bai = "results/alignment/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
    conda:
        config["conda_envs"]["aligners"]
    threads: 3
    log:
        "log/bam_indexing/{sample}.log"
    shell:
        "samtools index -@ {threads} {input.bam} "

rule fastqc_alignment:
    input:
        bam = aligned_reads
    output:
        html = "results/QC/alignment/FastQC/{sample}_Aligned.sortedByCoord.out_fastqc.html",
        zip = "results/QC/alignment/FastQC/{sample}_Aligned.sortedByCoord.out_fastqc.zip",
    conda:
        config["conda_envs"]["qc"]
    threads: 5
    params: 
        outdir = lambda wildcards, output : os.path.dirname(output.html)
    log:
        "log/QC/alignment/FastQC/{sample}_alignment_fastqc.log"
    shell:
        "mkdir -p {params.outdir} &&"
        "fastqc --outdir {params.outdir} --threads {threads} {input.bam} 2> {log} "

rule qualimap_bamqc:
    input:
        bam = aligned_reads
    output:
        qmap_report = "results/QC/alignment/qualimap/bamqc/{sample}/qualimapReport.html",
        genome_res = "results/QC/alignment/qualimap/bamqc/{sample}/genome_results.txt"
    conda:
        config["conda_envs"]["qualimap"]
    threads: 3
    resources:
        mem_mb=26624
    params:
        qmap_genome = config["qualimap"]["genome"],
        annotation = config["annotation"],
        outdir = "results/QC/alignment/qualimap/bamqc/{sample}",
        mem = config["qualimap"]["mem"]
    log: 
        "log/QC/alignment/qualimap/bamqc/{sample}_qualimap_bamqc.log"
    benchmark:
        "benchmarks/{sample}_qualimap_bamqc.bmk"
    shell: """
        qualimap bamqc -bam {input.bam} -gd {params.qmap_genome} -gff {params.annotation} \
            -hm 3 -nr 1000 -nt {threads} --outdir {params.outdir} -p strand-specific-forward \
            --java-mem-size={params.mem} 2> {log}
    """

rule qualimap_multi_bamqc:
    input:
        expand("results/QC/alignment/qualimap/bamqc/{sample}/qualimapReport.html", sample = config["sample"])
    output:
        "results/QC/alignment/qualimap/multi_bamqc/multisampleBamQcReport.html"
    conda:
        config["conda_envs"]["qualimap"]
    threads: 3
    resources:
        mem_mb=26624
    params: 
        qmap_input = "metadata/qualimap_multi_bamqc_input.txt",
        outdir = "results/QC/alignment/qualimap/multi_bamqc"
    log: 
        "log/QC/alignment/qualimap/mutli_bamqc/qualimap_multi_bamqc.log"
    benchmark:
        "benchmarks/qualimap_multi_bamqc.bmk"
    shell:
        "qualimap multi-bamqc -d {params.qmap_input} --outdir {params.outdir} 2> {log} "


rule qualimap_rnaseq:
    input:  
        bam = aligned_reads
    output: 
        "results/QC/alignment/qualimap/rnaseq/{sample}/qualimapReport.html",
        "results/QC/alignment/qualimap/rnaseq/{sample}/rnaseq_qc_results.txt"
    conda:
        config["conda_envs"]["qualimap"]
    params:
        annotation = config["annotation"],
        outdir = "results/QC/alignment/qualimap/rnaseq/{sample}",
        mem = config["qualimap"]["mem"]
    log:
        "log/QC/alignment/qualimap/rnaseq/{sample}_qualiamp_rnaseq.log"
    benchmark:
        "benchmarks/{sample}_qualimap_rnaseq.bmk"
    shell:"""
        qualimap rnaseq -bam {input.bam} -gtf {params.annotation} \
            -outdir {params.outdir} -p strand-specific-forward \
            --java-mem-size={params.mem} 2> {log}
    """

rule rseqc_strand:
    input: 
        bam = aligned_reads
    output: 
        rseqc_out = "results/QC/alignment/rseqc/{sample}_strandiness.txt"
    conda: 
        config["conda_envs"]["rseqc"]
    params:
        bed_file = config["reference_bed"]
    log:
        "log/QC/alignment/rseqc/{sample}_rseqc.log"
    shell:
        "infer_experiment.py -i {input.bam} -r {params.bed_file} > {output.rseqc_out} 2> {log} "

rule samtools_stats_flagstat:
    input:
        bam = aligned_reads
    output:
        samtools_stats = "results/QC/alignment/samtools_stats/{sample}_Aligned.sortedByCoord.out.bam.stats",
        samtools_flagstat = "results/QC/alignment/samtools_stats/{sample}_Aligned.sortedByCoord.out.bam.flagstat"
    conda:
        config["conda_envs"]["aligners"]
    log:
        log_stats = "log/QC/alignment/samtools/{sample}_samtools_stats.log",
        log_flagstat = "log/QC/alignment/samtools/{sample}_samtools_flagstat.log"
    shell:"""
        samtools stats {input.bam} > {output.samtools_stats} 2> {log.log_stats} &&
        samtools flagstat {input.bam} > {output.samtools_flagstat} 2> {log.log_flagstat}
    """

rule feature_counts:
    input:
        bam = expand(aligned_reads, sample = config["sample"])
    output: 
        feature_table = "results/feature_counts/counts.tsv" 
    conda: 
        config["conda_envs"]["quantification"]
    params:
        annotations = config["annotation"]
    log:
        "log/featureCounts/featureCounts.log"
    benchmark:
        "benchmarks/subread_featureCounts.bmk"
    shell:"""
        featureCounts -a {params.annotations} -O -F GTF -t gene -g gene_id \
            --extraAttributes gene_name,transcript_name -s 1 -T 15 \
            -o {output.feature_table} {input.bam} 2> {log}
    """

rule multiqc:
    input: 
        seqs_QC_raw = expand("results/QC/raw/qc_per_lane/{sample}/{sample}_{seq_lane}_fastqc.html", sample = config["sample"], seq_lane = config["seq_lane"]),
        seqs_QC_trim = expand("results/QC/trimmed/qc_per_lane/{sample}/{sample}_{seq_lane}_trimmed_fastqc.html", sample = config["sample"], seq_lane = config["seq_lane"]),
        fastq_screen_txt = expand("results/QC/fastq_screen/{sample}/{sample}_fastq_screen.txt", sample = config["sample"]),
        fastq_screen_png = expand("results/QC/fastq_screen/{sample}/{sample}_fastq_screen.png", sample = config["sample"]),
        bam_QC = expand("results/QC/alignment/FastQC/{sample}_Aligned.sortedByCoord.out_fastqc.html", sample = config["sample"]),
        qmap_bamqc_html = expand("results/QC/alignment/qualimap/bamqc/{sample}/qualimapReport.html", sample = config["sample"]),
        qmap_bamqc_txt = expand("results/QC/alignment/qualimap/bamqc/{sample}/genome_results.txt", sample = config["sample"]),
        qmap_rnaseq_html = expand("results/QC/alignment/qualimap/rnaseq/{sample}/qualimapReport.html", sample = config["sample"]),
        qmap_rnaseq_txt = expand("results/QC/alignment/qualimap/rnaseq/{sample}/rnaseq_qc_results.txt", sample = config["sample"]),
        samtools_stats = expand("results/QC/alignment/samtools_stats/{sample}_Aligned.sortedByCoord.out.bam.stats", sample = config["sample"]),
        samtools_flagstat = expand("results/QC/alignment/samtools_stats/{sample}_Aligned.sortedByCoord.out.bam.flagstat",  sample = config["sample"]),
        rseqc_strandiness = expand("results/QC/alignment/rseqc/{sample}_strandiness.txt", sample = config["sample"])
    output:
        multiqc = "results/MultiQC/multiqc_report.html"
    conda: 
        config["conda_envs"]["qc"]
    params: 
        outdir = lambda wildcards, output : os.path.dirname(output.multiqc)
    log:
        log = "log/MultiQC/multiqc_report.log",
    shell: 
        "multiqc {input} -o {params.outdir} 2> {log.log} "

## Let stablish specific rules to deal with UMIs.
if UMIs:

    rule umi_extract:
        input:
            "data/{sample}_{seq_lane}.fastq.gz"
        output:
            "results/umi_extract/{sample}_{seq_lane}.fastq.gz"
        conda:
            config["conda_envs"]["umi_tools"]
        threads: 3
        resources:
            mem_mb=15000
        params:
            pattern = lambda wildcards: config["umi_processing"]["pattern"]
        log:
            "log/umi_extract/{sample}_{seq_lane}.log"
        benchmark:
            "benchmarks/{sample}_{seq_lane}_umi_tools_extract.bmk"
        shell: """
            umi_tools extract --stdin={input} \
                --extract-method regex --bc-pattern={params.pattern} \
                --log={log} --stdout={output}
        """

    rule umi_dedup:
        input:
            bam = "results/alignment/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
            bam_bai = "results/alignment/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
        output:
            dedup = "results/alignment/dedup/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
        conda:
            config["conda_envs"]["umi_tools"]
        threads: 3
        resources:
            mem_mb=10000
        params:
            stats = "results/dedup/alignments/{sample}"
        log:
            "log/dedup/{sample}.log"
        benchmark:
            "benchmarks/{sample}_umi_tools_dedup.bmk"
        shell:"""
            umi_tools dedup -I {input.bam} --log={log} -S {output.dedup} --output-stats={params.stats} \
                --method=unique --multimapping-detection-method=NH
        """