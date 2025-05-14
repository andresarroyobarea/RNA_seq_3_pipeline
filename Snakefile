configfile: "config/config.yaml"

import os

UMIs = config["umi_processing"]["enabled"]

rule all:
    input:
        "results/feature_counts/counts_raw.tsv",
        "results/QC/alignment/qualimap/multi_bamqc/multisampleBamQcReport.html",
        "results/MultiQC/multiqc_report.html"

rule fastqc:
    input: 
        fastq = "data/{seqs_state}/{sample}_{seq_lane}_{seqs_state}.fastq.gz"
    output:
        html = "results/QC/qc_per_lane/{seqs_state}/{sample}_{seq_lane}_{seqs_state}_fastqc.html",
        zip = "results/QC/qc_per_lane/{seqs_state}/{sample}_{seq_lane}_{seqs_state}_fastqc.zip",
    log:
        fastqc = "log/QC/qc_per_lane/{seqs_state}/{sample}_{seq_lane}_{seqs_state}_fastqc.log",
    params: 
        outdir = lambda wildcards, output: os.path.dirname(output.html)
    conda:
        config["conda_envs"]["qc"]
    threads: 2
    shell:
        "mkdir -p {params.outdir} &&"
        "fastqc --outdir {params.outdir} --threads {threads} {input.fastq} 2> {log.fastqc} "

rule bbduk_se:
    input:
        sample = ["data/ + dir_in_trim + /{sample}_{seq_lane}_raw.fastq.gz"]
    output:
        trimmed = "data/trimmed/{sample}/{sample}_{seq_lane}_trimmed.fastq.gz",
        singleton = "data/trimmed/{sample}/{sample}_{seq_lane}_single.fastq.gz",
        discarded = "data/trimmed/{sample}/{sample}_{seq_lane}_discarded.fastq.gz",
        stats = "data/trimmed/{sample}/{sample}_{seq_lane}_stats.txt",
    log:
        "log/bbduk/{sample}_{seq_lane}_bbduk.log"
    conda: 
        config["conda_envs"]["rna_seq_3"]
    threads: 2
    params:
        adapters = "resources/trim_files/adapters.fa.gz",
        polyA = "resources/trim_files/polyA.fa.gz"
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

rule merged_fastq:
    input: 
        lambda wildcards: expand(
            "data/trimmed/{sample}/{sample}_{seq_lane}_trimmed.fastq.gz",
            sample = wildcards.sample,
            seq_lane = config["seq_lane"]
        )
    output: 
        fastq_merged = "data/merged/{sample}.fastq.gz"
    conda: 
        config["conda_envs"]["rna_seq_3"]
    log:
        "log/merged/{sample}_merging.log"
    shell:
        "cat {input} > {output.fastq_merged} 2> {log}"

rule fastq_screen:
    input: 
        fastq_merged = expand("data/merged/{sample}.fastq.gz", sample = config["sample"])
    output: 
        fastq_screen_txt = "results/QC/fastq_screen/{sample}_fastq_screen.txt",
        fastq_screen_png = "results/QC/fastq_screen/{sample}_fastq_screen.png"
    params:
        fastq_screen_config = config["fastq_screen_conf"],
        aligner = config["fastq_screen_aling"],
        outdir = "results/QC/fastq_screen/"
    conda:
        config["conda_envs"]["rna_seq_3_v2"]
    threads: 1
    resources:
        mem_mb=28728
    log:
        log = "log/QC/fastq_screen/{sample}_fastq_screen.log"
    shell:"""
        fastq_screen {input.fastq_merged} --aligner {params.aligner} \
            --conf {params.fastq_screen_config} --outdir {params.outdir} \
            -threads {threads} 2> {log}
    """

rule alignment:
    input: 
        fastq_merged = "data/merged/{sample}.fastq.gz"
    output: 
        bam = "results/alignment/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    params:
        genome_index = config["genome_index"],
        outdir = lambda wildcards, output: os.path.dirname(output.bam)
    conda:
        config["conda_envs"]["rna_seq_3"]
    threads: 5
    resources:
        mem_mb=18432
    log:
        "log/alignment/{sample}_alignment.log" 
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
    log:
        "log/bam_indexing/{sample}.log"
    threads: 3
    resources:
        mem_mb=3
    conda:
        config["conda_envs"]["rna_seq_3"]
    shell:
        "samtools index -@ {threads} {input.bam} "

rule fastqc_alignment:
    input:
        bam = "results/alignment/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        html = "results/QC/alignment/FastQC/{sample}_Aligned.sortedByCoord.out_fastqc.html",
        zip = "results/QC/alignment/FastQC/{sample}_Aligned.sortedByCoord.out_fastqc.zip",
    log:
        "log/QC/alignment/FastQC/{sample}_alignment_fastqc.log"
    params: 
        outdir = lambda wildcards, output : os.path.dirname(output.html)
    conda:
        config["conda_envs"]["qc"]
    threads: 5
    shell:
        "mkdir -p {params.outdir} &&"
        "fastqc --outdir {params.outdir} --threads {threads} {input.bam} 2> {log} "

rule qualimap_bamqc:
    input:
        bam = "results/alignment/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        qmap_report = "results/QC/alignment/qualimap/bamqc/{sample}/qualimapReport.html",
        genome_res = "results/QC/alignment/qualimap/bamqc/{sample}/genome_results.txt",
    params:
        qmap_genome = config["qualimap"]["genome"],
        annotation = config["annotation"],
        outdir = "results/QC/alignment/qualimap/bamqc/{sample}",
        mem = config["qualimap"]["mem"]
    resources:
        mem_mb=26624
    log: 
        "log/QC/alignment/qualimap/bamqc/{sample}_qualimap_bamqc.log"
    conda:
        config["conda_envs"]["rna_seq_3_v2"]
    threads: 3
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
    params: 
        qmap_input = "metadata/qualimap_multi_bamqc_input.txt",
        outdir = "results/QC/alignment/qualimap/multi_bamqc"
    resources:
        mem_mb=26624
    log: 
        "log/QC/alignment/qualimap/mutli_bamqc/qualimap_multi_bamqc.log"
    conda:
        config["conda_envs"]["rna_seq_3_v2"]
    threads: 3    
    shell:
        "qualimap multi-bamqc -d {params.qmap_input} --outdir {params.outdir} 2> {log} "


rule qualimap_rnaseq:
    input:  
        bam = "results/alignment/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    output: 
        "results/QC/alignment/qualimap/rnaseq/{sample}/qualimapReport.html",
        "results/QC/alignment/qualimap/rnaseq/{sample}/rnaseq_qc_results.txt"
    params:
        annotation = config["annotation"],
        outdir = "results/QC/alignment/qualimap/rnaseq/{sample}",
        mem = config["qualimap"]["mem"]
    log:
        "log/QC/alignment/qualimap/rnaseq/{sample}_qualiamp_rnaseq.log"
    conda:
        config["conda_envs"]["rna_seq_3_v2"]
    shell:"""
        qualimap rnaseq -bam {input.bam} -gtf {params.annotation} \
            -outdir {params.outdir} -p strand-specific-forward \
            --java-mem-size={params.mem} 2> {log}
    """

rule rseqc_strand:
    input: 
        bam = "results/alignment/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    output: 
        rseqc_out = "results/QC/alignment/rseqc/{sample}_strandiness.txt"
    params:
        bed_file = config["reference_bed"]
    conda: 
        config["conda_envs"]["rna_seq_3_v2"]
    log:
        "log/QC/alignment/rseqc/{sample}_rseqc.log"
    shell:
        "infer_experiment.py -i {input.bam} -r {params.bed_file} > {output.rseqc_out} 2> {log} "

rule samtools_stats_flagstat:
    input:
        bam = "results/alignment/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        samtools_stats = "results/QC/alignment/samtools_stats/{sample}_Aligned.sortedByCoord.out.bam.stats",
        samtools_flagstat = "results/QC/alignment/samtools_stats/{sample}_Aligned.sortedByCoord.out.bam.flagstat"
    conda:
        config["conda_envs"]["rna_seq_3"]
    log:
        log_stats = "log/QC/alignment/samtools/{sample}_samtools_stats.log",
        log_flagstat = "log/QC/alignment/samtools/{sample}_samtools_flagstat.log"
    shell:"""
        samtools stats {input.bam} > {output.samtools_stats} 2> {log.log_stats} &&
        samtools flagstat {input.bam} > {output.samtools_flagstat} 2> {log.log_flagstat}
    """

rule feature_counts:
    input:
        bam = expand("results/alignment/{sample}/{sample}_Aligned.sortedByCoord.out.bam", sample = config["sample"])
    output: 
        feature_table = "results/feature_counts/counts_raw.tsv" 
    params:
        annotations = config["annotation"]
    conda: 
        config["conda_envs"]["rna_seq_3"]
    log:
        "log/featureCounts/featureCounts.log"
    shell:"""
        featureCounts -a {params.annotations} -O -F GTF -t gene -g gene_id \
            --extraAttributes gene_name,transcript_name -s 1 -T 15 \
            -o {output.feature_table} {input.bam} 2> {log}
    """

rule multiqc:
    input: 
        seqs_QC = expand("results/QC/qc_per_lane/{seqs_state}/{sample}_{seq_lane}_{seqs_state}_fastqc.html", sample = config["sample"], seq_lane = config["seq_lane"], seqs_state = config["seqs_state"]),
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
    params: 
        outdir = lambda wildcards, output : os.path.dirname(output.multiqc)
    conda: 
        config["conda_envs"]["qc"]
    log:
        log = "log/MultiQC/multiqc_report.log",
    shell: 
        "multiqc {input} -o {params.outdir} 2> {log.log} "

## Let stablish where the trimming input will be saved.
if UMIs:
    dir_in_trim = "umi_extract"
else:
    dir_in_trim = "raw"




## Let stablish specific rules to deal with UMIs.
if UMIs:

    rule umi_tools_extract:
        input:
            "data/raw/{sample}_{seq_lane}_raw.fastq.gz"
        output:
            "data/umi_extract/{sample}_{seq_lane}_raw.fastq.gz"
        conda:
            config["conda_envs"]["rnrna_seq_3_v2"]
        threads: 3
        resources:
            mem_mb=15000
        params:
            config["umi_processing"]["pattern"]
        log:
            "log/umi_extract/{sample}_{seq_lane}.log"
        shell: """
            umi_tools extract --stdin={input.fastq} \
                --extract-method regex --bc-pattern={params} \
                --log={log} --stdout={output}
        """



# AÑADIR REGLA CONDICIONAL DE UMI-TOOLS ---> UMI-TOOLS DEDUPLICATION AFTER ALIGMENT.
#rule umi_tools_dedup: