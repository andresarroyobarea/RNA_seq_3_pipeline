configfile: "config/config.yaml"

# ---- FastQC and MultiQC results ------
#bam_sorted = expand("results/alignment/{sample}_{condition}/{sample}_{condition}Aligned.sortedByCoord.out.bam", sample = config["sample"], condition = config["condition"]),
feature_table = "results/feature_counts/counts_raw.tsv" 
multiqc_report = "results/MultiQC/multiqc_report.html"

rule all:
    input:
        multiqc_report,
        feature_table

rule fastqc_raw_trim:
    input:
        "data/{seqs_state}/{sample}_{condition}_{seq_lane}_R1_001_{seqs_state}.fastq.gz",
    output:
        html = "results/QC/qc_per_lane/{seqs_state}/{sample}_{condition}_{seq_lane}_R1_001_{seqs_state}_fastqc.html",
        zip = "results/QC/qc_per_lane/{seqs_state}/{sample}_{condition}_{seq_lane}_R1_001_{seqs_state}_fastqc.zip",
    log:
        log = "log/QC/qc_per_lane/{seqs_state}/{sample}_{condition}_{seq_lane}_R1_001_fastqc_{seqs_state}.log",
    params: 
        outdir = "results/QC/qc_per_lane/{seqs_state}",
    conda:
        config["conda_envs"]["qc"]
    threads: 2
    shell:
        "mkdir -p {params.outdir} &&"
        "fastqc --outdir {params.outdir} --threads {threads} {input} 2> {log.log} "

rule bbduk_se:
    input:
        sample = "data/raw/{sample}_{condition}_{seq_lane}_R1_001_raw.fastq.gz",
        adapters = "resources/trim_files/adapters.fa.gz",
        polyA = "resources/trim_files/polyA.fa.gz"
    output:
        trimmed = "data/trimmed/{sample}_{condition}_{seq_lane}_R1_001_trimmed.fastq.gz",
        singleton = "data/trimmed/{sample}_{condition}_{seq_lane}_R1_001_single.fastq.gz",
        discarded = "data/trimmed/{sample}_{condition}_{seq_lane}_R1_001_discarded.fastq.gz",
        stats = "data/trimmed/{sample}_{condition}_{seq_lane}_R1_001_stats.txt",
    log:
        "log/bbduk/{sample}_{condition}_{seq_lane}_R1_001.log"
    conda: 
        config["conda_envs"]["rna_seq_3"]
    threads: 2
    shell:
        """
        bbduk.sh in={input.sample} \
            out={output.trimmed} \
            outs={output.singleton} \
            outm={output.discarded} \
            stats={output.stats} \
            ref={input.adapters},{input.polyA} \
            k=13 ktrim=r mink=5 qtrim=r trimq=20 minlength=20 > {log} 2>&1
        """

rule merged_fastq:
    input: 
        lambda wildcards: expand(
            "data/trimmed/{sample}_{condition}_{seq_lane}_R1_001_trimmed.fastq.gz",
            sample = wildcards.sample,
            condition = wildcards.condition,
            seq_lane = config["seq_lane"]
        )
    output: 
        fastq_merged = "data/merged/{sample}_{condition}_merged.fastq.gz"
    log:
        "log/merged/{sample}_{condition}_merged.log"
    shell:
        "cat {input} > {output.fastq_merged} 2> {log}"

# rule fastq_screen:
#    input: 
#        fastq_merged = "data/merged/{sample}_{condition}_merged.fastq.gz"
#    output: 
#        fastq_screen_txt = "results/QC/fastq_screen/{sample}_{condition}/{sample}_{condition}_merged_screen.txt",
#        fastq_screen_png = "results/QC/fastq_screen/{sample}_{condition}/{sample}_{condition}_merged_screen.png"
#    params:
#        fastq_screen_config = config["fastq_screen_conf"],
#        aligner = config["fastq_screen_aling"],
#        outdir = "results/QC/fastq_screen/{sample}_{condition}"
#    conda:
#        config["conda_envs"]["rna_seq_3_v2"]
#    threads: 1
#    resources:
#        mem_mb=28728
#    log:
#        log = "log/QC/fastq_screen/{sample}_{condition}_fastq_screen.log"
#    shell:"""
#        fastq_screen {input.fastq_merged} --aligner {params.aligner} \
#            --conf {params.fastq_screen_config} --outdir {params.outdir} \
#            -threads {threads} 2> {log}
#    """

rule alignment:
    input: 
        "data/merged/{sample}_{condition}_merged.fastq.gz"
    output: 
        bam_sorted = "results/alignment/{sample}_{condition}/{sample}_{condition}Aligned.sortedByCoord.out.bam"
    params:
        genome_index = config["genome_index"],
        out_name = "results/alignment/{sample}_{condition}/{sample}_{condition}"
    conda:
        config["conda_envs"]["rna_seq_3"]
    threads: 5
    resources:
        mem_mb=18432
    log:
        "log/alignment/{sample}_{condition}_STAR.log" 
    shell: """
        STAR --runThreadN {threads} --genomeDir {params.genome_index} --genomeLoad LoadAndKeep --readFilesIn {input} \
            --readFilesCommand gunzip -c --outFilterType BySJout --outFilterMultimapNmax 25 --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.3 \
            --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --limitBAMsortRAM 12000000000 \
            --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.out_name} 2> {log}
    """

rule fastqc_alignment:
    input:
        "results/alignment/{sample}_{condition}/{sample}_{condition}Aligned.sortedByCoord.out.bam"
    output:
        html = "results/QC/alignment/FastQC/{sample}_{condition}Aligned.sortedByCoord.out_fastqc.html",
        zip = "results/QC/alignment/FastQC/{sample}_{condition}Aligned.sortedByCoord.out_fastqc.zip",
    log:
        log = "log/QC/alignment/FastQC/{sample}_{condition}_bam_fastqc.log"
    params: 
        outdir = "results/QC/alignment/FastQC"
    conda:
        config["conda_envs"]["qc"]
    threads: 5
    shell:
        "mkdir -p {params.outdir} &&"
        "fastqc --outdir {params.outdir} --threads {threads} {input} 2> {log.log} "

# rule qualimap_bamqc:
#    input:
#        bam_sorted = "results/alignment/{sample}_{condition}/{sample}_{condition}Aligned.sortedByCoord.out.bam"
#    output:
#        "results/QC/alignment/qualimap/bamqc/{sample}_{condition}/qualimapReport.html",
#        "results/QC/alignment/qualimap/bamqc/{sample}_{condition}/genome_results.txt",
#    params:
#        qmap_genome = config["qualimap"]["genome"],
#        annotation = config["annotation"],
#        outdir = "results/QC/alignment/qualimap/bamqc/{sample}_{condition}",
#        mem = config["qualimap"]["mem"]
#    resources:
#        mem_mb=26624
#    log: 
#        "log/QC/alignment/qualimap/bamqc/{sample}_{condition}_qualimap_bamqc.log"
#    conda:
#        config["conda_envs"]["rna_seq_3_v2"]
#    threads: 3
#    shell: """
#        qualimap bamqc -bam {input.bam_sorted} -gd {params.qmap_genome} -gff {params.annotation} \
#            -hm 3 -nr 1000 -nt {threads} --outdir {params.outdir} -p strand-specific-forward \
#            --java-mem-size={params.mem} 2> {log}
#    """

#rule qualimap_rnaseq:
#    input:  
#        bam_sorted = "results/alignment/{sample}_{condition}/{sample}_{condition}Aligned.sortedByCoord.out.bam"
#    output: 
#        "results/QC/alignment/qualimap/rnaseq/{sample}_{condition}/qualimapReport.html",
#        "results/QC/alignment/qualimap/rnaseq/{sample}_{condition}/rnaseq_qc_results.txt",
#    params:
#        annotation = config["annotation"],
#        outdir = "results/QC/alignment/qualimap/rnaseq/{sample}_{condition}",
#        mem = config["qualimap"]["mem"]
#    log:
#        "log/QC/alignment/qualimap/rnaseq/{sample}_{condition}_qualiamp_rnaseq.log"
#    conda:
#        config["conda_envs"]["rna_seq_3_v2"]
#    shell:"""
#        qualimap rnaseq -bam {input.bam_sorted} -gtf {params.annotation} \
#            -outdir {params.outdir} -p strand-specific-forward \
#            --java-mem-size={params.mem} 2> {log}
#    """

rule rseqc_strand:
    input: 
        bam_sorted = "results/alignment/{sample}_{condition}/{sample}_{condition}Aligned.sortedByCoord.out.bam"
    output: 
        "results/QC/alignment/rseqc/{sample}_{condition}_strandiness.txt"
    params:
        bed_file = config["reference_bed"]
    conda: 
        config["conda_envs"]["rna_seq_3_v2"]
    log:
        "log/QC/alignment/rseqc/{sample}_{condition}_strandiness.log"
    shell:
        "infer_experiment.py -i {input.bam_sorted} -r {params.bed_file} > {output} 2> {log} "

rule samtools_stats_flagstat:
    input:
        bam_sorted = "results/alignment/{sample}_{condition}/{sample}_{condition}Aligned.sortedByCoord.out.bam"
    output:
        samtools_stats = "results/QC/alignment/samtools_stats/{sample}_{condition}Aligned.sortedByCoord.out.bam.stats",
        samtools_flagstat = "results/QC/alignment/samtools_stats/{sample}_{condition}Aligned.sortedByCoord.out.bam.flagstat"
    conda:
        config["conda_envs"]["rna_seq_3"]
    log:
        log_stats = "log/QC/alignment/samtools/{sample}_{condition}_samtools_stats.log",
        log_flagstat = "log/QC/alignment/samtools/{sample}_{condition}_samtools_flagstat.log"
    shell:"""
        samtools stats {input.bam_sorted} > {output.samtools_stats} 2> {log.log_stats} &&
        samtools flagstat {input.bam_sorted} > {output.samtools_flagstat} 2> {log.log_flagstat}
    """

rule feature_counts:
    input:
        bam_sorted = expand("results/alignment/{sample}_{condition}/{sample}_{condition}Aligned.sortedByCoord.out.bam", sample = config["sample"], condition = config["condition"])
    output: 
        feature_table = "results/feature_counts/counts_raw.tsv" 
    params:
        annotations = config["annotation"]
    conda: 
        config["conda_envs"]["rna_seq_3"]
    log:
        "log/featureCounts/counts_raw.log"
    shell:"""
        featureCounts -a {params.annotations} -O -F GTF -t gene -g gene_id \
            --extraAttributes gene_name,transcript_name -s 1 -T 15 \
            -o {output.feature_table} {input.bam_sorted} 2> {log}
    """

rule multiqc:
    input: 
        seqs_QC = expand("results/QC/qc_per_lane/{seqs_state}/{sample}_{condition}_{seq_lane}_R1_001_{seqs_state}_fastqc.html", sample = config["sample"], condition = config["condition"], seq_lane = config["seq_lane"], seqs_state = config["seqs_state"]),
        #fastq_screen_txt = expand("results/QC/fastq_screen/{sample}_{condition}/{sample}_{condition}_merged_screen.txt", sample = config["sample"], condition = config["condition"]),
        #fastq_screen_png = expand("results/QC/fastq_screen/{sample}_{condition}/{sample}_{condition}_merged_screen.png", sample = config["sample"], condition = config["condition"]),
        bam_QC = expand("results/QC/alignment/FastQC/{sample}_{condition}Aligned.sortedByCoord.out_fastqc.html", sample = config["sample"], condition = config["condition"]),
        #qmap_bamqc_html = expand("results/QC/alignment/qualimap/bamqc/{sample}_{condition}/qualimapReport.html", sample = config["sample"], condition = config["condition"]),
        #qmap_bamqc_txt = expand("results/QC/alignment/qualimap/bamqc/{sample}_{condition}/genome_results.txt", sample = config["sample"], condition = config["condition"]),
        #qmap_rnaseq_html = expand("results/QC/alignment/qualimap/rnaseq/{sample}_{condition}/qualimapReport.html", sample = config["sample"], condition = config["condition"]),
        #qmap_rnaseq_txt = expand("results/QC/alignment/qualimap/rnaseq/{sample}_{condition}/rnaseq_qc_results.txt", sample = config["sample"], condition = config["condition"]),
        samtools_stats = expand("results/QC/alignment/samtools_stats/{sample}_{condition}Aligned.sortedByCoord.out.bam.stats", sample = config["sample"], condition = config["condition"]),
        samtools_flagstat = expand("results/QC/alignment/samtools_stats/{sample}_{condition}Aligned.sortedByCoord.out.bam.flagstat",  sample = config["sample"], condition = config["condition"]),
        rseqc_strandiness = expand("results/QC/alignment/rseqc/{sample}_{condition}_strandiness.txt", sample = config["sample"], condition = config["condition"])
    output:
        output = "results/MultiQC/multiqc_report.html"
    params: 
        outdir = "results/MultiQC",
    conda: 
        config["conda_envs"]["qc"]
    log:
        log = "log/MultiQC/multiqc_report.log",
    shell: 
        "multiqc {input} -o {params.outdir} 2> {log.log} "

