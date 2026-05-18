rule fastqc_raw:
    input: 
        fastq = "data/{sample}_{seq_lane}.fastq.gz"
    output:
        html = "results/QC/raw/{sample}/{sample}_{seq_lane}_fastqc.html",
        zip = "results/QC/raw/{sample}/{sample}_{seq_lane}_fastqc.zip"
    conda:
        config["conda_envs"]["qc"]
    threads: 
        get_resource(config, "fastqc", "threads")
    resources:
        mem_mb = get_resource(config, "fastqc", "mem_mb"),
        runtime = get_resource(config, "fastqc", "runtime")
    params: 
        outdir = lambda wildcards, output: os.path.dirname(output.html),
        extra = config["parameters"]["fastqc"]["extra"]
    log:
        "log/QC/raw/fastqc/{sample}_{seq_lane}.log"
    benchmark:
        "benchmarks/QC/raw/fastqc/{sample}_{seq_lane}.bmk"
    shell: """
        mkdir -p {params.outdir} &&
        fastqc --outdir {params.outdir} \
            --threads {threads} \
            {input.fastq} \
            {params.extra} 2> {log} 
    """

rule fastqc_trim:
    input: 
        fastq = "results/trimmed/{sample}/{sample}_{seq_lane}_trimmed.fastq.gz"
    output:
        html = "results/QC/trimmed/fastqc/{sample}/{sample}_{seq_lane}_trimmed_fastqc.html",
        zip = "results/QC/trimmed/fastqc/{sample}/{sample}_{seq_lane}_trimmed_fastqc.zip",
    conda:
        config["conda_envs"]["qc"]
    threads: 
        get_resource(config, "fastqc", "threads")
    resources:
        mem_mb = get_resource(config, "fastqc", "mem_mb"),
        runtime = get_resource(config, "fastqc", "runtime")
    params: 
        outdir = lambda wildcards, output: os.path.dirname(output.html),
        extra = config["parameters"]["fastqc"]["extra"]
    log:
        "log/QC/trimmed/fastqc/{sample}_{seq_lane}.log"
    benchmark:
        "benchmarks/QC/trimmed/fastqc/{sample}_{seq_lane}.bmk"
    shell: """
        mkdir -p {params.outdir} &&
        fastqc --outdir {params.outdir} \
            --threads {threads} \
            {input.fastq} \
            {params.extra} 2> {log} 
    """


rule fastq_screen_files:
    input: 
        fastq = "results/trimmed/{sample}/{sample}_{seq_lane}_trimmed.fastq.gz"
    output: 
        fastq_screen_txt = "results/QC/trimmed/fastq_screen/{sample}/{sample}_{seq_lane}_trimmed_screen.txt",
        fastq_screen_png = "results/QC/trimmed/fastq_screen/{sample}/{sample}_{seq_lane}_trimmed_screen.png"
    conda:
        config["conda_envs"]["fastq_screen"]
    threads: 
        get_resource(config, "fastq_screen", "threads")
    resources:
        mem_mb = get_resource(config, "fastq_screen", "mem_mb"),
        runtime = get_resource(config, "fastq_screen", "runtime")
    params: # TODO: Config as input
        fscreen_config = config["parameters"]["fastq_screen"]["config"],
        aligner = config["parameters"]["fastq_screen"]["aligner"],
        outdir = lambda wildcards, output: os.path.dirname(output.fastq_screen_txt),
        extra = config["parameters"]["fastq_screen"]["extra"]
    log:
        "log/QC/trimmed/fastq_screen/{sample}_{seq_lane}.log"
    benchmark:
        "benchmarks/QC/trimmed/fastq_screen/{sample}_{seq_lane}.bmk"
    shell:"""
        fastq_screen {input.fastq} \
            --aligner {params.aligner} \
            --conf {params.fscreen_config} \
            --outdir {params.outdir} \
            -threads {threads} 2> {log}
    """

rule fastqc_merged:
    input: 
        fastq = "results/merged/{sample}.fastq.gz"
    output:
        html = "results/QC/merged/FastQC/{sample}/{sample}_fastqc.html",
        zip = "results/QC/merged/FastQC/{sample}/{sample}_fastqc.zip"
    conda:
        config["conda_envs"]["qc"]
    threads: 
        get_resource(config, "fastqc", "threads")
    resources:
        mem_mb = get_resource(config, "fastqc", "mem_mb"),
        runtime = get_resource(config, "fastqc", "runtime")
    params: 
        outdir =lambda wildcards, output: os.path.dirname(output.html),
        extra = config["parameters"]["fastqc"]["extra"]
    log:
        "log/QC/merged/fastqc/{sample}.log"
    benchmark:
        "benchmarks/QC/merged/fastqc/{sample}.bmk"
    shell: """
        mkdir -p {params.outdir} &&
        fastqc --outdir {params.outdir} \
            --threads {threads} \
            {input.fastq} \
            {params.extra} 2> {log} 
    """

rule fastq_screen_merged:
    input: 
        fastq_merged = "results/merged/{sample}.fastq.gz"
    output: 
        fastq_screen_txt = "results/QC/merged/fastq_screen/{sample}/{sample}_screen.txt",
        fastq_screen_png = "results/QC/merged/fastq_screen/{sample}/{sample}_screen.png",
        #fastq_screen_html = "results/QC/merged/fastq_screen/{sample}/{sample}_screen.html"
    conda:
        config["conda_envs"]["fastq_screen"]
    threads: 
        get_resource(config, "fastq_screen", "threads")
    resources:
        mem_mb = get_resource(config, "fastq_screen", "mem_mb"),
        runtime = get_resource(config, "fastq_screen", "runtime")
    params:
        fscreen_config = config["parameters"]["fastq_screen"]["config"],
        aligner = config["parameters"]["fastq_screen"]["aligner"],
        outdir = lambda wildcards, output: os.path.dirname(output.fastq_screen_txt),
        extra = config["parameters"]["fastq_screen"]["extra"]
    log:
        "log/QC/merged/fastq_screen/{sample}.log"
    benchmark:
        "benchmarks/QC/merged/fastq_screen/{sample}.bmk"
     shell:"""
        fastq_screen {input.fastq} \
            --aligner {params.aligner} \
            --conf {params.fscreen_config} \
            --outdir {params.outdir} \
            -threads {threads} 2> {log}
    """


rule fastqc_alignment:
    input:
        bam = "results/alignment/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        html = "results/QC/alignment/FastQC/{sample}/{sample}_Aligned.sortedByCoord.out_fastqc.html",
        zip = "results/QC/alignment/FastQC/{sample}/{sample}_Aligned.sortedByCoord.out_fastqc.zip"
    conda:
        config["conda_envs"]["qc"]
    threads: get_resource(config, "fastqc", "threads")
    resources:
        mem_mb = get_resource(config, "fastqc", "mem_mb"),
        runtime = get_resource(config, "fastqc", "runtime")
    params: 
        outdir = lambda wildcards, output : os.path.dirname(output.html)
    log:
        "log/QC/alignment/FastQC/{sample}_alignment_fastqc.log"
    shell:
        "mkdir -p {params.outdir} &&"
        "fastqc --outdir {params.outdir} --threads {threads} {input.bam} 2> {log} "

rule qualimap_bamqc:
    input:
        bam = "results/alignment/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        qmap_report = "results/QC/alignment/qualimap/bamqc/{sample}/qualimapReport.html",
        genome_res = "results/QC/alignment/qualimap/bamqc/{sample}/genome_results.txt"
    conda:
        config["conda_envs"]["qualimap"]
    threads: get_resource(config, "qualimap", "threads")
    resources:
        mem_mb = get_resource(config, "qualimap", "mem_mb"),
        runtime = get_resource(config, "qualimap", "runtime")
    params:
        qmap_genome = config["qualimap"]["genome"],
        annotation = config["annotation"],
        outdir = lambda wildcards, output: os.path.dirname(output.qmap_report),
        mem = f"{get_resource(config, 'qualimap', 'mem_mb') // 1024}G"
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
        qmap_report = "results/QC/alignment/qualimap/multi_bamqc/multisampleBamQcReport.html"
    conda:
        config["conda_envs"]["qualimap"]
    threads: get_resource(config, "qualimap", "threads")
    resources:
        mem_mb = get_resource(config, "qualimap", "mem_mb"),
        runtime = get_resource(config, "qualimap", "runtime")
    params: 
        qmap_input = "metadata/qualimap_multi_bamqc_input.txt",
        outdir = lambda wildcards, output: os.path.dirname(output.qmap_report)
    log: 
        "log/QC/alignment/qualimap/mutli_bamqc/qualimap_multi_bamqc.log"
    benchmark:
        "benchmarks/qualimap_multi_bamqc.bmk"
    shell:
        "qualimap multi-bamqc -d {params.qmap_input} --outdir {params.outdir} 2> {log} "


rule qualimap_rnaseq:
    input:  
        bam = "results/alignment/{sample}_Aligned.sortedByCoord.out.bam"
    output: 
        qmap_report = "results/QC/alignment/qualimap/rnaseq/{sample}/qualimapReport.html",
        qmap_res = "results/QC/alignment/qualimap/rnaseq/{sample}/rnaseq_qc_results.txt"
    conda:
        config["conda_envs"]["qualimap"]
    threads: get_resource(config, "qualimap", "threads")
    resources:
        mem_mb = get_resource(config, "qualimap", "mem_mb"),
        runtime = get_resource(config, "qualimap", "runtime")
    params:
        annotation = config["annotation"],
        outdir = lambda wildcards, output: os.path.dirname(output.qmap_report),
        mem = f"{get_resource(config, 'qualimap', 'mem_mb') // 1024}G"
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
        bam = "results/alignment/{sample}_Aligned.sortedByCoord.out.bam"
    output: 
        rseqc_out = "results/QC/alignment/rseqc/{sample}_strandiness.txt"
    conda: 
        config["conda_envs"]["rseqc"]
    threads: get_resource(config, "rseqc_strand", "threads")
    resources:
        mem_mb = get_resource(config, "rseqc_strand", "mem_mb"),
        runtime = get_resource(config, "rseqc_strand", "runtime")
    params:
        bed_file = config["reference_bed"]
    log:
        "log/QC/alignment/rseqc/{sample}_rseqc.log"
    shell:
        "infer_experiment.py -i {input.bam} -r {params.bed_file} > {output.rseqc_out} 2> {log} "

rule samtools_stats_flagstat:
    input:
        bam = "results/alignment/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        samtools_stats = "results/QC/alignment/samtools_stats/{sample}_Aligned.sortedByCoord.out.bam.stats",
        samtools_flagstat = "results/QC/alignment/samtools_stats/{sample}_Aligned.sortedByCoord.out.bam.flagstat"
    conda:
        config["conda_envs"]["aligners"]
    threads: get_resource(config, "default", "threads")
    resources:
        mem_mb = get_resource(config, "default", "mem_mb"),
        runtime = get_resource(config, "default", "runtime")
    log:
        log_stats = "log/QC/alignment/samtools/{sample}_samtools_stats.log",
        log_flagstat = "log/QC/alignment/samtools/{sample}_samtools_flagstat.log"
    shell:"""
        samtools stats {input.bam} > {output.samtools_stats} 2> {log.log_stats} &&
        samtools flagstat {input.bam} > {output.samtools_flagstat} 2> {log.log_flagstat}
    """

rule multiqc_trimmed:
    input: 
        seqs_QC_trim = expand("results/QC/trimmed/fastqc/{sample}/{sample}_{seq_lane}_trimmed_fastqc.html", sample = config["sample"], seq_lane = config["seq_lane"]),
        fastq_screen_txt = expand("results/QC/trimmed/fastq_screen/{sample}/{sample}_{seq_lane}_trimmed_screen.txt", sample = config["sample"], seq_lane = config["seq_lane"]),
        fastq_screen_png = expand("results/QC/trimmed/fastq_screen/{sample}/{sample}_{seq_lane}_trimmed_screen.png", sample = config["sample"], seq_lane = config["seq_lane"])
    output:
        multiqc = "results/QC/MultiQC/trimmed/multiqc_report.html"
    conda: 
        config["conda_envs"]["qc"]
    threads: get_resource(config, "default", "threads")
    resources:
        mem_mb = get_resource(config, "default", "mem_mb"),
        runtime = get_resource(config, "default", "runtime")
    params:
        inputdir = ["results/QC/trimmed", "results/trimmed", "log/bbduk"],
        outdir = lambda wildcards, output : os.path.dirname(output.multiqc)
    log:
        log = "log/QC/MultiQC/trimmed/multiqc_report.log",
    shell: 
        "multiqc {params.inputdir} -o {params.outdir} 2> {log.log} "

rule multiqc_merge:
    input:
        seq_QC_merged = expand("results/QC/merged/FastQC/{sample}/{sample}_fastqc.html", sample = config["sample"]),
        fastq_screen_merged_txt = expand("results/QC/merged/fastq_screen/{sample}/{sample}_screen.txt", sample = config["sample"]),
        fastq_screen_merged_png = expand("results/QC/merged/fastq_screen/{sample}/{sample}_screen.png", sample = config["sample"]),
        bam_QC = expand("results/QC/alignment/FastQC/{sample}/{sample}_Aligned.sortedByCoord.out_fastqc.html", sample = config["sample"]),
        umi_extract = expand("log/umitools/extract/{sample}_{seq_lane}.log", sample = config["sample"], seq_lane = config["seq_lane"]),
        umi_dedup = expand("log/umitools/dedup/{sample}.log", sample = config["sample"]),
        qmap_bamqc_html = expand("results/QC/alignment/qualimap/bamqc/{sample}/qualimapReport.html", sample = config["sample"]),
        qmap_bamqc_txt = expand("results/QC/alignment/qualimap/bamqc/{sample}/genome_results.txt", sample = config["sample"]),
        qmap_rnaseq_html = expand("results/QC/alignment/qualimap/rnaseq/{sample}/qualimapReport.html", sample = config["sample"]),
        qmap_rnaseq_txt = expand("results/QC/alignment/qualimap/rnaseq/{sample}/rnaseq_qc_results.txt", sample = config["sample"]),
        qmap_multi = "results/QC/alignment/qualimap/multi_bamqc/multisampleBamQcReport.html", 
        samtools_stats = expand("results/QC/alignment/samtools_stats/{sample}_Aligned.sortedByCoord.out.bam.stats", sample = config["sample"]),
        samtools_flagstat = expand("results/QC/alignment/samtools_stats/{sample}_Aligned.sortedByCoord.out.bam.flagstat",  sample = config["sample"]),
        rseqc_strandiness = expand("results/QC/alignment/rseqc/{sample}_strandiness.txt", sample = config["sample"])
    output:
        multiqc = "results/QC/MultiQC/merged/multiqc_report.html"
    conda: 
        config["conda_envs"]["qc"]
    threads: get_resource(config, "default", "threads")
    resources:
        mem_mb = get_resource(config, "default", "mem_mb"),
        runtime = get_resource(config, "default", "runtime")
    params: 
        outdir = lambda wildcards, output : os.path.dirname(output.multiqc)
    log:
        log = "log/QC/MultiQC/multiqc_report_global.log",
    shell: 
        "multiqc {input} -o {params.outdir} 2> {log.log} "
