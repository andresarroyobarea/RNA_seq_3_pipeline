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
    threads: 
        get_resource(config, "fastqc", "threads")
    resources:
        mem_mb = get_resource(config, "fastqc", "mem_mb"),
        runtime = get_resource(config, "fastqc", "runtime")
    params: 
        outdir = lambda wildcards, output : os.path.dirname(output.html),
        extra = config["parameters"]["fastqc"]["extra"]
    log:
        "log/QC/alignment/fastqc/{sample}.log"
    benchmark: 
        "benchmarks/QC/alignment/fastqc/{sample}.bmk"
    shell: """
        mkdir -p {params.outdir} &&
        fastqc --outdir {params.outdir} \
            --threads {threads} \
            {input.bam} \
            {params.extra} 2> {log} 
    """

rule qualimap_bamqc:
    input:
        bam = "results/alignment/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        qmap_report = "results/QC/alignment/qualimap/bamqc/{sample}/qualimapReport.html",
        genome_res = "results/QC/alignment/qualimap/bamqc/{sample}/genome_results.txt"
    conda:
        config["conda_envs"]["qualimap"]
    threads: 
        get_resource(config, "qualimap", "threads")
    resources:
        mem_mb = get_resource(config, "qualimap", "mem_mb"),
        runtime = get_resource(config, "qualimap", "runtime")
    params:
        genome = config["qualimap"]["genome"],
        annotation = config["genome"]["annotation_gtf"],
        outdir = lambda wildcards, output: os.path.dirname(output.qmap_report),
        mem = f"{get_resource(config, 'qualimap', 'mem_mb') // 1024}G"
        extra_single = config["parameters"]["qualimap"]["extra"]
    log: 
        "log/QC/alignment/qualimap/bamqc/{sample}.log"
    benchmark:
        "benchmarks/QC/alignment/qualimap/bamqc/{sample}.bmk"
    shell: """
        qualimap bamqc \
            -bam {input.bam} \
            -gd {params.genome} \ 
            -gff {params.annotation} \
            -nt {threads} \
            --outdir {params.outdir} \
            {params.extra_single} \
            --java-mem-size={params.mem} 2> {log}
    """

rule qualimap_multi_bamqc:
    input:
        expand("results/QC/alignment/qualimap/bamqc/{sample}/qualimapReport.html", sample = samples)
    output:
        qmap_report = "results/QC/alignment/qualimap/multi_bamqc/multisampleBamQcReport.html"
    conda:
        config["conda_envs"]["qualimap"]
    threads: 
        get_resource(config, "qualimap", "threads")
    resources:
        mem_mb = get_resource(config, "qualimap", "mem_mb"),
        runtime = get_resource(config, "qualimap", "runtime")
    params: 
        qmap_input = "metadata/qualimap_multi_bamqc_input.txt",
        outdir = lambda wildcards, output: os.path.dirname(output.qmap_report),
        extra_multi = config["parameters"]["qualimap"]["extra_multi"]
    log: 
        "log/QC/alignment/qualimap/mutli_bamqc/qualimap_multi_bamqc.log"
    benchmark:
        "benchmarks/QC/alignment/qualimap/multi_bamqc/qualimap_multi_bamqc.bmk"
    shell:"""
        qualimap multi-bamqc \
            -d {params.qmap_input} \
            --outdir {params.outdir} \
            {params.extra_multi} 2> {log} 
    """

rule qualimap_rnaseq:
    input:  
        bam = "results/alignment/{sample}_Aligned.sortedByCoord.out.bam"
    output: 
        qmap_report = "results/QC/alignment/qualimap/rnaseq/{sample}/qualimapReport.html",
        qmap_res = "results/QC/alignment/qualimap/rnaseq/{sample}/rnaseq_qc_results.txt"
    conda:
        config["conda_envs"]["qualimap"]
    threads: 
        get_resource(config, "qualimap", "threads")
    resources:
        mem_mb = get_resource(config, "qualimap", "mem_mb"),
        runtime = get_resource(config, "qualimap", "runtime")
    params:
        annotation = config["genome"]["annotation_gtf"],
        outdir = lambda wildcards, output: os.path.dirname(output.qmap_report),
        mem = f"{get_resource(config, 'qualimap', 'mem_mb') // 1024}G",
        extra_rnaseq = config["parameters"]["qualimap"]["extra_rnaseq"]
    log:
        "log/QC/alignment/qualimap/rnaseq/{sample}.log"
    benchmark:
        "benchmarks/QC/alignment/qualimap/rnaseq/{sample}.bmk"
    shell:"""
        qualimap rnaseq \
            -bam {input.bam} \
            -gtf {params.annotation} \
            -outdir {params.outdir} \ 
            --java-mem-size={params.mem} \
            {params.extra_rnaseq} \
            2> {log}
    """

rule rseqc_strand:
    input: 
        bam = "results/alignment/{sample}_Aligned.sortedByCoord.out.bam"
    output: 
        rseqc_out = "results/QC/alignment/rseqc/{sample}_strandiness.txt"
    conda: 
        config["conda_envs"]["rseqc"]
    threads: 
        get_resource(config, "rseqc_strand", "threads")
    resources:
        mem_mb = get_resource(config, "rseqc_strand", "mem_mb"),
        runtime = get_resource(config, "rseqc_strand", "runtime")
    params:
        annotation = config["genome"]["annotation_bed"],
        extra = config["parameters"]["rseqc_strand"]["extra"]
    log:
        "log/QC/alignment/rseqc/{sample}.log"
    benchmark:
        "benchmarks/QC/alignment/rseqc/{sample}.bmk"
    shell:"""
        infer_experiment.py -i {input.bam} \ 
            -r {params.annotation} \
            {params.extra} > {output.rseqc_out} 2> {log}
    """

rule samtools_qc:
    input:
        bam = "results/alignment/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        samtools_stats = "results/QC/alignment/samtools_qc/{sample}_Aligned.sortedByCoord.out.bam.stats",
        samtools_flagstat = "results/QC/alignment/samtools_qc/{sample}_Aligned.sortedByCoord.out.bam.flagstat"
    conda:
        config["conda_envs"]["aligners"]
    threads: 
        get_resource(config, "default", "threads")
    resources:
        mem_mb = get_resource(config, "default", "mem_mb"),
        runtime = get_resource(config, "default", "runtime")
    params:
        extra = config["parameters"]["samtools_qc"]["extra"]
    log:
        log_stats = "log/QC/alignment/samtools_qc/{sample}.log",
        log_flagstat = "log/QC/alignment/samtools_qc/{sample}.log"
    shell:"""
        samtools stats {input.bam} > {output.samtools_stats} 2> {log.log_stats} &&
        samtools flagstat {input.bam} > {output.samtools_flagstat} 2> {log.log_flagstat}
    """

rule multiqc_trimmed:
    input: 
        seqs_QC_trim = expand("results/QC/trimmed/fastqc/{sample}/{sample}_{seq_lane}_trimmed_fastqc.html", sample = samples, seq_lane = lanes),
        fastq_screen_txt = expand("results/QC/trimmed/fastq_screen/{sample}/{sample}_{seq_lane}_trimmed_screen.txt", sample = samples, seq_lane = lanes),
        fastq_screen_png = expand("results/QC/trimmed/fastq_screen/{sample}/{sample}_{seq_lane}_trimmed_screen.png", sample = samples, seq_lane = lanes)
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
        outdir = lambda wildcards, output : os.path.dirname(output.multiqc),
        extra = config["parameters"]["multiqc"]["extra"]
    log:
        "log/QC/MultiQC/trimmed/multiqc_report.log",
    shell: """
        multiqc {params.inputdir} \
            -o {params.outdir} \
            {params.extra} 2> {log} 
    """

rule multiqc_merge:
    input:
        seq_QC_merged = expand("results/QC/merged/FastQC/{sample}/{sample}_fastqc.html", sample = samples),
        fastq_screen_merged_txt = expand("results/QC/merged/fastq_screen/{sample}/{sample}_screen.txt", sample = samples),
        fastq_screen_merged_png = expand("results/QC/merged/fastq_screen/{sample}/{sample}_screen.png", sample = samples),
        bam_QC = expand("results/QC/alignment/FastQC/{sample}/{sample}_Aligned.sortedByCoord.out_fastqc.html", sample = samples),
        umi_extract = expand("log/umitools/extract/{sample}_{seq_lane}.log", sample = samples, seq_lane = lanes),
        umi_dedup = expand("log/umitools/dedup/{sample}.log", sample = samples),
        qmap_bamqc_html = expand("results/QC/alignment/qualimap/bamqc/{sample}/qualimapReport.html", sample = samples),
        qmap_bamqc_txt = expand("results/QC/alignment/qualimap/bamqc/{sample}/genome_results.txt", sample = samples),
        qmap_rnaseq_html = expand("results/QC/alignment/qualimap/rnaseq/{sample}/qualimapReport.html", sample = samples),
        qmap_rnaseq_txt = expand("results/QC/alignment/qualimap/rnaseq/{sample}/rnaseq_qc_results.txt", sample = samples),
        qmap_multi = "results/QC/alignment/qualimap/multi_bamqc/multisampleBamQcReport.html", 
        samtools_stats = expand("results/QC/alignment/samtools_qc/{sample}_Aligned.sortedByCoord.out.bam.stats", sample = samples),
        samtools_flagstat = expand("results/QC/alignment/samtools_qc/{sample}_Aligned.sortedByCoord.out.bam.flagstat",  sample = samples),
        rseqc_strandiness = expand("results/QC/alignment/rseqc/{sample}_strandiness.txt", sample = samples)
    output:
        multiqc = "results/QC/MultiQC/merged/multiqc_report.html"
    conda: 
        config["conda_envs"]["qc"]
    threads: 
        get_resource(config, "default", "threads")
    resources:
        mem_mb = get_resource(config, "default", "mem_mb"),
        runtime = get_resource(config, "default", "runtime"),
        extra = config["parameters"]["multiqc"]["extra"]
    params: 
        outdir = lambda wildcards, output : os.path.dirname(output.multiqc)
    log:
        "log/QC/MultiQC/multiqc_report_global.log",
    shell: """
        multiqc {input} \ 
            -o {params.outdir} \ 
            {params.extra} 2> {log} 
    """