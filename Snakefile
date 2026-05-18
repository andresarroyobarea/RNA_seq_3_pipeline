configfile: "config/config.yaml"

import glob
import os
import sys
from workflow.utils.utils import get_resource

# ---- Config and global variables ---- #
# Load samples
units = pd.read_table(config["units"], dtype=str).set_index(["sample", "lane"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

# Samples
samples=units['sample'].unique()

# Lanes
lanes=units['lane'].unique()

## Let read if there are UMIs in the sequences.
UMIs = config["umi_processing"]["enabled"]

## Let stablish where the trimming input will be saved.
if UMIs:
    dir_in_trim = "results/umi_extract"
else:
    dir_in_trim = "data"

## Let stablish the path for BAM files in each case
if UMIs:
    aligned_reads = "results/alignment/dedup/{sample}_Aligned.sortedByCoord.out.bam"
else:
    aligned_reads = "results/alignment/{sample}_Aligned.sortedByCoord.out.bam"

# ----------- RULE MODULES -----------
include: "workflow/rules/preprocess.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/concat_fastq.smk"
include: "workflow/rules/trimming.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/quantification.smk"

rule all:
    input:
        "results/feature_counts/counts.tsv",
        "results/QC/alignment/qualimap/multi_bamqc/multisampleBamQcReport.html",
        "results/QC/MultiQC/raw/multiqc_report.html",
        "results/QC/MultiQC/trimmed/multiqc_report.html",
        "results/QC/MultiQC/merged/multiqc_report.html"
