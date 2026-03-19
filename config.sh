#!/bin/bash
#  config.sh — Edit these paths before running the pipeline

#  Reference files (only things that truly need absolute paths) 
GENOMEFASTA="/path/to/genome.fa"
GENOMEGTF="/path/to/annotation.gtf"
HISATINDEX="/path/to/hisat_index/prefix"   # prefix, not a directory

#  Resources 
THREADS=16

#  Conda environment names 
ENV_HISAT="hisat2_env"
ENV_SAMTOOLS="samtools_env"          # rename if you have a dedicated samtools env
ENV_FEATURECOUNTS="featurecounts"
ENV_R="r_env"                    # env that has R + DESeq2 + ggplot2

#  DEG analysis 
# Path to a two-column TSV: sample_name <TAB> condition
# sample_name must match the BAM basename (without _sorted.bam)
# Example:
#   Sample_A1   control
#   Sample_A2   control
#   Sample_B1   treatment
METADATA="/path/to/metadata.tsv"

# Reference condition for DESeq2 contrast (must match a value in METADATA)
REF_CONDITION="control"

# Adjusted p-value and log2FC cutoffs for DEG calling
PADJ_CUTOFF=0.05
LFC_CUTOFF=1
