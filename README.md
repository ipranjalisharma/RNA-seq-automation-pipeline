# RNA-seq Analysis Pipeline

A generic, end-to-end paired-end RNA-seq pipeline covering raw QC through to differential expression analysis. Drop your raw reads into the project folder, fill in `config.sh`, and run.

---

## Table of Contents

- [Overview](#overview)
- [Directory Structure](#directory-structure)
- [Dependencies](#dependencies)
- [Setup](#setup)
- [Running the Pipeline](#running-the-pipeline)
- [Pipeline Steps](#pipeline-steps)
- [Outputs](#outputs)
- [Metadata File Format](#metadata-file-format)
- [Customisation](#customisation)
- [Troubleshooting](#troubleshooting)

---

## Overview

```
Raw FASTQ → QC → Trimming → Re-QC → Alignment → BAM processing → Counting → Normalisation → DEG Analysis
```

| Step | Tool | Output |
|------|------|--------|
| Raw QC | FastQC + MultiQC | HTML reports |
| Trimming | fastp | Trimmed FASTQ |
| Trimmed QC | FastQC + MultiQC | HTML reports |
| Alignment | HISAT2 | SAM files |
| BAM processing | samtools | Sorted BAM + flagstat |
| Aligned QC | FastQC + MultiQC | HTML reports |
| Counting | featureCounts | Raw count matrix |
| Normalisation | DESeq2 (R) | Normalised count tables + plots |
| DEG analysis | DESeq2 (R) | Results tables + volcano plots |

---

## Directory Structure

Place the three pipeline files in a folder together. Your project folder should look like this before running:

```
project/
├── config.sh                  # ← You edit this
├── rnaseq_pipeline.sh         # Main pipeline script
├── deseq2_analysis.R          # Normalisation + DEG script
├── metadata.tsv               # Sample → condition mapping (you create this)
└── raws/                      # ← Put your raw FASTQ files here
    ├── Sample_A1_1.fastq.gz
    ├── Sample_A1_2.fastq.gz
    ├── Sample_B1_1.fastq.gz
    └── Sample_B1_2.fastq.gz
```
NOTE: WORKDIR is automatically set to the project directory when running the pipeline. All output folders are created automatically during the run.

---

## Dependencies

### Tools (must be in PATH or in conda environments)

| Tool | Version tested | Conda env (default name) |
|------|---------------|--------------------------|
| FastQC | ≥ 0.12 | base or any |
| MultiQC | ≥ 1.19 | base or any |
| fastp | ≥ 0.23 | base or any |
| HISAT2 | ≥ 2.2 | `hisat2_env` |
| samtools | ≥ 1.18 | `samtools_env` |
| featureCounts (Subread) | ≥ 2.0 | `featurecounts` |
| R | ≥ 4.3 | `r_env` |

> Conda environment names are just defaults. Change them in `config.sh` to match your setup.

##System Requirements

RAM: 16 GB minimum (32 GB recommended for human genome)

Storage: Allow 2-3× raw data size for intermediate files

CPU: Multi-core recommended (set THREADS appropriately)


### R packages

Install once inside your R environment:

```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "ggplot2", "ggrepel",
                       "dplyr", "tidyr", "pheatmap",
                       "RColorBrewer", "ashr"))
```

### Reference files required

- Genome FASTA (e.g. `Homo_sapiens.GRCh38.dna.toplevel.fa`)
- Genome annotation GTF (e.g. `Homo_sapiens.GRCh38.115.gtf`)
- Pre-built HISAT2 index (see [Building a HISAT2 index](#building-a-hisat2-index) below)

---

## Setup

### 1. Edit config.sh

Open `config.sh` and fill in the paths and settings for your project:

```bash
# Absolute paths to reference files
GENOMEFASTA="/path/to/genome.fa"
GENOMEGTF="/path/to/annotation.gtf"
HISATINDEX="/path/to/hisat_index/human_genome"  # index prefix, not directory

# CPU threads to use
THREADS=16

# Conda environment names (match what you have installed)
ENV_HISAT="hisat2_env"
ENV_SAMTOOLS="star_env"
ENV_FEATURECOUNTS="featurecounts"
ENV_R="r_env"

# Metadata file path and DEG thresholds
METADATA="/path/to/metadata.tsv"
REF_CONDITION="control"
PADJ_CUTOFF=0.05
LFC_CUTOFF=1
```

### 2. Prepare the metadata file

Create a tab-separated file with no header — two columns: sample name and condition. The sample name must exactly match the FASTQ basename (without `_1.fastq.gz`).

```
Sample_A1    control
Sample_A2    control
Sample_A3    control
Sample_B1    treatment
Sample_B2    treatment
Sample_B3    treatment
```
Can create using this bash command:
``` bash
printf "Sample_A1\tcontrol\nSample_B1\ttreatment\n" > metadata.tsv
```
Save it anywhere and point `METADATA` in `config.sh` to it.

### 3. Name your FASTQ files

Raw files must follow the `{sample}_1.fastq.gz` / `{sample}_2.fastq.gz` naming convention for paired-end reads. If your files are named differently, rename them before running:

```bash
# Example rename
mv SRR12345_R1.fastq.gz Sample_A1_1.fastq.gz
mv SRR12345_R2.fastq.gz Sample_A1_2.fastq.gz
```

### Building a HISAT2 index

If you don't already have a HISAT2 index, build one like this:

```bash
conda activate hisat2_env
mkdir -p hisat_index
hisat2-build -p 16 /path/to/genome.fa hisat_index/human_genome
```

Point `HISATINDEX` in `config.sh` to `hisat_index/human_genome` (the prefix, not the folder).

---

## Running the Pipeline

```bash
# From inside the project directory
bash rnaseq_pipeline.sh

# Or point to the project directory explicitly
bash rnaseq_pipeline.sh /path/to/project
```

The pipeline will print progress for each step and write per-tool logs to `logs/`. If any step fails, the pipeline stops immediately (`set -euo pipefail`) rather than silently continuing with bad data.

### Running individual steps

If a step fails and you want to re-run from a specific point, you can comment out the earlier steps in `rnaseq_pipeline.sh` (each step is clearly delimited by a `STEP N` header). Alternatively, run just the R analysis directly:

```bash
source config.sh
export COUNT_FILE="$WORKDIR/counts/raw_counts.txt"
export METADATA_FILE="$METADATA"
export NORMAL_DIR="$WORKDIR/normal"
export DEG_DIR="$WORKDIR/deg"
export REF_CONDITION PADJ_CUTOFF LFC_CUTOFF
Rscript deseq2_analysis.R
```

---

## Pipeline Steps

### Step 1 — Raw QC
FastQC is run on every raw FASTQ file, then MultiQC aggregates all reports into `qcraw/raw_multiqc_report.html`. Check this before proceeding to confirm read quality.

### Step 2 — Trimming
fastp trims adapter sequences, removes low-quality bases (Phred < 20), discards short reads (< 36 bp), and trims polyG tails. Both paired reads are processed together. Per-sample HTML reports are written to `trimming/`.

### Step 3 — Trimmed QC
FastQC and MultiQC are repeated on the trimmed reads for confirmation.

### Step 4 — Alignment
HISAT2 aligns trimmed paired-end reads to the reference genome. Alignment summary logs (mapping rate, concordant pairs, etc.) are written per sample to `logs/`.

### Step 5 — BAM processing
samtools converts each SAM to a coordinate-sorted BAM, indexes it, and writes a flagstat report. SAM files are deleted after successful conversion to save disk space.

### Step 6 — Aligned QC
FastQC and MultiQC run on the sorted BAMs for a final alignment quality check.

### Step 7 — Feature counting
featureCounts quantifies reads over gene features from the GTF annotation. All samples are counted together into a single matrix (`counts/raw_counts.txt`).

### Step 8 — Normalisation (R / DESeq2)
Three normalised count tables are produced:

| File | Method | Best used for |
|------|--------|---------------|
| `deseq2_normalised_counts.csv` | DESeq2 size factors | Between-sample comparison |
| `vst_counts.csv` | Variance stabilising transform | Heatmaps, PCA, clustering |
| `tpm_counts.csv` | TPM | Comparing expression levels within a sample |

### Step 9 — DEG analysis (R / DESeq2)
DESeq2 is run for every pairwise contrast among conditions. Results use `lfcShrink` with the `ashr` method for shrunken, unbiased fold-change estimates. Genes are called significant at `padj < PADJ_CUTOFF` and `|log2FC| > LFC_CUTOFF`.

---

## Outputs

```
project/
├── logs/                          # Per-tool log files
├── qcraw/
│   └── raw_multiqc_report.html
├── trimming/
│   ├── *_1.fq.gz / *_2.fq.gz     # Trimmed reads
│   └── *_fastp.html               # Per-sample fastp report
├── qctrimmed/
│   └── trimmed_multiqc_report.html
├── alignment/
│   ├── *_sorted.bam
│   ├── *_sorted.bam.bai
│   └── *_flagstat.txt
├── qcaligned/
│   └── aligned_multiqc_report.html
├── counts/
│   └── raw_counts.txt
├── normal/
│   ├── deseq2_normalised_counts.csv
│   ├── vst_counts.csv
│   ├── tpm_counts.csv
│   ├── PCA_plot.pdf / .png
│   ├── sample_distance_heatmap.pdf
│   └── top50_variable_genes_heatmap.pdf
└── deg/
    ├── {contrast}_all_results.csv
    ├── {contrast}_significant_DEGs.csv
    ├── {contrast}_volcano.pdf / .png
    └── ...  (one set per contrast)
```

### DEG results columns

| Column | Description |
|--------|-------------|
| `gene_id` | Ensembl gene ID |
| `baseMean` | Mean normalised count across all samples |
| `log2FoldChange` | Shrunken log2 fold change (ashr) |
| `lfcSE` | Standard error of the fold change |
| `pvalue` | Wald test p-value |
| `padj` | Benjamini-Hochberg adjusted p-value |
| `regulation` | `Up`, `Down`, or `Not significant` |

---

## Metadata File Format

- Tab-separated, no header row
- Column 1: sample name (must match FASTQ basename exactly)
- Column 2: condition label
- Minimum 2 samples per condition recommended; 3+ required for robust DEG calling

```
Sample_ctrl_1	control
Sample_ctrl_2	control
Sample_ctrl_3	control
Sample_treat_1	treatment
Sample_treat_2	treatment
Sample_treat_3	treatment
```

For more than two conditions, all pairwise contrasts are run automatically.

---

## Customisation

**Change DEG thresholds** — edit `PADJ_CUTOFF` and `LFC_CUTOFF` in `config.sh` and re-run the R script only.

**Change fastp parameters** — edit the `fastp` block in `rnaseq_pipeline.sh`. Common additions:

```bash
--cut_front              # trim low-quality bases from 5' end
--cut_tail               # trim low-quality bases from 3' end
--n_base_limit 5         # discard reads with > 5 N bases
```

**Add a gene name column to DEG results** — if you have a gene ID → symbol mapping table, add this to `deseq2_analysis.R` after `res_df` is created:

```r
symbols <- read.table("gene_id_to_symbol.tsv", header=FALSE,
                      col.names=c("gene_id","gene_name"))
res_df <- left_join(res_df, symbols, by="gene_id")
```

**Stranded libraries** — add `-s 1` (forward) or `-s 2` (reverse) to the `featureCounts` call in `rnaseq_pipeline.sh`.

---

## Troubleshooting

**Pipeline stops with "conda activate" error**
Ensure `conda init bash` has been run and your shell has been restarted, or that `$(conda info --base)/etc/profile.d/conda.sh` exists.

**featureCounts reports very low assignment rate**
Check that the GTF chromosome naming (e.g. `1` vs `chr1`) matches the genome FASTA. They must be consistent.

**DESeq2 error: "all samples have zero counts for gene X"**
This is normal — the low-count filter (`rowSums >= 10`) should catch most of these. If the error persists, check that the featureCounts column names were stripped correctly (look at the top of the `deseq2.log` file).

**Volcano plot has no labelled genes**
Either no genes pass the thresholds, or the `gene_id` column contains only Ensembl IDs with no symbol. Try relaxing `PADJ_CUTOFF` or `LFC_CUTOFF` in `config.sh` to confirm the analysis ran correctly.

**Out of memory during HISAT2 alignment**
Reduce `THREADS` in `config.sh`, or request more RAM. HISAT2 typically requires 8–16 GB for a human genome index.

---

## Citation

If you use this pipeline, please cite the underlying tools:

- **FastQC** — Andrews, S., 2017. FastQC: a quality control tool for high throughput sequence data. 2010 [online]
- **MultiQC** — Ewels, P., Magnusson, M., Lundin, S. and Käller, M., 2016. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), pp.3047-3048.
- **fastp** — Chen, S., Zhou, Y., Chen, Y. and Gu, J., 2018. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), pp.i884-i890.
- **HISAT2** — Kim, D., Paggi, J.M., Park, C., Bennett, C. and Salzberg, S.L., 2019. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nature biotechnology, 37(8), pp.907-915.
- **samtools** — Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R. and 1000 Genome Project Data Processing Subgroup, 2009. The sequence alignment/map format and SAMtools. bioinformatics, 25(16), pp.2078-2079.
- **featureCounts** — Liao, Y., Smyth, G.K. and Shi, W., 2014. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7), pp.923-930.
- **DESeq2** — Love, M.I., Huber, W. and Anders, S., 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome biology, 15(12), p.550.
