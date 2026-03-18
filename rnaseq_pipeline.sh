#!/bin/bash
#  rnaseq_pipeline.sh
#  Generic paired-end RNA-seq pipeline
#  Usage:  bash rnaseq_pipeline.sh [/optional/working/dir]
set -euo pipefail

#Working directory: use argument or current directory 
WORKDIR="${1:-$(pwd)}"
cd "$WORKDIR"
echo "Working directory: $WORKDIR"

#  Load config 
CONFIG="$WORKDIR/config.sh"
if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: config.sh not found in $WORKDIR. Copy and edit config.sh first."
    exit 1
fi
source "$CONFIG"

# Validate required reference files
for f in "$GENOMEFASTA" "$GENOMEGTF" "$METADATA"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Required file not found: $f"
        exit 1
    fi
done

# Directory layout
RAW="$WORKDIR/raws"
TRIM="$WORKDIR/trimming"
ALIGN="$WORKDIR/alignment"
COUNT="$WORKDIR/counts"
NORMAL="$WORKDIR/normal"
DEG="$WORKDIR/deg"
LOGS="$WORKDIR/logs"

mkdir -p "$TRIM" "$ALIGN" "$COUNT" "$NORMAL" "$DEG" "$LOGS" \
         "$WORKDIR/qcraw" "$WORKDIR/qctrimmed" "$WORKDIR/qcaligned"

if [[ ! -d "$RAW" ]]; then
    echo "ERROR: Raw reads directory not found: $RAW"
    exit 1
fi

#Helper: activate conda env
activate_env() {
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$1"
    echo "Activated conda environment: $1"
}

#Step 1: QC on raw reads 
echo "STEP 1: QC on raw reads"

ls "$RAW"/*.fastq.gz | while read -r i; do
    echo "Running FastQC: $(basename "$i")"
    fastqc "$i" --outdir "$WORKDIR/qcraw" -t "$THREADS" \
        2>> "$LOGS/fastqc_raw.log"
done
multiqc "$WORKDIR/qcraw" --outdir "$WORKDIR/qcraw" \
    -n raw_multiqc_report --force 2>> "$LOGS/multiqc_raw.log"
echo "Raw QC done."

# Step 2: Trimming with fastp 

echo "STEP 2: Trimming with fastp"

ls "$RAW"/*_1.fastq.gz | while read -r i; do
    base=$(basename "$i" _1.fastq.gz)
    j="${i/_1.fastq.gz/_2.fastq.gz}"

    if [[ ! -f "$j" ]]; then
        echo "WARNING: Paired file not found for $base, skipping."
        continue
    fi

    echo "Trimming: $base"
    fastp \
        -i "$i"  -I "$j" \
        -o "$TRIM/${base}_1.fq.gz" \
        -O "$TRIM/${base}_2.fq.gz" \
        -h "$TRIM/${base}_fastp.html" \
        -j "$TRIM/${base}_fastp.json" \
        --thread "$THREADS" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --trim_poly_g \
        --length_required 36 \
        2>> "$LOGS/fastp_${base}.log"
    echo "Trimming done: $base"
done
echo "All trimming done."

echo "Running FastQC on trimmed reads..."
ls "$TRIM"/*.fq.gz | while read -r i; do
    fastqc "$i" --outdir "$WORKDIR/qctrimmed" -t "$THREADS" \
        2>> "$LOGS/fastqc_trimmed.log"
done
multiqc "$WORKDIR/qctrimmed" --outdir "$WORKDIR/qctrimmed" \
    -n trimmed_multiqc_report --force 2>> "$LOGS/multiqc_trimmed.log"

# Step 3: Alignment with hisat2
echo "STEP 3: Alignment with HISAT2"

activate_env "$ENV_HISAT"

ls "$TRIM"/*_1.fq.gz | while read -r i; do
    j="${i/_1.fq.gz/_2.fq.gz}"
    base=$(basename "$i" _1.fq.gz)

    echo "Aligning: $base"
    hisat2 -p "$THREADS" \
        -x "$HISATINDEX" \
        -1 "$i" -2 "$j" \
        -S "$ALIGN/${base}.sam" \
        --summary-file "$LOGS/hisat2_${base}.log" \
        2>> "$LOGS/hisat2_${base}.log"
    echo "Alignment done: $base"
done
echo "All alignments done."

#  Step 4: SAM -> sorted BAM -> index -> flagstat

echo "STEP 4: SAM → sorted BAM + QC"

activate_env "$ENV_SAMTOOLS"

ls "$ALIGN"/*.sam | while read -r i; do
    base=$(basename "$i" .sam)
    echo "Processing: $base"

    samtools sort -@ "$THREADS" -o "$ALIGN/${base}_sorted.bam" "$i" \
        2>> "$LOGS/samtools_${base}.log"
    samtools index "$ALIGN/${base}_sorted.bam" \
        2>> "$LOGS/samtools_${base}.log"
    samtools flagstat "$ALIGN/${base}_sorted.bam" \
        > "$ALIGN/${base}_flagstat.txt"

    echo "Done: $base"
done

echo "Removing SAM files..."
rm "$ALIGN"/*.sam

echo "Running FastQC + MultiQC on BAM files..."
ls "$ALIGN"/*_sorted.bam | while read -r i; do
    fastqc "$i" --outdir "$WORKDIR/qcaligned" -t "$THREADS" \
        2>> "$LOGS/fastqc_aligned.log"
done
multiqc "$WORKDIR/qcaligned" --outdir "$WORKDIR/qcaligned" \
    -n aligned_multiqc_report --force 2>> "$LOGS/multiqc_aligned.log"
echo "BAM QC done."

# Step 5: Feature counts 

echo "STEP 5: featureCounts"

activate_env "$ENV_FEATURECOUNTS"

featureCounts \
    -p -T "$THREADS" \
    -a "$GENOMEGTF" \
    -o "$COUNT/raw_counts.txt" \
    "$ALIGN"/*_sorted.bam \
    2>> "$LOGS/featurecounts.log"

echo "featureCounts done."

#  Step 6: Normalization + DEG with DESeq2 (R) 

echo "STEP 6: Normalization and DEG analysis (DESeq2)"
activate_env "$ENV_R"

# Pass all config variables as env vars into the R script
export COUNT_FILE="$COUNT/raw_counts.txt"
export METADATA_FILE="$METADATA"
export NORMAL_DIR="$NORMAL"
export DEG_DIR="$DEG"
export REF_CONDITION
export PADJ_CUTOFF
export LFC_CUTOFF

RSCRIPT_PATH="$WORKDIR/deseq2_analysis.R"
if [[ ! -f "$RSCRIPT_PATH" ]]; then
    echo "ERROR: deseq2_analysis.R not found in $WORKDIR"
    exit 1
fi

Rscript "$RSCRIPT_PATH" 2>&1 | tee "$LOGS/deseq2.log"
echo "Normalization and DEG analysis done."

#  Done

echo "Pipeline complete!"
echo "  Raw counts    : $COUNT/raw_counts.txt"
echo "  Norm counts   : $NORMAL/"
echo "  DEG results   : $DEG/"
echo "  QC reports    : qcraw / qctrimmed / qcaligned"
echo "  Logs          : $LOGS/"

