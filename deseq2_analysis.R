#  deseq2_analysis.R
#  Normalization (DESeq2 size factors + TPM) and DEG analysis
#  All inputs come from environment variables set by the pipeline

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
    library(dplyr)
    library(tidyr)
    library(pheatmap)
    library(RColorBrewer)
})

# Read config from environment variables
count_file    <- Sys.getenv("COUNT_FILE")
metadata_file <- Sys.getenv("METADATA_FILE")
normal_dir    <- Sys.getenv("NORMAL_DIR")
deg_dir       <- Sys.getenv("DEG_DIR")
ref_condition <- Sys.getenv("REF_CONDITION")
padj_cutoff   <- as.numeric(Sys.getenv("PADJ_CUTOFF", unset = "0.05"))
lfc_cutoff    <- as.numeric(Sys.getenv("LFC_CUTOFF",  unset = "1"))

cat("DESeq2 Analysis \n")
cat("Count file     :", count_file, "\n")
cat("Metadata file  :", metadata_file, "\n")
cat("Reference cond :", ref_condition, "\n")
cat("padj cutoff    :", padj_cutoff, "\n")
cat("LFC  cutoff    :", lfc_cutoff, "\n\n")

dir.create(normal_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(deg_dir,    showWarnings = FALSE, recursive = TRUE)

#  1. Load raw counts 
# featureCounts output: skip the comment line, columns 1 + 7 onward
raw <- read.table(count_file, header = TRUE, sep = "\t",
                  skip = 1, row.names = 1, check.names = FALSE)

# Keep only count columns (cols 6 onward after Geneid; drop Chr,Start,End,Strand,Length)
count_data <- raw[, 6:ncol(raw)]

# Clean column names: strip path + _sorted.bam suffix
colnames(count_data) <- sub(".*\\/", "",          colnames(count_data))  # remove path
colnames(count_data) <- sub("_sorted\\.bam$", "", colnames(count_data))  # remove suffix

cat("Samples detected:\n")
print(colnames(count_data))

# 2. Load metadata
meta <- read.table(metadata_file, header = FALSE, sep = "\t",
                   col.names = c("sample", "condition"))
rownames(meta) <- meta$sample

# Ensure sample order matches count matrix
if (!all(colnames(count_data) %in% rownames(meta))) {
    missing <- setdiff(colnames(count_data), rownames(meta))
    stop("Samples in count matrix not found in metadata: ",
         paste(missing, collapse = ", "))
}
meta <- meta[colnames(count_data), , drop = FALSE]
meta$condition <- relevel(factor(meta$condition), ref = ref_condition)

cat("\nSample metadata:\n")
print(meta)

# 3. Build DESeq2 object 
dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData   = meta,
    design    = ~ condition
)

# Pre-filter: keep genes with at least 10 counts in total
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
cat("\nGenes after low-count filter:", nrow(dds), "\n")

#  4. Run DESeq2 
cat("\nRunning DESeq2...\n")
dds <- DESeq(dds)

#  5. Normalization outputs
cat("\nExporting normalised counts...\n")

## 5a. DESeq2 size-factor normalised counts
norm_counts <- counts(dds, normalized = TRUE)
write.csv(
    as.data.frame(norm_counts),
    file = file.path(normal_dir, "deseq2_normalised_counts.csv")
)

## 5b. VST (variance-stabilising transformation) — best for visualisation
vst_mat <- assay(vst(dds, blind = FALSE))
write.csv(
    as.data.frame(vst_mat),
    file = file.path(normal_dir, "vst_counts.csv")
)

## 5c. TPM — needs gene lengths from featureCounts column 5
gene_lengths <- raw[rownames(dds), 5]   # "Length" column

rpk        <- norm_counts / (gene_lengths / 1e3)
tpm        <- sweep(rpk, 2, colSums(rpk) / 1e6, "/")
write.csv(
    as.data.frame(tpm),
    file = file.path(normal_dir, "tpm_counts.csv")
)
cat("Saved: deseq2_normalised_counts.csv | vst_counts.csv | tpm_counts.csv\n")

#  6. DEG results 
conditions <- levels(meta$condition)
contrasts_list <- combn(conditions, 2, simplify = FALSE)

all_degs <- list()

for (pair in contrasts_list) {
    cond_test <- pair[pair != ref_condition]
    cond_ref  <- ref_condition
    # If neither is the reference, just use as-is
    if (length(cond_test) == 0) {
        cond_test <- pair[1]; cond_ref <- pair[2]
    }

    contrast_name <- paste0(cond_test, "_vs_", cond_ref)
    cat("\nContrast:", contrast_name, "\n")

    res <- results(dds,
                   contrast  = c("condition", cond_test, cond_ref),
                   alpha     = padj_cutoff)
    res <- lfcShrink(dds,
                     contrast  = c("condition", cond_test, cond_ref),
                     res       = res,
                     type      = "ashr")

    res_df <- as.data.frame(res) %>%
        tibble::rownames_to_column("gene_id") %>%
        arrange(padj)

    # Classify
    res_df$regulation <- "Not significant"
    res_df$regulation[res_df$padj < padj_cutoff & res_df$log2FoldChange >  lfc_cutoff] <- "Up"
    res_df$regulation[res_df$padj < padj_cutoff & res_df$log2FoldChange < -lfc_cutoff] <- "Down"

    sig_df <- res_df %>% filter(regulation != "Not significant")

    cat("  Up  :", sum(res_df$regulation == "Up"), "\n")
    cat("  Down:", sum(res_df$regulation == "Down"), "\n")

    # Save full + significant results
    write.csv(res_df,
              file = file.path(deg_dir, paste0(contrast_name, "_all_results.csv")),
              row.names = FALSE)
    write.csv(sig_df,
              file = file.path(deg_dir, paste0(contrast_name, "_significant_DEGs.csv")),
              row.names = FALSE)

    all_degs[[contrast_name]] <- res_df

    #  Volcano plot 
    top_genes <- res_df %>%
        filter(regulation != "Not significant") %>%
        arrange(padj) %>%
        slice_head(n = 20)

    p_volcano <- ggplot(res_df,
                        aes(x = log2FoldChange, y = -log10(padj),
                            colour = regulation)) +
        geom_point(alpha = 0.6, size = 1.5) +
        geom_text_repel(data = top_genes,
                        aes(label = gene_id),
                        size = 3, max.overlaps = 20) +
        scale_colour_manual(values = c("Up" = "#E64B35",
                                       "Down" = "#4DBBD5",
                                       "Not significant" = "grey70")) +
        geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
                   linetype = "dashed", colour = "grey40") +
        geom_hline(yintercept = -log10(padj_cutoff),
                   linetype = "dashed", colour = "grey40") +
        labs(title   = paste("Volcano:", contrast_name),
             x       = expression(log[2]~"Fold Change"),
             y       = expression(-log[10]~"(adjusted p-value)"),
             colour  = "Regulation") +
        theme_bw(base_size = 13) +
        theme(plot.title = element_text(hjust = 0.5))

    ggsave(file.path(deg_dir, paste0(contrast_name, "_volcano.pdf")),
           p_volcano, width = 8, height = 6)
    ggsave(file.path(deg_dir, paste0(contrast_name, "_volcano.png")),
           p_volcano, width = 8, height = 6, dpi = 300)
}

#  7. Exploratory plots (all samples) 
cat("\nGenerating exploratory plots...\n")

## PCA
pca_data <- plotPCA(vst(dds, blind = TRUE), intgroup = "condition",
                    returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, colour = condition, label = name)) +
    geom_point(size = 4) +
    geom_text_repel(size = 3) +
    xlab(paste0("PC1 (", pct_var[1], "% variance)")) +
    ylab(paste0("PC2 (", pct_var[2], "% variance)")) +
    labs(title = "PCA — VST counts") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(normal_dir, "PCA_plot.pdf"), p_pca, width = 7, height = 5)
ggsave(file.path(normal_dir, "PCA_plot.png"), p_pca, width = 7, height = 5, dpi = 300)

## Sample distance heatmap
samp_dists  <- dist(t(vst_mat))
samp_mat    <- as.matrix(samp_dists)
anno_col    <- data.frame(condition = meta$condition, row.names = rownames(meta))
colours     <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pdf(file.path(normal_dir, "sample_distance_heatmap.pdf"), width = 7, height = 6)
pheatmap(samp_mat,
         clustering_distance_rows = samp_dists,
         clustering_distance_cols = samp_dists,
         annotation_col           = anno_col,
         col                      = colours,
         main                     = "Sample Distance Heatmap")
dev.off()

## Top variable genes heatmap (top 50 by variance across samples)
top_var_genes <- head(order(rowVars(vst_mat), decreasing = TRUE), 50)
heatmap_mat   <- vst_mat[top_var_genes, ]
heatmap_mat   <- heatmap_mat - rowMeans(heatmap_mat)  # center

pdf(file.path(normal_dir, "top50_variable_genes_heatmap.pdf"), width = 8, height = 10)
pheatmap(heatmap_mat,
         annotation_col = anno_col,
         show_rownames  = TRUE,
         fontsize_row   = 7,
         main           = "Top 50 Variable Genes (VST, centred)")
dev.off()

cat("\n=== Analysis complete ===\n")
cat("Normalised counts :", normal_dir, "\n")
cat("DEG results       :", deg_dir,    "\n")
