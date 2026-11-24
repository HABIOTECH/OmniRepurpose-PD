#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
})

base_dir <- "path/to/project_root"
pb_dir <- file.path(base_dir, "results/pseudobulk")
out_dir <- file.path(base_dir, "results/de_deseq2_r")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# All cell types present in pseudobulk
all_cell_types <- c(
  "Microglia", "Astrocytes", "Oligodendrocytes", "OPCs", 
  "Excitatory", "Inhibitory", "Endothelial_cells", "Pericytes", 
  "Ependymal", "GABA", "DaNs"
)

summary_list <- list()

for (ct in all_cell_types) {
  cat(sprintf("\n=== Running R DESeq2 for: %s ===\n", ct))
  counts_file <- file.path(pb_dir, sprintf("counts_%s.tsv", ct))
  coldata_file <- file.path(pb_dir, sprintf("coldata_%s.tsv", ct))
  
  if (!file.exists(counts_file) || !file.exists(coldata_file)) {
    cat(sprintf("Skipping %s (files not found)\n", ct))
    next
  }
  
  counts <- as.matrix(read.table(counts_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE))
  coldata <- read.table(coldata_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
  
  coldata <- coldata[colnames(counts), ]
  coldata$diagnosis <- factor(coldata$diagnosis, levels = c("Control", "PD"))
  coldata$sex <- factor(coldata$sex)
  coldata$age <- as.numeric(coldata$age)
  # Center age to improve GLM convergence
  coldata$age <- coldata$age - mean(coldata$age)
  
  # Pre-filter low count genes: at least 10 counts across donors
  keep <- rowSums(counts) >= 10
  counts <- counts[keep, , drop=FALSE]
  
  if (nrow(counts) < 50) {
    cat(sprintf("Skipping %s: too few expressed genes (%d)\n", ct, nrow(counts)))
    next
  }
  
  # For DaNs (only 4 samples), use ~ diagnosis to avoid rank deficiency
  if (nrow(coldata) <= 4) {
    cat(sprintf("Note: %s has only %d donors. Using design = ~ diagnosis\n", ct, nrow(coldata)))
    design_formula <- ~ diagnosis
  } else {
    design_formula <- ~ sex + age + diagnosis
  }
  
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = coldata,
                                design = design_formula)
  
  dds <- DESeq(dds, quiet=TRUE)
  res <- results(dds, contrast=c("diagnosis", "PD", "Control"))
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$cell_type <- ct
  
  out_csv <- file.path(out_dir, sprintf("deseq2_R_%s.tsv", ct))
  write.table(res_df, file=out_csv, sep="\t", quote=FALSE, row.names=FALSE)
  
  sig_nom <- sum(!is.na(res_df$pvalue) & res_df$pvalue < 0.05)
  sig_fdr <- sum(!is.na(res_df$padj) & res_df$padj < 0.1)
  
  cat(sprintf("Completed %s: tested %d genes | %d nominal (p<0.05) | %d FDR (padj<0.1)\n",
              ct, nrow(res_df), sig_nom, sig_fdr))
}
