#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
})

base_dir <- "/sc/kzd/home/abdullah/bioinformatics/PK-multiomics"
pb_dir <- file.path(base_dir, "results/pseudobulk")
out_dir <- file.path(base_dir, "results/de_deseq2_r")

cell_types <- c("Microglia", "Astrocytes", "Oligodendrocytes", "OPCs", "Excitatory", "Inhibitory")

for (ct in cell_types) {
  cat(sprintf("\n=== [R DESeq2] Processing: %s ===\n", ct))
  counts_file <- file.path(pb_dir, sprintf("counts_%s.tsv", ct))
  coldata_file <- file.path(pb_dir, sprintf("coldata_%s.tsv", ct))
  
  if (!file.exists(counts_file) || !file.exists(coldata_file)) {
    cat(sprintf("Skipping %s (missing files)\n", ct))
    next
  }
  
  counts <- as.matrix(read.table(counts_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE))
  coldata <- read.table(coldata_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
  
  coldata <- coldata[colnames(counts), ]
  coldata$diagnosis <- factor(coldata$diagnosis, levels = c("Control", "PD"))
  coldata$sex <- factor(coldata$sex)
  coldata$age <- as.numeric(coldata$age)
  
  keep <- rowSums(counts) >= 10
  counts <- counts[keep, ]
  
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = coldata,
                                design = ~ sex + age + diagnosis)
  
  dds <- DESeq(dds, quiet=TRUE)
  res <- results(dds, contrast=c("diagnosis", "PD", "Control"))
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  out_csv <- file.path(out_dir, sprintf("deseq2_R_%s.tsv", ct))
  write.table(res_df, file=out_csv, sep="\t", quote=FALSE, row.names=FALSE)
  
  sig_nom <- sum(!is.na(res_df$pvalue) & res_df$pvalue < 0.05)
  sig_fdr <- sum(!is.na(res_df$padj) & res_df$padj < 0.1)
  cat(sprintf("R DESeq2 %s: %d tested | %d nominal (p<0.05) | %d FDR (padj<0.1)\n",
              ct, nrow(res_df), sig_nom, sig_fdr))
}
