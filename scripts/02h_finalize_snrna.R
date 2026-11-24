#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
})

base_dir <- "path/to/project_root"
results_dir <- file.path(base_dir, "results")
pb_dir <- file.path(results_dir, "pseudobulk")
qc_dir <- file.path(results_dir, "qc")
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== [1/2] Generating Pseudobulk PCA Diagnostic Plots ===\n")
# Check PCA across the top abundant cell types to verify diagnosis vs technical variance
test_cts <- c("Oligodendrocytes", "Astrocytes", "Microglia", "Excitatory", "Inhibitory")

pdf(file.path(qc_dir, "pseudobulk_pca_diagnostics.pdf"), width=10, height=7)
par(mfrow=c(2, 3))

for (ct in test_cts) {
  counts_file <- file.path(pb_dir, sprintf("counts_%s.tsv", ct))
  coldata_file <- file.path(pb_dir, sprintf("coldata_%s.tsv", ct))
  
  if (!file.exists(counts_file)) next
  
  counts <- as.matrix(read.table(counts_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE))
  coldata <- read.table(coldata_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
  coldata <- coldata[colnames(counts), ]
  coldata$diagnosis <- factor(coldata$diagnosis, levels = c("Control", "PD"))
  
  # Filter low counts
  keep <- rowSums(counts) >= 10
  counts <- counts[keep, ]
  
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ diagnosis)
  vsd <- vst(dds, blind = TRUE)
  
  pca_data <- plotPCA(vsd, intgroup = c("diagnosis", "sex"), returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  
  plot(pca_data$PC1, pca_data$PC2, 
       col = ifelse(pca_data$diagnosis == "PD", "#D95F02", "#1B9E77"),
       pch = ifelse(pca_data$sex == "M", 16, 17),
       cex = 1.8,
       xlab = paste0("PC1: ", percentVar[1], "% variance"),
       ylab = paste0("PC2: ", percentVar[2], "% variance"),
       main = paste0(ct, " (PCA)"))
  legend("topright", legend = c("Control", "PD", "Male", "Female"),
         col = c("#1B9E77", "#D95F02", "black", "black"),
         pch = c(15, 15, 16, 17), bty = "n", cex = 0.8)
}
dev.off()
cat(sprintf("PCA diagnostics saved to: %s\n", file.path(qc_dir, "pseudobulk_pca_diagnostics.pdf")))

cat("\n=== [2/2] Mapping Gene Symbols to results/snrna_gene_scores.tsv ===\n")
scores_file <- file.path(results_dir, "snrna_gene_scores.tsv")
scores <- read.table(scores_file, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)

# Use AnnotationDbi / org.Hs.eg.db if installed, or biomaRt, or fallback offline map
has_org <- suppressWarnings(requireNamespace("org.Hs.eg.db", quietly = TRUE))

if (has_org) {
  library(org.Hs.eg.db)
  cat("Querying org.Hs.eg.db for Ensembl -> HGNC Symbol...\n")
  clean_ens <- gsub("\\..*$", "", scores$ensembl_id)
  syms <- mapIds(org.Hs.eg.db, keys = clean_ens, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  scores$symbol <- ifelse(!is.na(syms), syms, scores$ensembl_id)
} else {
  cat("org.Hs.eg.db not detected; attempting biomaRt query...\n")
  has_bm <- suppressWarnings(requireNamespace("biomaRt", quietly = TRUE))
  if (has_bm) {
    tryCatch({
      library(biomaRt)
      ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      clean_ens <- unique(gsub("\\..*$", "", scores$ensembl_id))
      bm_res <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                      filters = "ensembl_gene_id", values = clean_ens, mart = ensembl)
      sym_map <- setNames(bm_res$hgnc_symbol, bm_res$ensembl_gene_id)
      clean_ids <- gsub("\\..*$", "", scores$ensembl_id)
      mapped <- sym_map[clean_ids]
      scores$symbol <- ifelse(!is.na(mapped) & mapped != "", mapped, scores$ensembl_id)
    }, error = function(e) {
      cat("BiomaRt offline or timed out; keeping existing symbol column.\n")
    })
  } else {
    cat("No annotation package installed in deseq2_env. Symbol column left intact.\n")
  }
}

write.table(scores, file = scores_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Updated snrna_gene_scores.tsv with mapped symbols (%d total genes).\n", nrow(scores)))
cat("Preview top 10:\n")
print(head(scores[, c("ensembl_id", "symbol", "score", "direction", "effect", "p", "padj", "notes")], 10))
