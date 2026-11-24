#!/usr/bin/env python3
"""
Calculate gene-by-cell-type specificity matrix from GSE157783 for MAGMA gene-property analysis.
Specificity definition: Mean expression in cell type C / Sum of mean expression across all cell types.
Output: results/specificity.txt (rows = Ensembl gene IDs, columns = cell types).
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

def main():
    base_dir = "/sc/kzd/home/abdullah/bioinformatics/PK-multiomics"
    h5ad_file = os.path.join(base_dir, "GSE157783_IPDCO_midbrain.h5ad")
    out_file  = os.path.join(base_dir, "results/specificity.txt")

    print("[1/3] Loading AnnData...")
    adata = sc.read_h5ad(h5ad_file)
    print(f"      AnnData loaded: {adata.shape}")

    # Standard normalization for specificity: target 10,000 counts per cell
    print("[2/3] Normalizing counts (CP10k)...")
    sc.pp.normalize_total(adata, target_sum=1e4)

    # Calculate average expression per cell type
    print("      Calculating mean expression per cell type...")
    cell_types = adata.obs["cell_ontology"].unique()
    mean_expr = pd.DataFrame(index=adata.var_names, columns=cell_types, dtype=float)

    for ct in cell_types:
        subset = adata[adata.obs["cell_ontology"] == ct]
        # Mean across cells for each gene
        mean_expr[ct] = np.asarray(subset.X.mean(axis=0)).ravel()

    # Calculate specificity ratio: expr_c / sum(expr_all)
    print("[3/3] Computing specificity proportion and exporting...")
    row_sums = mean_expr.sum(axis=1)
    
    # Avoid divide-by-zero for non-expressed genes
    valid_genes = row_sums > 0
    specificity = pd.DataFrame(0.0, index=mean_expr.index, columns=mean_expr.columns)
    specificity.loc[valid_genes] = mean_expr.loc[valid_genes].div(row_sums[valid_genes], axis=0)

    # Save format expected by MAGMA (tab-delimited, ENSG in first column named GENE)
    specificity.index.name = "GENE"
    specificity.to_csv(out_file, sep="\t")
    print(f"Done! Specificity matrix written to: {out_file}")
    print(f"      Shape: {specificity.shape[0]} genes x {specificity.shape[1]} cell types")

if __name__ == "__main__":
    main()
