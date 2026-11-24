#!/usr/bin/env python3
"""
Pseudobulk aggregation: sum raw integer counts per (donor, cell_type).
Generates count matrices and sample metadata tables ready for DESeq2/edgeR/pydeseq2.
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

def main():
    base_dir = "path/to/project_root"
    h5ad_file = os.path.join(base_dir, "GSE157783_IPDCO_midbrain.h5ad")
    out_dir   = os.path.join(base_dir, "results/pseudobulk")
    os.makedirs(out_dir, exist_ok=True)

    print("[1/3] Loading AnnData (raw counts)...")
    adata = ad.read_h5ad(h5ad_file)

    cell_types = adata.obs["cell_ontology"].unique()
    donors     = adata.obs["patient"].unique()
    
    print(f"      Found {len(cell_types)} cell types and {len(donors)} donors.")

    donor_meta = adata.obs[["patient", "diagnosis", "sex", "age", "pmi"]].drop_duplicates().set_index("patient")

    print("[2/3] Aggregating raw counts by (donor x cell_type)...")
    summary_rows = []

    for ct in cell_types:
        ct_clean = ct.replace(" ", "_").replace("+", "pos")
        ct_cells = adata[adata.obs["cell_ontology"] == ct]
        
        # Check donor representation
        donor_counts = ct_cells.obs["patient"].value_counts()
        # Keep donors with at least 5 nuclei for this cell type
        valid_donors = donor_counts[donor_counts >= 5].index.tolist()

        if len(valid_donors) < 4:
            print(f"      Skipping {ct}: only {len(valid_donors)} donors have >= 5 cells.")
            continue

        pb_counts = pd.DataFrame(index=adata.var_names, columns=valid_donors, dtype=int)
        
        for donor in valid_donors:
            donor_ct_cells = ct_cells[ct_cells.obs["patient"] == donor]
            # Sum raw counts
            summed = np.asarray(donor_ct_cells.X.sum(axis=0)).ravel()
            pb_counts[donor] = np.round(summed).astype(int)

        # Write cell-type count matrix and matching coldata
        counts_path = os.path.join(out_dir, f"counts_{ct_clean}.tsv")
        coldata_path = os.path.join(out_dir, f"coldata_{ct_clean}.tsv")

        pb_counts.to_csv(counts_path, sep="\t")
        sub_meta = donor_meta.loc[valid_donors]
        sub_meta["n_cells"] = donor_counts.loc[valid_donors]
        sub_meta.to_csv(coldata_path, sep="\t")

        summary_rows.append({
            "cell_type": ct,
            "n_donors": len(valid_donors),
            "n_ctrl": sum(sub_meta["diagnosis"] == "Control"),
            "n_pd": sum(sub_meta["diagnosis"] == "PD"),
            "total_cells": int(sub_meta["n_cells"].sum())
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(os.path.join(out_dir, "pseudobulk_summary.tsv"), sep="\t", index=False)
    print("[3/3] Done aggregating pseudobulks!")
    print(summary_df.to_string(index=False))

if __name__ == "__main__":
    main()
