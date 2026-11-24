#!/usr/bin/env python3
import sys
import os
import time
import pandas as pd
import numpy as np
from scipy import sparse
import anndata as ad

def main():
    base_dir = "path/to/project_root"
    cell_file = os.path.join(base_dir, "GSE157783_IPDCO_hg_midbrain_cell/IPDCO_hg_midbrain_cell.tsv")
    gene_file = os.path.join(base_dir, "GSE157783_IPDCO_hg_midbrain_genes/IPDCO_hg_midbrain_genes.tsv")
    umi_file  = os.path.join(base_dir, "GSE157783_IPDCO_hg_midbrain_UMI/IPDCO_hg_midbrain_UMI.tsv")
    out_h5ad  = os.path.join(base_dir, "GSE157783_IPDCO_midbrain.h5ad")

    print("[1/4] Reading cell & gene metadata...")
    obs = pd.read_csv(cell_file, sep="\t")
    obs.set_index("barcode", inplace=True)

    # Annotate clinical covariates directly from donor table
    donor_meta = {
        "C1":  {"diagnosis": "Control", "sex": "M", "age": 67, "pmi": 10.5},
        "C2":  {"diagnosis": "Control", "sex": "M", "age": 74, "pmi": 14.5},
        "C3":  {"diagnosis": "Control", "sex": "F", "age": 79, "pmi": 17.5},
        "C4":  {"diagnosis": "Control", "sex": "F", "age": 82, "pmi": 16.0},
        "C5":  {"diagnosis": "Control", "sex": "M", "age": 60, "pmi": 9.0},
        "C6":  {"diagnosis": "Control", "sex": "M", "age": 71, "pmi": 21.0},
        "PD1": {"diagnosis": "PD",      "sex": "M", "age": 72, "pmi": 11.5},
        "PD2": {"diagnosis": "PD",      "sex": "F", "age": 83, "pmi": 12.0},
        "PD3": {"diagnosis": "PD",      "sex": "M", "age": 76, "pmi": 13.0},
        "PD4": {"diagnosis": "PD",      "sex": "M", "age": 81, "pmi": 18.0},
        "PD5": {"diagnosis": "PD",      "sex": "M", "age": 79, "pmi": 14.0},
    }
    
    donor_df = pd.DataFrame.from_dict(donor_meta, orient="index")
    obs = obs.join(donor_df, on="patient")

    var = pd.read_csv(gene_file, sep="\t")
    var.set_index("gene", inplace=True)
    var.sort_values(by="row", inplace=True)

    n_cells = len(obs)
    n_genes = len(var)
    print(f"       Total cells: {n_cells}, Total genes: {n_genes}")

    print("[2/4] Streaming dense Gene x Cell matrix into sparse CSR (Cell x Gene)...")
    t0 = time.time()
    chunksize = 2500
    row_chunks = []

    for i, chunk in enumerate(pd.read_csv(umi_file, sep="\t", chunksize=chunksize, dtype=np.float32)):
        sp_chunk = sparse.csr_matrix(chunk.values)
        row_chunks.append(sp_chunk)
        processed = min((i + 1) * chunksize, n_genes)
        print(f"       Chunk {i+1}: processed {processed} / {n_genes} genes ({time.time() - t0:.1f}s)")

    print("       Stacking rows and transposing...")
    gene_by_cell = sparse.vstack(row_chunks, format="csr")
    cell_by_gene = gene_by_cell.T.tocsr()

    print("[3/4] Assembling AnnData object...")
    adata = ad.AnnData(X=cell_by_gene, obs=obs, var=var)
    adata.raw = adata  # keep raw counts preserved

    print(f"[4/4] Writing compressed .h5ad to: {out_h5ad}")
    adata.write_h5ad(out_h5ad, compression="gzip")
    print(f"Done! Shape: {adata.shape} in {time.time() - t0:.1f} seconds.")

if __name__ == "__main__":
    main()
