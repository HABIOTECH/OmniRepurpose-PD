#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

base_dir = "path/to/project_root"
pb_dir = os.path.join(base_dir, "results/pseudobulk")
out_dir = os.path.join(base_dir, "results/de_pydeseq2_py")
os.makedirs(out_dir, exist_ok=True)

cell_types = ["Microglia", "Astrocytes", "Oligodendrocytes", "OPCs", "Excitatory", "Inhibitory"]

for ct in cell_types:
    print(f"\n=== [PyDESeq2] Processing: {ct} ===")
    counts_file = os.path.join(pb_dir, f"counts_{ct}.tsv")
    coldata_file = os.path.join(pb_dir, f"coldata_{ct}.tsv")
    
    if not (os.path.exists(counts_file) and os.path.exists(coldata_file)):
        continue
        
    counts = pd.read_csv(counts_file, sep="\t", index_col=0).T  # Samples x Genes
    coldata = pd.read_csv(coldata_file, sep="\t", index_col=0)
    
    # Pre-filter low count genes (matching R DESeq2: sum >= 10)
    keep = counts.sum(axis=0) >= 10
    counts = counts.loc[:, keep].astype(int)
    coldata = coldata.loc[counts.index]
    
    # Center age to improve GLM convergence
    coldata["age"] = coldata["age"].astype(float)
    coldata["age"] = coldata["age"] - coldata["age"].mean()
    coldata["sex"] = coldata["sex"].astype(str)
    coldata["diagnosis"] = coldata["diagnosis"].astype(str)
    
    dds = DeseqDataSet(
        counts=counts,
        metadata=coldata,
        design_factors=["sex", "age", "diagnosis"],
        refit_cooks=True,
        quiet=True
    )
    dds.deseq2()
    
    stat_res = DeseqStats(dds, contrast=["diagnosis", "PD", "Control"], quiet=True)
    stat_res.summary()
    res_df = stat_res.results_df
    res_df["gene"] = res_df.index
    
    out_csv = os.path.join(out_dir, f"pydeseq2_py_{ct}.tsv")
    res_df.to_csv(out_csv, sep="\t", index=False)
    
    sig_nom = (res_df["pvalue"] < 0.05).sum()
    sig_fdr = (res_df["padj"] < 0.1).sum()
    print(f"PyDESeq2 {ct}: {len(res_df)} tested | {sig_nom} nominal (p<0.05) | {sig_fdr} FDR (padj<0.1)")
