#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr

base_dir = "path/to/project_root"
r_dir = os.path.join(base_dir, "results/de_deseq2_r")
py_dir = os.path.join(base_dir, "results/de_pydeseq2_py")
out_dir = os.path.join(base_dir, "results/de_comparison")

cell_types = ["Microglia", "Astrocytes", "Oligodendrocytes", "OPCs", "Excitatory", "Inhibitory"]
rows = []

for ct in cell_types:
    r_file = os.path.join(r_dir, f"deseq2_R_{ct}.tsv")
    py_file = os.path.join(py_dir, f"pydeseq2_py_{ct}.tsv")
    
    if not (os.path.exists(r_file) and os.path.exists(py_file)):
        continue
        
    df_r = pd.read_csv(r_file, sep="\t").set_index("gene")
    df_py = pd.read_csv(py_file, sep="\t").set_index("gene")
    
    common = df_r.index.intersection(df_py.index)
    df_r = df_r.loc[common].dropna(subset=["log2FoldChange", "pvalue", "stat"])
    df_py = df_py.loc[common].dropna(subset=["log2FoldChange", "pvalue", "stat"])
    valid = df_r.index.intersection(df_py.index)
    
    sub_r = df_r.loc[valid]
    sub_py = df_py.loc[valid]
    
    r_lfc, _ = pearsonr(sub_r["log2FoldChange"], sub_py["log2FoldChange"])
    rho_lfc, _ = spearmanr(sub_r["log2FoldChange"], sub_py["log2FoldChange"])
    r_stat, _ = pearsonr(sub_r["stat"], sub_py["stat"])
    
    dir_concordance = np.mean(np.sign(sub_r["log2FoldChange"]) == np.sign(sub_py["log2FoldChange"])) * 100
    
    sig_r_nom = set(sub_r[sub_r["pvalue"] < 0.05].index)
    sig_py_nom = set(sub_py[sub_py["pvalue"] < 0.05].index)
    inter_nom = len(sig_r_nom & sig_py_nom)
    union_nom = len(sig_r_nom | sig_py_nom)
    jaccard_nom = inter_nom / union_nom if union_nom > 0 else 1.0
    
    sig_r_fdr = set(sub_r[sub_r["padj"] < 0.1].index)
    sig_py_fdr = set(sub_py[sub_py["padj"] < 0.1].index)
    inter_fdr = len(sig_r_fdr & sig_py_fdr)
    union_fdr = len(sig_r_fdr | sig_py_fdr)
    jaccard_fdr = inter_fdr / union_fdr if union_fdr > 0 else 1.0
    
    rows.append({
        "cell_type": ct,
        "n_genes": len(valid),
        "r_log2FC": round(r_lfc, 4),
        "rho_log2FC": round(rho_lfc, 4),
        "r_stat": round(r_stat, 4),
        "direction_agreement_%": round(dir_concordance, 2),
        "R_nom": len(sig_r_nom),
        "Py_nom": len(sig_py_nom),
        "overlap_nom": inter_nom,
        "jaccard_nom": round(jaccard_nom, 4),
        "R_fdr0.1": len(sig_r_fdr),
        "Py_fdr0.1": len(sig_py_fdr),
        "overlap_fdr": inter_fdr,
        "jaccard_fdr": round(jaccard_fdr, 4)
    })

res_df = pd.DataFrame(rows)
res_df.to_csv(os.path.join(out_dir, "deseq2_vs_pydeseq2_summary.tsv"), sep="\t", index=False)

print("\n" + "="*85)
print("     COMPARISON SUMMARY: R DESeq2 vs PyDESeq2 across Cell Types")
print("="*85)
print(res_df.to_string(index=False))
