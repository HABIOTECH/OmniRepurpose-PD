#!/usr/bin/env python3
"""
Generate results/snrna_gene_scores.tsv according to the pipeline contract:
ensembl_id, symbol, score, direction, effect, p, padj, notes
For each gene, select the strongest association (-log10 p) across cell types,
and record the cell type in notes.
"""

import os
import glob
import numpy as np
import pandas as pd

base_dir = "/sc/kzd/home/abdullah/bioinformatics/PK-multiomics"
r_dir = os.path.join(base_dir, "results/de_deseq2_r")
out_file = os.path.join(base_dir, "results/snrna_gene_scores.tsv")

all_files = glob.glob(os.path.join(r_dir, "deseq2_R_*.tsv"))
all_dfs = []

for fpath in all_files:
    ct_name = os.path.basename(fpath).replace("deseq2_R_", "").replace(".tsv", "")
    df = pd.read_csv(fpath, sep="\t")
    df = df.dropna(subset=["pvalue", "log2FoldChange"])
    df["cell_type"] = ct_name
    all_dfs.append(df)

combined = pd.concat(all_dfs, ignore_index=True)
combined["pvalue"] = combined["pvalue"].clip(lower=1e-300)
combined["score"] = -np.log10(combined["pvalue"])

# Select top hit per gene across all cell types
best_hits = combined.sort_values(by="score", ascending=False).groupby("gene").first().reset_index()

# Attempt to map gene symbols via mygene or fallback to ensembl_id
ensembl_list = best_hits["gene"].tolist()
symbol_map = {}
try:
    import mygene
    mg = mygene.MyGeneInfo()
    res = mg.querymany(ensembl_list, scopes='ensembl.gene', fields='symbol', species='human', as_dataframe=True)
    if 'symbol' in res.columns:
        symbol_map = res['symbol'].dropna().to_dict()
except Exception:
    pass

contract_df = pd.DataFrame()
contract_df["ensembl_id"] = best_hits["gene"]
contract_df["symbol"] = best_hits["gene"].map(lambda g: symbol_map.get(g, g))
contract_df["score"] = best_hits["score"].round(4)
contract_df["direction"] = np.sign(best_hits["log2FoldChange"]).astype(int)
contract_df["effect"] = best_hits["log2FoldChange"].round(4)
contract_df["p"] = best_hits["pvalue"]
contract_df["padj"] = best_hits["padj"].fillna(1.0)
contract_df["notes"] = "celltype=" + best_hits["cell_type"] + ";baseMean=" + best_hits["baseMean"].round(1).astype(str)

contract_df.to_csv(out_file, sep="\t", index=False)
print(f"\n[Contract Deliverable Generated]: {out_file}")
print(f"Total unique genes covered: {len(contract_df)}")
print(f"Top 15 strongest midbrain associations across cell types:")
print(contract_df.head(15)[["ensembl_id", "symbol", "score", "direction", "effect", "p", "padj", "notes"]].to_string(index=False))
