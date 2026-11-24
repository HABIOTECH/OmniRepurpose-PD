#!/usr/bin/env python3
import os
import pandas as pd

base_dir = "path/to/project_root"
snrna_file = os.path.join(base_dir, "results/snrna_gene_scores.tsv")

df = pd.read_csv(snrna_file, sep="\t")

# Re-sort descending by score
df = df.sort_values(by="score", ascending=False)
df.to_csv(snrna_file, sep="\t", index=False)

print("="*90)
print("TOP 20 MOST SIGNIFICANT HITS ACROSS ALL MIDBRAIN CELL TYPES")
print("="*90)
print(df.head(20)[["ensembl_id", "score", "direction", "effect", "p", "padj", "notes"]].to_string(index=False))

# Check known key Parkinson's genes
# SNCA: ENSG00000145335, LRRK2: ENSG00000188906, GBA1: ENSG00000177628, 
# PRKN: ENSG00000185345, PINK1: ENSG00000158828, MAPT: ENSG00000186868, TH: ENSG00000180172
pd_genes = {
    "SNCA": "ENSG00000145335",
    "LRRK2": "ENSG00000188906",
    "GBA1": "ENSG00000177628",
    "PRKN": "ENSG00000185345",
    "PINK1": "ENSG00000158828",
    "MAPT": "ENSG00000186868",
    "TH": "ENSG00000180172"
}

print("\n" + "="*90)
print("KNOWN PD GENES CHECK IN MIDBRAIN snRNA-seq RESULTS")
print("="*90)
pd_rows = []
for sym, ens in pd_genes.items():
    match = df[df["ensembl_id"] == ens]
    if len(match) > 0:
        row = match.iloc[0].to_dict()
        row["gene_name"] = sym
        pd_rows.append(row)

if pd_rows:
    pd_df = pd.DataFrame(pd_rows)
    print(pd_df[["gene_name", "ensembl_id", "score", "direction", "effect", "p", "padj", "notes"]].to_string(index=False))
else:
    print("No matching PD genes found in snRNA-seq table.")
