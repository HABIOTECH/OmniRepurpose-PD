#!/usr/bin/env python3
import pandas as pd

base = "/sc/kzd/home/abdullah/bioinformatics/PK-multiomics"
map_file = f"{base}/data/ref/id_map.tsv"
snrna_file = f"{base}/results/snrna_gene_scores.tsv"

print("Loading ID map and snRNA results...")
id_map = pd.read_csv(map_file, sep="\t").set_index("ensembl_id")
snrna = pd.read_csv(snrna_file, sep="\t")

# Map symbols and biotypes
clean_ids = snrna["ensembl_id"].str.replace(r"\..*$", "", regex=True)
mapped_symbols = clean_ids.map(id_map["symbol"])
mapped_biotypes = clean_ids.map(id_map["gene_biotype"])

# Update symbol column with official HGNC symbols where available
snrna["symbol"] = mapped_symbols.fillna(snrna["ensembl_id"])
snrna["gene_biotype"] = mapped_biotypes.fillna("unknown")

# Keep sorted descending by evidence score
snrna = snrna.sort_values(by="score", ascending=False)
snrna.to_csv(snrna_file, sep="\t", index=False)

print("\nSuccessfully updated results/snrna_gene_scores.tsv!")
print("Top 15 strongest hits with official symbols:")
print(snrna[["ensembl_id", "symbol", "gene_biotype", "score", "direction", "effect", "p", "padj", "notes"]].head(15).to_string(index=False))

print("\nChecking landmark PD genes:")
pd_symbols = ["SNCA", "LRRK2", "GBA1", "PRKN", "PINK1", "MAPT"]
print(snrna[snrna["symbol"].isin(pd_symbols)][["ensembl_id", "symbol", "gene_biotype", "score", "direction", "effect", "p", "padj", "notes"]].to_string(index=False))
