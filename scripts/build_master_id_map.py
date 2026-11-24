#!/usr/bin/env python3
"""
Extract comprehensive human gene annotations from 10x Genomics GRCh38 GTF
to construct data/ref/id_map.tsv per OmniRepurpose-PD pipeline specifications.
"""

import gzip
import re
import os
import pandas as pd

gtf_path = "path/to/references/genes.gtf.gz"
out_path = "path/to/project_root/data/ref/id_map.tsv"

print(f"Reading and parsing GTF: {gtf_path}")

gene_records = []
re_gene_id   = re.compile(r'gene_id "([^"]+)"')
re_gene_name = re.compile(r'gene_name "([^"]+)"')
re_biotype   = re.compile(r'gene_type "([^"]+)"|gene_biotype "([^"]+)"')

with gzip.open(gtf_path, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        
        feature = parts[2]
        if feature != "gene":
            continue
        
        chrom  = parts[0]
        start  = parts[3]
        end    = parts[4]
        strand = parts[6]
        attrs  = parts[8]
        
        gid_m  = re_gene_id.search(attrs)
        gname_m = re_gene_name.search(attrs)
        bio_m  = re_biotype.search(attrs)
        
        if not gid_m:
            continue
        
        gene_id = gid_m.group(1)
        symbol = gname_m.group(1) if gname_m else gene_id
        biotype = (bio_m.group(1) or bio_m.group(2)) if bio_m else "unknown"
        
        gene_records.append({
            "ensembl_id": gene_id,
            "symbol": symbol,
            "gene_biotype": biotype,
            "chr": chrom,
            "start": int(start),
            "end": int(end),
            "strand": strand,
            "genome_build": "GRCh38"
        })

df_map = pd.DataFrame(gene_records).drop_duplicates(subset=["ensembl_id"])
df_map.to_csv(out_path, sep="\t", index=False)

print(f"Done! Master ID Map saved to: {out_path}")
print(f"Total annotated human genes: {len(df_map)}")
print("\nPreview of landmark PD genes:")
pd_genes = ["SNCA", "LRRK2", "GBA1", "PRKN", "PINK1", "MAPT"]
print(df_map[df_map["symbol"].isin(pd_genes)])
