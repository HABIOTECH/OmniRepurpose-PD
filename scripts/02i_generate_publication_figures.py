#!/usr/bin/env python3
"""
Generate publication-quality figures for GSE157783 snRNA-seq midbrain analysis.
Outputs 300 DPI PNGs and vector PDFs in results/figures/.
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Setup publication style
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 15,
    'figure.dpi': 300,
    'axes.edgecolor': '#333333',
    'axes.linewidth': 0.8,
    'grid.color': '#E0E0E0',
    'grid.linestyle': '--',
    'grid.linewidth': 0.5,
})

base_dir = "path/to/project_root"
fig_dir = os.path.join(base_dir, "results/figures")
os.makedirs(fig_dir, exist_ok=True)

# -------------------------------------------------------------------------
# Figure 1: Cohort Composition & Nuclei per Cell Type
# -------------------------------------------------------------------------
print("Generating Figure 1: Cell Type & Donor Composition...")
summary_file = os.path.join(base_dir, "results/pseudobulk/pseudobulk_summary.tsv")
if os.path.exists(summary_file):
    df_sum = pd.read_csv(summary_file, sep="\t").sort_values(by="total_cells", ascending=True)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5), gridspec_kw={'width_ratios': [1.2, 1]})
    
    # 1A: Total nuclei per cell type
    bars = ax1.barh(df_sum["cell_type"], df_sum["total_cells"], color="#2B5C8F", edgecolor="#1C3D5E", alpha=0.9)
    for bar in bars:
        w = bar.get_width()
        ax1.text(w + 300, bar.get_y() + bar.get_height()/2, f"{w:,}", va='center', ha='left', fontsize=9, color="#222222")
    ax1.set_xlabel("Number of Nuclei")
    ax1.set_title("A. Midbrain Nuclei Distribution (Kamath et al. GSE157783)", fontweight="bold")
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.set_xlim(0, max(df_sum["total_cells"]) * 1.18)
    
    # 1B: Donor Representation (Control vs PD)
    y_pos = np.arange(len(df_sum))
    p1 = ax2.barh(y_pos, df_sum["n_ctrl"], color="#1B9E77", label="Control (N=6)", alpha=0.85)
    p2 = ax2.barh(y_pos, df_sum["n_pd"], left=df_sum["n_ctrl"], color="#D95F02", label="PD Case (N=5)", alpha=0.85)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(df_sum["cell_type"])
    ax2.set_xlabel("Donors Represented (>= 5 nuclei)")
    ax2.set_title("B. Donor Coverage per Subpopulation", fontweight="bold")
    ax2.legend(loc="lower right", frameon=True)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_xlim(0, 12)
    
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "Fig1_celltype_composition.png"), dpi=300)
    plt.savefig(os.path.join(fig_dir, "Fig1_celltype_composition.pdf"))
    plt.close()

# -------------------------------------------------------------------------
# Figure 2: Differential Expression Discovery Yield across Cell Types
# -------------------------------------------------------------------------
print("Generating Figure 2: DE Discovery Yield...")
r_dir = os.path.join(base_dir, "results/de_deseq2_r")
de_files = glob.glob(os.path.join(r_dir, "deseq2_R_*.tsv"))

de_summary = []
for f in de_files:
    ct = os.path.basename(f).replace("deseq2_R_", "").replace(".tsv", "")
    df = pd.read_csv(f, sep="\t")
    n_nom = (df["pvalue"] < 0.05).sum()
    n_fdr = (df["padj"] < 0.10).sum()
    de_summary.append({"cell_type": ct, "nominal": n_nom, "fdr": n_fdr, "total": len(df)})

df_de = pd.DataFrame(de_summary).sort_values(by="fdr", ascending=True)

fig, ax = plt.subplots(figsize=(8.5, 6))
y = np.arange(len(df_de))
bar_nom = ax.barh(y - 0.2, df_de["nominal"], height=0.38, color="#A6CEE3", label="Nominal (p < 0.05)", edgecolor="#1F78B4")
bar_fdr = ax.barh(y + 0.2, df_de["fdr"], height=0.38, color="#E31A1C", label="FDR Significant (padj < 0.10)", edgecolor="#800000")

for bar in bar_fdr:
    w = bar.get_width()
    if w > 0:
        ax.text(w + 1.5, bar.get_y() + bar.get_height()/2, f"{int(w)}", va='center', ha='left', fontsize=9.5, fontweight="bold", color="#800000")

ax.set_yticks(y)
ax.set_yticklabels(df_de["cell_type"])
ax.set_xlabel("Number of Differentially Expressed Genes (PD vs Control)")
ax.set_title("Midbrain snRNA-seq Differential Expression Yield by Cell Type", fontweight="bold", pad=12)
ax.legend(loc="lower right", frameon=True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "Fig2_de_discovery_yield.png"), dpi=300)
plt.savefig(os.path.join(fig_dir, "Fig2_de_discovery_yield.pdf"))
plt.close()

# -------------------------------------------------------------------------
# Figure 3: Hallmark Parkinson's Disease Genes Profiling
# -------------------------------------------------------------------------
print("Generating Figure 3: Hallmark PD Genes...")
scores_file = os.path.join(base_dir, "results/snrna_gene_scores.tsv")
if os.path.exists(scores_file):
    df_scores = pd.read_csv(scores_file, sep="\t")
    pd_genes = ["SNCA", "LRRK2", "GBA1", "PRKN", "PINK1", "MAPT"]
    df_pd = df_scores[df_scores["symbol"].isin(pd_genes)].copy()
    
    if len(df_pd) > 0:
        fig, ax = plt.subplots(figsize=(8, 4.5))
        df_pd = df_pd.sort_values(by="score", ascending=True)
        colors = ["#D95F02" if x < 0 else "#2B5C8F" for x in df_pd["effect"]]
        
        bars = ax.barh(df_pd["symbol"], df_pd["score"], color=colors, edgecolor="#222222", alpha=0.85)
        ax.axvline(-np.log10(0.05), color="gray", linestyle="--", linewidth=1, label="Nominal significance (p = 0.05)")
        ax.axvline(-np.log10(0.01), color="darkred", linestyle=":", linewidth=1, label="FDR threshold (p = 0.01)")
        
        for bar, (_, row) in zip(bars, df_pd.iterrows()):
            ct_note = row["notes"].split(";")[0].replace("celltype=", "")
            eff = f"log2FC={row['effect']:+.2f} ({ct_note})"
            ax.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height()/2, eff, va='center', ha='left', fontsize=9.5)
            
        ax.set_xlabel(r"Evidence Score ($-\log_{10} p$)")
        ax.set_title("Expression Alteration of Hallmark PD Risk Genes in Human Midbrain", fontweight="bold", pad=12)
        ax.legend(loc="lower right", frameon=True)
        ax.set_xlim(0, max(df_pd["score"]) * 1.55)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, "Fig3_hallmark_pd_genes.png"), dpi=300)
        plt.savefig(os.path.join(fig_dir, "Fig3_hallmark_pd_genes.pdf"))
        plt.close()

# -------------------------------------------------------------------------
# Figure 4: Volcano Plot for Dopaminergic Neurons (DaNs)
# -------------------------------------------------------------------------
print("Generating Figure 4: Volcano Plot for DaNs...")
dan_file = os.path.join(r_dir, "deseq2_R_DaNs.tsv")
id_map_file = os.path.join(base_dir, "data/ref/id_map.tsv")

if os.path.exists(dan_file):
    df_dan = pd.read_csv(dan_file, sep="\t").dropna(subset=["pvalue", "log2FoldChange"])
    if os.path.exists(id_map_file):
        id_map = pd.read_csv(id_map_file, sep="\t").set_index("ensembl_id")["symbol"].to_dict()
        df_dan["symbol"] = df_dan["gene"].map(lambda g: id_map.get(g, g))
    else:
        df_dan["symbol"] = df_dan["gene"]
        
    df_dan["neg_log10_p"] = -np.log10(df_dan["pvalue"].clip(lower=1e-20))
    
    fig, ax = plt.subplots(figsize=(8, 6.5))
    
    # Classify points
    sig = (df_dan["padj"] < 0.10) & (df_dan["log2FoldChange"].abs() >= 1.0)
    up = sig & (df_dan["log2FoldChange"] > 0)
    down = sig & (df_dan["log2FoldChange"] < 0)
    ns = ~sig
    
    ax.scatter(df_dan.loc[ns, "log2FoldChange"], df_dan.loc[ns, "neg_log10_p"], c="#B0B0B0", s=14, alpha=0.4, label="Not significant")
    ax.scatter(df_dan.loc[up, "log2FoldChange"], df_dan.loc[up, "neg_log10_p"], c="#E41A1C", s=24, alpha=0.85, label="Upregulated (FDR < 0.1, |log2FC| >= 1)")
    ax.scatter(df_dan.loc[down, "log2FoldChange"], df_dan.loc[down, "neg_log10_p"], c="#377EB8", s=24, alpha=0.85, label="Downregulated (FDR < 0.1, |log2FC| >= 1)")
    
    ax.axhline(-np.log10(0.05), color="gray", linestyle="--", linewidth=0.8)
    ax.axvline(-1.0, color="gray", linestyle=":", linewidth=0.8)
    ax.axvline(1.0, color="gray", linestyle=":", linewidth=0.8)
    
    # Label top hits and SNCA
    label_candidates = df_dan[df_dan["symbol"] == "SNCA"]
    top_candidates = df_dan[sig].sort_values(by="neg_log10_p", ascending=False).head(12)
    to_label = pd.concat([label_candidates, top_candidates]).drop_duplicates(subset=["gene"])
    
    for _, r in to_label.iterrows():
        ax.annotate(r["symbol"], (r["log2FoldChange"], r["neg_log10_p"]),
                    xytext=(5, 4), textcoords="offset points",
                    fontsize=9, fontweight="bold", color="#111111")
        
    ax.set_xlabel(r"$\log_2$ Fold Change (PD vs Control)")
    ax.set_ylabel(r"$-\log_{10}(P\text{-value})$")
    ax.set_title("Differential Expression in Midbrain Dopaminergic Neurons (DaNs)", fontweight="bold", pad=12)
    ax.legend(loc="upper left", frameon=True, fontsize=9.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "Fig4_volcano_DaNs.png"), dpi=300)
    plt.savefig(os.path.join(fig_dir, "Fig4_volcano_DaNs.pdf"))
    plt.close()

print(f"\nAll publication figures successfully created in: {fig_dir}")
print("Figures available in both 300 DPI PNG and vector PDF formats.")
