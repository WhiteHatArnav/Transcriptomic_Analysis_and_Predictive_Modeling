import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

BLACK_DIR = "/Users/arnavjoshi/Desktop/BlackSamplesPCA/PRAD_PCA"
WHITE_DIR = "/Users/arnavjoshi/Desktop/WhiteSamplesPCA/PRAD_PCA"
COMB_DIR  = "/Users/arnavjoshi/Desktop/PCAWhiteBlackCombined/PRAD"

SCORES_FILE = "PRAD_per_sample_PCA_scores.csv"
BLACK_SCORES = os.path.join(BLACK_DIR, SCORES_FILE)
WHITE_SCORES = os.path.join(WHITE_DIR, SCORES_FILE)

# For single PCA we need the expression files
BLACK_INPUT = "/Users/arnavjoshi/Desktop/PRADPredictiveAnalytics/prad_black_survival_data.csv"
WHITE_INPUT = "/Users/arnavjoshi/Desktop/PRADPredictiveAnalytics/prad_white_survival_data.csv"

def safe_load_scores(path):
    df = pd.read_csv(path)
    for col in ["PC1","PC2","OS"]:
        if col not in df.columns:
            raise ValueError(f"{col} missing in {path}")
    df["OS"] = pd.to_numeric(df["OS"], errors="coerce").fillna(0).astype(int)
    return df

def load_expr(path):
    df = pd.read_csv(path)
    if df.columns[0] != "sample" and df.columns[0] != "Patient ID":
        df = df.rename(columns={df.columns[0]: "sample"})
    df.columns = df.columns.str.strip()
    if "OS" not in df.columns or "OS.time" not in df.columns:
        raise ValueError(f"Missing OS/OS.time in {path}")
    df = df.dropna(subset=["OS","OS.time"]).copy()
    df["OS"] = pd.to_numeric(df["OS"], errors="coerce").fillna(0).astype(int)
    return df

def plot_overlay_separate_pca():
    os.makedirs(COMB_DIR, exist_ok=True)
    bd = safe_load_scores(BLACK_SCORES)
    wd = safe_load_scores(WHITE_SCORES)

    masks = {
        "White, OS=0": wd["OS"]==0,
        "White, OS=1": wd["OS"]==1,
        "Black, OS=0": bd["OS"]==0,
        "Black, OS=1": bd["OS"]==1,
    }
    plt.figure(figsize=(7.6, 6.4))
    for label, mask in masks.items():
        df = wd if "White" in label else bd
        if mask.any():
            plt.scatter(
                df.loc[mask,"PC1"], df.loc[mask,"PC2"],
                c="#1f77b4" if "OS=0" in label else "#d62728",
                marker="o" if "White" in label else "^",
                alpha=0.9, edgecolor="k", linewidth=0.3, s=50,
                label=label
            )
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("PRAD: PC1 vs PC2 — Overlay of Separate PCAs")
    plt.legend(frameon=True, loc="best")
    plt.tight_layout()
    outpath = os.path.join(COMB_DIR, "PRAD_WhiteBlack_PC1_PC2_combined.png")
    plt.savefig(outpath, dpi=300)
    plt.close()
    print(f"Saved: {outpath}")

def plot_single_pca_merged():
    os.makedirs(COMB_DIR, exist_ok=True)
    bd = load_expr(BLACK_INPUT)
    wd = load_expr(WHITE_INPUT)

    gene_cols = [c for c in bd.columns if c not in ["sample","Patient ID","OS","OS.time"]]
    gene_cols = [c for c in gene_cols if c in wd.columns]

    bd["Race"]="Black"; wd["Race"]="White"
    comb = pd.concat([bd,wd], ignore_index=True)
    X = comb[gene_cols].values
    Xs = StandardScaler().fit_transform(X)

    pca = PCA(n_components=2, random_state=42)
    pcs = pca.fit_transform(Xs)
    comb["PC1"] = pcs[:,0]; comb["PC2"] = pcs[:,1]

    groups = [
        ("White, OS=0", (comb["Race"]=="White") & (comb["OS"]==0)),
        ("White, OS=1", (comb["Race"]=="White") & (comb["OS"]==1)),
        ("Black, OS=0", (comb["Race"]=="Black") & (comb["OS"]==0)),
        ("Black, OS=1", (comb["Race"]=="Black") & (comb["OS"]==1)),
    ]

    plt.figure(figsize=(7.6, 6.4))
    for label, m in groups:
        if m.any():
            plt.scatter(
                comb.loc[m,"PC1"], comb.loc[m,"PC2"],
                c="#1f77b4" if "OS=0" in label else "#d62728",
                marker="o" if "White" in label else "^",
                alpha=0.9, edgecolor="k", linewidth=0.3, s=50, label=label
            )
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}% var)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}% var)")
    plt.title("PRAD: PC1 vs PC2 — Single PCA on Combined Cohort")
    plt.legend(frameon=True, loc="best")
    plt.tight_layout()
    outpath = os.path.join(COMB_DIR, "PRAD_WhiteBlack_PC1_PC2_singlePCA.png")
    plt.savefig(outpath, dpi=300)
    plt.close()
    print(f"Saved: {outpath}")

def main():
    plot_overlay_separate_pca()
    plot_single_pca_merged()

if __name__ == "__main__":
    main()
