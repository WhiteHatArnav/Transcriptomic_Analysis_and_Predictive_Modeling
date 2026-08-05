#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# -----------------------------
# Folder definitions
# -----------------------------
RACE_PCA_DIRS = {
    "BRCA": {
        "black": "/Users/arnavjoshi/Desktop/BlackSamplesPCA/BRCA_PCA",
        "white": "/Users/arnavjoshi/Desktop/WhiteSamplesPCA/BRCA_PCA",
        "outfile_dir": "/Users/arnavjoshi/Desktop/PCAWhiteBlackCombined/BRCA"
    },
    "LIHC": {
        "black": "/Users/arnavjoshi/Desktop/BlackSamplesPCA/LIHC_PCA",
        "white": "/Users/arnavjoshi/Desktop/WhiteSamplesPCA/LIHC_PCA",
        "outfile_dir": "/Users/arnavjoshi/Desktop/PCAWhiteBlackCombined/LIHC"
    },
    "LUAD": {
        "black": "/Users/arnavjoshi/Desktop/BlackSamplesPCA/LUAD_PCA",
        "white": "/Users/arnavjoshi/Desktop/WhiteSamplesPCA/LUAD_PCA",
        "outfile_dir": "/Users/arnavjoshi/Desktop/PCAWhiteBlackCombined/LUAD"
    },
}

# -----------------------------
# Load function
# -----------------------------
def load_expression_matrix(path):
    df = pd.read_csv(path)
    df.rename(columns={df.columns[0]: "Patient ID"}, inplace=True)
    df.columns = df.columns.str.strip()
    if "OS" not in df.columns or "OS.time" not in df.columns:
        raise ValueError(f"Missing OS/OS.time in {path}")
    df = df.dropna(subset=["OS", "OS.time"])
    df["OS"] = pd.to_numeric(df["OS"], errors="coerce").fillna(0).astype(int)
    return df

# -----------------------------
# Combined PCA function
# -----------------------------
def combined_pca_plot(cancer, black_path, white_path, outfile_dir):
    os.makedirs(outfile_dir, exist_ok=True)

    black_df = load_expression_matrix(black_path)
    white_df = load_expression_matrix(white_path)

    # Select common features (genes)
    gene_cols = [c for c in black_df.columns if c not in ["Patient ID", "OS", "OS.time"]]
    gene_cols = [c for c in gene_cols if c in white_df.columns]

    # Add race label
    black_df["Race"] = "Black"
    white_df["Race"] = "White"

    combined_df = pd.concat([black_df, white_df], axis=0, ignore_index=True)
    X = combined_df[gene_cols].values
    X_scaled = StandardScaler().fit_transform(X)

    pca = PCA(n_components=2, random_state=42)
    pcs = pca.fit_transform(X_scaled)

    combined_df["PC1"] = pcs[:, 0]
    combined_df["PC2"] = pcs[:, 1]

    # 4 groups
    mask_w0 = (combined_df["Race"] == "White") & (combined_df["OS"] == 0)
    mask_w1 = (combined_df["Race"] == "White") & (combined_df["OS"] == 1)
    mask_b0 = (combined_df["Race"] == "Black") & (combined_df["OS"] == 0)
    mask_b1 = (combined_df["Race"] == "Black") & (combined_df["OS"] == 1)

    plt.figure(figsize=(7.6, 6.4))

    if mask_w0.any():
        plt.scatter(
            combined_df.loc[mask_w0, "PC1"], combined_df.loc[mask_w0, "PC2"],
            c="#1f77b4", alpha=0.8, edgecolor="k", linewidth=0.3, s=45,
            marker="o", label="White, OS=0"
        )
    if mask_w1.any():
        plt.scatter(
            combined_df.loc[mask_w1, "PC1"], combined_df.loc[mask_w1, "PC2"],
            c="#d62728", alpha=0.95, edgecolor="k", linewidth=0.3, s=55,
            marker="o", label="White, OS=1"
        )
    if mask_b0.any():
        plt.scatter(
            combined_df.loc[mask_b0, "PC1"], combined_df.loc[mask_b0, "PC2"],
            c="#1f77b4", alpha=0.95, edgecolor="k", linewidth=0.3, s=55,
            marker="^", label="Black, OS=0"
        )
    if mask_b1.any():
        plt.scatter(
            combined_df.loc[mask_b1, "PC1"], combined_df.loc[mask_b1, "PC2"],
            c="#d62728", alpha=1.0, edgecolor="k", linewidth=0.3, s=65,
            marker="^", label="Black, OS=1"
        )

    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}% variance)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}% variance)")
    plt.title(f"{cancer}: Combined PCA (White + Black, PC1 vs PC2)")
    plt.legend(frameon=True, loc="best")
    plt.tight_layout()

    outpath = os.path.join(outfile_dir, f"{cancer}_WhiteBlack_PC1_PC2_singlePCA.png")
    plt.savefig(outpath, dpi=300)
    plt.close()
    print(f"Saved combined PCA 2D: {outpath}")

# -----------------------------
# Main
# -----------------------------
def main():
    for cancer, d in RACE_PCA_DIRS.items():
        black_path = os.path.join(d["black"], f"{cancer}_per_sample_PCA_scores.csv")
        white_path = os.path.join(d["white"], f"{cancer}_per_sample_PCA_scores.csv")
        combined_pca_plot(cancer, black_path, white_path, d["outfile_dir"])

if __name__ == "__main__":
    main()
