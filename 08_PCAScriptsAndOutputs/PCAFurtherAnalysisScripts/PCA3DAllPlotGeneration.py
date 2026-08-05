# This single script:
#
# 1. Overlays 3D plots from separate PCAs (reads each race’s PCA scores).
# 2. Performs a single combined PCA across both races and plots the 4 OS subgroups.
# 3. Generates separate 3D plots for each race (saved to their original race-specific PCA folders).
#
# All output folders and filenames follow the existing structure and conventions.


import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# -----------------------------
# Folder Definitions
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

FILENAME = "{cancer}_per_sample_PCA_scores.csv"

# -----------------------------
# Utility Functions
# -----------------------------
def safe_load_scores(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing file: {path}")
    df = pd.read_csv(path)
    for col in ["PC1", "PC2", "PC3", "OS"]:
        if col not in df.columns:
            raise ValueError(f"{col} missing in {path}")
    df["OS"] = pd.to_numeric(df["OS"], errors="coerce").fillna(0).astype(int)
    return df

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
# Plot Helpers
# -----------------------------
def plot_separate_race(df, cancer, race, outdir):
    """Creates single 3D PCA plot per race from its PCA scores."""
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")

    os0 = (df["OS"] == 0)
    os1 = (df["OS"] == 1)

    if os0.any():
        ax.scatter(df.loc[os0, "PC1"], df.loc[os0, "PC2"], df.loc[os0, "PC3"],
                   c="#1f77b4", alpha=0.9, edgecolor="k", linewidth=0.3,
                   s=40, marker="o", label=f"{race}, OS=0")
    if os1.any():
        ax.scatter(df.loc[os1, "PC1"], df.loc[os1, "PC2"], df.loc[os1, "PC3"],
                   c="#d62728", alpha=1.0, edgecolor="k", linewidth=0.3,
                   s=55, marker="^", label=f"{race}, OS=1")

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_zlabel("PC3")
    ax.set_title(f"{cancer} {race} PCA (3D)")
    ax.legend()
    plt.tight_layout()

    outpath = os.path.join(outdir, f"{cancer}_{race}_3D_PCA.png")
    plt.savefig(outpath, dpi=300)
    plt.close()
    print(f"Saved separate 3D PCA for {cancer} {race}: {outpath}")

def plot_combined_overlay(white_df, black_df, cancer, outfile_dir):
    """Overlays existing PCA plots (from separate PCAs)."""
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")

    # 4 subgroups: Race x OS
    masks = {
        "White, OS=0": (white_df["OS"] == 0),
        "White, OS=1": (white_df["OS"] == 1),
        "Black, OS=0": (black_df["OS"] == 0),
        "Black, OS=1": (black_df["OS"] == 1),
    }

    colors = {
        "White, OS=0": "#1f77b4",
        "White, OS=1": "#d62728",
        "Black, OS=0": "#1f77b4",
        "Black, OS=1": "#d62728",
    }
    markers = {
        "White, OS=0": "o",
        "White, OS=1": "o",
        "Black, OS=0": "^",
        "Black, OS=1": "^",
    }

    for label, mask in masks.items():
        df = white_df if "White" in label else black_df
        if mask.any():
            ax.scatter(df.loc[mask, "PC1"], df.loc[mask, "PC2"], df.loc[mask, "PC3"],
                       c=colors[label], marker=markers[label],
                       alpha=0.9, edgecolor="k", linewidth=0.3, s=45, label=label)

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_zlabel("PC3")
    ax.set_title(f"{cancer}: Combined Plot (Separate PCAs Overlay)")
    ax.legend()
    plt.tight_layout()

    outpath = os.path.join(outfile_dir, f"{cancer}_WhiteBlack_SeparatePCA_3D.png")
    plt.savefig(outpath, dpi=300)
    plt.close()
    print(f"Saved combined overlay (Separate PCA): {outpath}")

def plot_single_combined_pca(black_path, white_path, cancer, outfile_dir):
    """Single PCA across both races, 3D plot (4 groups)."""
    os.makedirs(outfile_dir, exist_ok=True)
    black_df = load_expression_matrix(black_path)
    white_df = load_expression_matrix(white_path)

    gene_cols = [c for c in black_df.columns if c not in ["Patient ID", "OS", "OS.time"]]
    gene_cols = [c for c in gene_cols if c in white_df.columns]

    black_df["Race"] = "Black"
    white_df["Race"] = "White"

    combined_df = pd.concat([black_df, white_df], ignore_index=True)
    X = combined_df[gene_cols].values
    X_scaled = StandardScaler().fit_transform(X)
    pca = PCA(n_components=3, random_state=42)
    pcs = pca.fit_transform(X_scaled)

    combined_df["PC1"] = pcs[:, 0]
    combined_df["PC2"] = pcs[:, 1]
    combined_df["PC3"] = pcs[:, 2]

    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")

    groups = [
        ("White, OS=0", (combined_df["Race"] == "White") & (combined_df["OS"] == 0)),
        ("White, OS=1", (combined_df["Race"] == "White") & (combined_df["OS"] == 1)),
        ("Black, OS=0", (combined_df["Race"] == "Black") & (combined_df["OS"] == 0)),
        ("Black, OS=1", (combined_df["Race"] == "Black") & (combined_df["OS"] == 1)),
    ]

    for label, mask in groups:
        if mask.any():
            ax.scatter(combined_df.loc[mask, "PC1"], combined_df.loc[mask, "PC2"], combined_df.loc[mask, "PC3"],
                       c="#1f77b4" if "OS=0" in label else "#d62728",
                       marker="o" if "White" in label else "^",
                       alpha=0.9, edgecolor="k", linewidth=0.3, s=50, label=label)

    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}% var)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}% var)")
    ax.set_zlabel(f"PC3 ({pca.explained_variance_ratio_[2]*100:.2f}% var)")
    ax.set_title(f"{cancer}: Single PCA (Combined 3D)")
    ax.legend()
    plt.tight_layout()

    outpath = os.path.join(outfile_dir, f"{cancer}_WhiteBlack_singlePCA_3D.png")
    plt.savefig(outpath, dpi=300)
    plt.close()
    print(f"Saved combined single PCA 3D: {outpath}")

# -----------------------------
# Main Routine
# -----------------------------
def main():
    for cancer, d in RACE_PCA_DIRS.items():
        os.makedirs(d["outfile_dir"], exist_ok=True)
        black_path = os.path.join(d["black"], FILENAME.format(cancer=cancer))
        white_path = os.path.join(d["white"], FILENAME.format(cancer=cancer))

        # Separate 3D plots per race
        black_df = safe_load_scores(black_path)
        white_df = safe_load_scores(white_path)
        plot_separate_race(black_df, cancer, "Black", d["black"])
        plot_separate_race(white_df, cancer, "White", d["white"])

        # Overlay of separate PCA results
        plot_combined_overlay(white_df, black_df, cancer, d["outfile_dir"])

        # Single PCA (joint)
        plot_single_combined_pca(black_path, white_path, cancer, d["outfile_dir"])

if __name__ == "__main__":
    main()
