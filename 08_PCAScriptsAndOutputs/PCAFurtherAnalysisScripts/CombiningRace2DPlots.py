import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -----------------------------
# Inputs: race-specific PCA score folders (produced earlier)
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
# Helper
# -----------------------------
def safe_load_scores(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing PCA scores file: {path}")
    df = pd.read_csv(path)
    for col in ["PC1", "PC2", "OS"]:
        if col not in df.columns:
            raise ValueError(f"{col} missing in {path}")
    # Normalize OS to {0,1} if possible
    try:
        df["OS"] = pd.to_numeric(df["OS"], errors="coerce")
    except Exception:
        pass
    df = df.dropna(subset=["OS", "PC1", "PC2"]).copy()
    df["OS"] = df["OS"].astype(int)
    return df

# -----------------------------
# Main
# -----------------------------
def main():
    for cancer, d in RACE_PCA_DIRS.items():
        os.makedirs(d["outfile_dir"], exist_ok=True)

        black_path = os.path.join(d["black"], FILENAME.format(cancer=cancer))
        white_path = os.path.join(d["white"], FILENAME.format(cancer=cancer))

        black_df = safe_load_scores(black_path)
        white_df = safe_load_scores(white_path)

        # Masks for four groups
        wb1 = (white_df["OS"] == 1)
        wb0 = (white_df["OS"] == 0)
        bb1 = (black_df["OS"] == 1)
        bb0 = (black_df["OS"] == 0)

        plt.figure(figsize=(7.6, 6.4))

        # Colors by OS, markers by Race
        # OS=1 -> red; OS=0 -> blue
        # White -> 'o'; Black -> '^'
        # Plot in an order that keeps legend clean and points visible
        # White, OS=0
        if wb0.any():
            plt.scatter(
                white_df.loc[wb0, "PC1"], white_df.loc[wb0, "PC2"],
                c="#1f77b4", alpha=0.85, edgecolor="k", linewidth=0.3,
                marker="o", s=45, label="White, OS=0"
            )
        # White, OS=1
        if wb1.any():
            plt.scatter(
                white_df.loc[wb1, "PC1"], white_df.loc[wb1, "PC2"],
                c="#d62728", alpha=0.95, edgecolor="k", linewidth=0.3,
                marker="o", s=55, label="White, OS=1"
            )
        # Black, OS=0
        if bb0.any():
            plt.scatter(
                black_df.loc[bb0, "PC1"], black_df.loc[bb0, "PC2"],
                c="#1f77b4", alpha=0.95, edgecolor="k", linewidth=0.3,
                marker="^", s=55, label="Black, OS=0"
            )
        # Black, OS=1
        if bb1.any():
            plt.scatter(
                black_df.loc[bb1, "PC1"], black_df.loc[bb1, "PC2"],
                c="#d62728", alpha=1.0, edgecolor="k", linewidth=0.3,
                marker="^", s=65, label="Black, OS=1"
            )

        plt.xlabel("PC1")
        plt.ylabel("PC2")
        plt.title(f"{cancer}: PC1 vs PC2 — White vs Black by OS")
        lgd = plt.legend(frameon=True, loc="best")
        for lh in lgd.legendHandles:
            try:
                lh.set_sizes([60.0])
            except Exception:
                pass
        plt.tight_layout()

        outpath = os.path.join(d["outfile_dir"], f"{cancer}_WhiteBlack_PC1_PC2_combined.png")
        plt.savefig(outpath, dpi=300)
        plt.close()
        print(f"Saved 2D overlay (4 groups): {outpath}")

if __name__ == "__main__":
    main()
