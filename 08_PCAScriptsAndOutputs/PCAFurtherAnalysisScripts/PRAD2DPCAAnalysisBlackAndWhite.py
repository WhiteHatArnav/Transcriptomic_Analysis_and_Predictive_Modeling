import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os

# ------------------
# Paths
# ------------------
INPUTS = {
    "Black": "/Users/arnavjoshi/Desktop/PRADPredictiveAnalytics/prad_black_survival_data.csv",
    "White": "/Users/arnavjoshi/Desktop/PRADPredictiveAnalytics/prad_white_survival_data.csv",
}
OUTPUTS = {
    "Black": "/Users/arnavjoshi/Desktop/BlackSamplesPCA/PRAD_PCA",
    "White": "/Users/arnavjoshi/Desktop/WhiteSamplesPCA/PRAD_PCA",
}

def run_one(race_label, input_path, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(input_path)
    # Ensure first col is a sample id
    if df.columns[0] != "sample" and df.columns[0] != "Patient ID":
        df = df.rename(columns={df.columns[0]: "sample"})
    if "sample" not in df.columns:
        df = df.rename(columns={df.columns[0]: "sample"})
    # Required columns
    for col in ["OS", "OS.time"]:
        if col not in df.columns:
            raise ValueError(f"Missing required column {col} in {input_path}")

    meta_cols = ["sample", "OS", "OS.time"]
    meta_df = df[meta_cols].copy()
    expr_df = df.drop(columns=meta_cols).copy()

    # Drop zero-variance genes
    expr_std = expr_df.std(axis=0, ddof=1)
    expr_df = expr_df.loc[:, expr_std > 0]

    # Drop samples with any missing expression
    mask_complete = expr_df.notna().all(axis=1)
    expr_df = expr_df.loc[mask_complete].reset_index(drop=True)
    meta_df = meta_df.loc[mask_complete].reset_index(drop=True)

    n_samples, n_genes = expr_df.shape
    if n_samples < 3:
        raise ValueError(f"Too few samples for PCA in {race_label}: {n_samples}")

    # Standardize
    X = expr_df.values
    X_scaled = StandardScaler(with_mean=True, with_std=True).fit_transform(X)

    # We want at least 3 PCs for 3D plots; cap at 20 or rank limit
    max_components = min(20, n_samples, n_genes)
    n_components = max(3, min(max_components, 20))
    pca = PCA(n_components=n_components, random_state=42)
    X_pca = pca.fit_transform(X_scaled)

    evr = pca.explained_variance_ratio_
    cev = np.cumsum(evr)

    # Survival-based coloring
    OS_series = meta_df["OS"].astype(int)
    OS_time_series = meta_df["OS.time"]
    survival_cut = OS_time_series.median()
    short_survival_flag = (OS_time_series <= survival_cut).astype(int)

    # Scree
    pcs = np.arange(1, n_components + 1)
    plt.figure()
    plt.bar(pcs, evr * 100.0)
    plt.xlabel("Principal Component")
    plt.ylabel("Variance Explained (%)")
    plt.title("Scree Plot: Variance Explained per PC")
    plt.xticks(pcs)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "PRAD_scree_plot_variance_per_PC.png"), dpi=300)
    plt.close()

    # Cumulative
    plt.figure()
    plt.plot(pcs, cev * 100.0, marker='o')
    plt.xlabel("Number of Principal Components")
    plt.ylabel("Cumulative Variance Explained (%)")
    plt.title("Cumulative Variance Explained")
    plt.ylim(0, 100)
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "PRAD_cumulative_variance_explained.png"), dpi=300)
    plt.close()

    # PC1 vs PC2 by OS
    plt.figure()
    scatter = plt.scatter(
        X_pca[:, 0], X_pca[:, 1],
        c=OS_series, cmap="coolwarm",
        alpha=0.85, edgecolor="k", linewidth=0.4
    )
    plt.xlabel(f"PC1 ({evr[0]*100:.2f}% var)")
    plt.ylabel(f"PC2 ({evr[1]*100:.2f}% var)")
    plt.title(f"PRAD {race_label}: PC1 vs PC2 by OS")
    cbar = plt.colorbar(scatter)
    cbar.set_label("OS (1=event, 0=censored)")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "PRAD_pc1_pc2_by_OS.png"), dpi=300)
    plt.close()

    # PC1 vs PC2 by survival time median
    plt.figure()
    scatter2 = plt.scatter(
        X_pca[:, 0], X_pca[:, 1],
        c=short_survival_flag, cmap="viridis",
        alpha=0.85, edgecolor="k", linewidth=0.4
    )
    plt.xlabel(f"PC1 ({evr[0]*100:.2f}% var)")
    plt.ylabel(f"PC2 ({evr[1]*100:.2f}% var)")
    plt.title(f"PRAD {race_label}: PC1 vs PC2 by Survival Time (median={survival_cut:.1f})")
    cbar2 = plt.colorbar(scatter2)
    cbar2.set_label(f"Short survival (<= {survival_cut:.1f})")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "PRAD_pc1_pc2_by_survivalTimeMedian.png"), dpi=300)
    plt.close()

    # Outliers via Mahalanobis on first k PCs
    k_for_outliers = min(10, n_components)
    X_k = X_pca[:, :k_for_outliers]
    center = X_k.mean(axis=0, keepdims=True)
    cov = np.cov(X_k, rowvar=False)
    ridge = 1e-6 * np.eye(cov.shape[0])
    inv_cov = np.linalg.pinv(cov + ridge)
    diff = X_k - center
    md2 = np.sum(diff.dot(inv_cov) * diff, axis=1)
    md = np.sqrt(md2)
    thr = md.mean() + 3 * md.std(ddof=1)
    is_outlier = md > thr

    plt.figure()
    plt.scatter(X_pca[~is_outlier, 0], X_pca[~is_outlier, 1],
                c="lightgray", alpha=0.6, edgecolor="none", label="Inlier")
    plt.scatter(X_pca[is_outlier, 0], X_pca[is_outlier, 1],
                c="red", alpha=0.95, edgecolor="k", linewidth=0.5, label="Outlier")
    plt.xlabel(f"PC1 ({evr[0]*100:.2f}% var)")
    plt.ylabel(f"PC2 ({evr[1]*100:.2f}% var)")
    plt.title(f"PRAD {race_label}: PC1 vs PC2 with Outliers")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "PRAD_pc1_pc2_outliers.png"), dpi=300)
    plt.close()

    # Save scores with PC1..PCk (incl. PC3 for 3D)
    scores_df = pd.DataFrame(X_pca, columns=[f"PC{i+1}" for i in range(n_components)])
    scores_df = pd.concat([meta_df.reset_index(drop=True), scores_df], axis=1)
    scores_df["short_survival_flag"] = short_survival_flag.values
    scores_df["mahalanobis_distance"] = md
    scores_df["is_outlier"] = is_outlier
    scores_df.to_csv(os.path.join(output_dir, "PRAD_per_sample_PCA_scores.csv"), index=False)

    # Summary
    summary = []
    summary.append(f"PCA Summary (PRAD {race_label})")
    summary.append(f"Samples after QC: {n_samples}")
    summary.append(f"Genes after variance filter: {n_genes}")
    summary.append("Explained variance first 5 PCs (%):")
    for i in range(min(5, n_components)):
        summary.append(f"  PC{i+1}: {evr[i]*100:.2f}")
    summary.append(f"Cumulative (first {n_components} PCs): {cev[-1]*100:.2f}%")
    summary.append(f"Outlier fraction: {is_outlier.mean()*100:.2f}%")
    with open(os.path.join(output_dir, "PRAD_PCA_summary.txt"), "w") as f:
        f.write("\n".join(summary))

    print(f"Done PRAD {race_label}: outputs -> {output_dir}")

def main():
    for race in ["Black", "White"]:
        run_one(race, INPUTS[race], OUTPUTS[race])

if __name__ == "__main__":
    main()
