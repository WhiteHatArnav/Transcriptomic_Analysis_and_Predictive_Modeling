import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import cdist
from sklearn.cluster import KMeans

# ------------------
# Paths
# ------------------
INPUT_PATH = "/Users/arnavjoshi/Desktop/BRCAPredictiveAnalytics/brca_black_survival_data.csv"
OUTPUT_DIR = "/Users/arnavjoshi/Desktop/BlackSamplesPCA"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ------------------
# Load data
# ------------------
df = pd.read_csv(INPUT_PATH)

# Required columns
required_cols = ["sample", "OS", "OS.time"]
for col in required_cols:
    if col not in df.columns:
        raise ValueError(f"Missing required column: {col}")

# Separate metadata and expression matrix
meta_cols = ["sample", "OS", "OS.time"]
meta_df = df[meta_cols].copy()
expr_df = df.drop(columns=meta_cols).copy()

# Drop genes with no variance
expr_std = expr_df.std(axis=0, ddof=1)
expr_df = expr_df.loc[:, expr_std > 0]

# Drop samples with missing data
mask_complete = expr_df.notna().all(axis=1)
expr_df = expr_df.loc[mask_complete]
meta_df = meta_df.loc[mask_complete].reset_index(drop=True)
expr_df = expr_df.reset_index(drop=True)

# ------------------
# Sample and gene counts
# ------------------
n_samples = expr_df.shape[0]
n_genes = expr_df.shape[1]

# ------------------
# Standardize features
# ------------------
scaler = StandardScaler(with_mean=True, with_std=True)
X_scaled = scaler.fit_transform(expr_df.values)

# ------------------
# PCA computation
# ------------------
max_components = min(20, n_samples, n_genes)
pca = PCA(n_components=max_components)
X_pca = pca.fit_transform(X_scaled)

explained_var_ratio = pca.explained_variance_ratio_
cumulative_var_ratio = np.cumsum(explained_var_ratio)

# ------------------
# Prepare coloring based on survival variables
# ------------------
OS_series = meta_df["OS"]
OS_time_series = meta_df["OS.time"]
survival_cut = OS_time_series.median()
short_survival_flag = (OS_time_series <= survival_cut).astype(int)

# ------------------
# Scree plot
# ------------------
plt.figure()
pcs = np.arange(1, max_components + 1)
plt.bar(pcs, explained_var_ratio * 100.0)
plt.xlabel("Principal Component")
plt.ylabel("Variance Explained (%)")
plt.title("Scree Plot: Variance Explained per PC")
plt.xticks(pcs)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "BRCA_scree_plot_variance_per_PC.png"), dpi=300)
plt.close()

# ------------------
# Cumulative variance
# ------------------
plt.figure()
plt.plot(pcs, cumulative_var_ratio * 100.0, marker='o')
plt.xlabel("Number of Principal Components")
plt.ylabel("Cumulative Variance Explained (%)")
plt.title("Cumulative Variance Explained")
plt.ylim(0, 100)
plt.grid(True, linestyle="--", alpha=0.4)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "BRCA_cumulative_variance_explained.png"), dpi=300)
plt.close()

# ------------------
# PC1 vs PC2 colored by OS
# ------------------
plt.figure()
scatter = plt.scatter(
    X_pca[:, 0], X_pca[:, 1],
    c=OS_series,
    cmap="coolwarm",
    alpha=0.8,
    edgecolor="k",
    linewidth=0.5
)
plt.xlabel(f"PC1 ({explained_var_ratio[0]*100:.2f}% var)")
plt.ylabel(f"PC2 ({explained_var_ratio[1]*100:.2f}% var)")
plt.title("PC1 vs PC2 colored by OS event status")
cbar = plt.colorbar(scatter)
cbar.set_label("OS (1=event/death,0=censored)")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "BRCA_pc1_pc2_by_OS.png"), dpi=300)
plt.close()

# ------------------
# PC1 vs PC2 colored by short vs long survival
# ------------------
plt.figure()
scatter2 = plt.scatter(
    X_pca[:, 0], X_pca[:, 1],
    c=short_survival_flag,
    cmap="viridis",
    alpha=0.8,
    edgecolor="k",
    linewidth=0.5
)
plt.xlabel(f"PC1 ({explained_var_ratio[0]*100:.2f}% var)")
plt.ylabel(f"PC2 ({explained_var_ratio[1]*100:.2f}% var)")
plt.title("PC1 vs PC2 colored by Short vs Long Survival")
cbar2 = plt.colorbar(scatter2)
cbar2.set_label(f"Short survival (<= median {survival_cut:.1f})")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "BRCA_pc1_pc2_by_survivalTimeMedian.png"), dpi=300)
plt.close()

# ------------------
# Outlier detection in PCA space
# ------------------
k_for_outliers = min(10, max_components)
X_k = X_pca[:, :k_for_outliers]
center = X_k.mean(axis=0, keepdims=True)
cov = np.cov(X_k, rowvar=False)
ridge = 1e-6 * np.eye(cov.shape[0])
inv_cov = np.linalg.pinv(cov + ridge)
diff = X_k - center
md2 = np.sum(diff.dot(inv_cov) * diff, axis=1)
md = np.sqrt(md2)
md_mean = md.mean()
md_std = md.std(ddof=1)
outlier_threshold = md_mean + 3 * md_std
is_outlier = md > outlier_threshold

# ------------------
# Outlier plot
# ------------------
plt.figure()
plt.scatter(
    X_pca[~is_outlier, 0], X_pca[~is_outlier, 1],
    c="lightgray", alpha=0.6, edgecolor="none", label="Inlier"
)
plt.scatter(
    X_pca[is_outlier, 0], X_pca[is_outlier, 1],
    c="red", alpha=0.9, edgecolor="k", linewidth=0.7, label="Outlier"
)
plt.xlabel(f"PC1 ({explained_var_ratio[0]*100:.2f}% var)")
plt.ylabel(f"PC2 ({explained_var_ratio[1]*100:.2f}% var)")
plt.title("PC1 vs PC2 with Outlier Flagging")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "BRCA_pc1_pc2_outliers.png"), dpi=300)
plt.close()

# ------------------
# Save PCA scores and metadata
# ------------------
scores_df = pd.DataFrame(X_pca, columns=[f"PC{i+1}" for i in range(max_components)])
scores_df = pd.concat([meta_df, scores_df], axis=1)
scores_df["short_survival_flag"] = short_survival_flag.values
scores_df["mahalanobis_distance"] = md
scores_df["is_outlier"] = is_outlier  # <-- FIXED: removed .values

scores_df.to_csv(os.path.join(OUTPUT_DIR, "BRCA_per_sample_PCA_scores.csv"), index=False)

# ------------------
# Cluster analysis for estimating synthetic sample scaling
# ------------------
inertias = []
cluster_range = range(1, 6)
for k in cluster_range:
    km = KMeans(n_clusters=k, n_init=10, random_state=0)
    km.fit(X_k)
    inertias.append(km.inertia_)

drops = []
for i in range(1, len(inertias)):
    prev_inertia = inertias[i-1]
    this_inertia = inertias[i]
    drop = (prev_inertia - this_inertia) / prev_inertia
    drops.append(drop)

if len(drops) == 0:
    best_k = 1
else:
    best_k = 1
    for i, d in enumerate(drops, start=2):
        if d < 0.10:
            best_k = i-1
            break
    else:
        best_k = cluster_range[np.argmax(drops)+1]

km_final = KMeans(n_clusters=best_k, n_init=10, random_state=0).fit(X_k)
cluster_labels = km_final.labels_
scores_df["cluster_label"] = cluster_labels
scores_df.to_csv(os.path.join(OUTPUT_DIR, "BRCA_per_sample_PCA_scores_with_clusters.csv"), index=False)

cluster_counts = scores_df["cluster_label"].value_counts().sort_index().to_dict()
outlier_fraction = is_outlier.mean()

# Oversampling heuristic
oversample_factor = 2.5
if best_k > 1:
    oversample_factor = oversample_factor / best_k
if outlier_fraction > 0.2:
    oversample_factor = oversample_factor * (0.2 / outlier_fraction)
oversample_factor = max(1.0, min(oversample_factor, 3.0))
max_synthetic = int(np.floor((oversample_factor - 1.0) * n_samples))

# ------------------
# Summary report
# ------------------
summary_lines = []
summary_lines.append("PCA / Synthetic Sample Guidance Summary (BRCA Black Cohort)")
summary_lines.append("--------------------------------------------------")
summary_lines.append(f"Number of real Black samples used (after QC): {n_samples}")
summary_lines.append(f"Number of genes (features) used after variance filter: {n_genes}")
summary_lines.append("")
summary_lines.append("Explained variance by first 5 PCs (%):")
for i in range(min(5, max_components)):
    summary_lines.append(f"  PC{i+1}: {explained_var_ratio[i]*100:.2f}%")
summary_lines.append(f"Cumulative variance (first {max_components} PCs): {cumulative_var_ratio[-1]*100:.2f}%")
summary_lines.append("")
summary_lines.append(f"Estimated number of clusters in Black samples (k-means elbow): {best_k}")
summary_lines.append(f"Cluster sizes: {cluster_counts}")
summary_lines.append(f"Outlier fraction (>3 SD Mahalanobis): {outlier_fraction*100:.2f}%")
summary_lines.append("")
summary_lines.append("Heuristic oversampling guidance:")
summary_lines.append(
    f"- We can plausibly generate up to ~{oversample_factor:.2f}x total samples "
    f"relative to our {n_samples} real samples without overwhelming the covariance structure."
)
summary_lines.append(
    f"- That corresponds to adding approximately {max_synthetic} synthetic samples "
    "before the representation may start to distort."
)
summary_lines.append("")
summary_lines.append("Interpretation notes:")
summary_lines.append("* If the PC1/PC2 plots show a single compact cluster, synthetic augmentation near that region is likely valid.")
summary_lines.append("* If multiple sub-clusters exist, oversampling should be done within each cluster separately.")
summary_lines.append("* Outliers should be examined carefully before inclusion in data augmentation.")
summary_lines.append("* After augmentation, rerunning LASSO-Cox and verifying convergence and stability is essential.")

with open(os.path.join(OUTPUT_DIR, "BRCA_PCA_summary.txt"), "w") as f:
    f.write("\n".join(summary_lines))

print("Done. All BRCA-prefixed outputs saved to:", OUTPUT_DIR)
