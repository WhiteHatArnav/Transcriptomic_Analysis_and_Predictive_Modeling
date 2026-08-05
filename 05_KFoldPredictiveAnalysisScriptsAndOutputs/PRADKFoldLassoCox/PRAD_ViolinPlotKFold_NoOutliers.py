import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# -----------------------------
# Paths (PRAD K-Fold outputs)
# -----------------------------
all_file = "/Users/arnavjoshi/Desktop/KFoldLassoCox/PRADKFoldLassoCox/PRAD_All_KFold10_RiskScores.csv"
white_file = "/Users/arnavjoshi/Desktop/KFoldLassoCox/PRADKFoldLassoCox/PRAD_White_KFold10_RiskScores.csv"
black_file = "/Users/arnavjoshi/Desktop/KFoldLassoCox/PRADKFoldLassoCox/PRAD_Black_KFold10_RiskScores.csv"

output_folder = "/Users/arnavjoshi/Desktop/KFoldLassoCox/PRADKFoldLassoCox"
os.makedirs(output_folder, exist_ok=True)

# -----------------------------
# Load & prepare
# -----------------------------
all_df   = pd.read_csv(all_file).copy()
white_df = pd.read_csv(white_file).copy()
black_df = pd.read_csv(black_file).copy()

# Harmonize column name to "Risk Score"
for df in (all_df, white_df, black_df):
    if "Risk Score" not in df.columns and "RiskScore_OutOfFold" in df.columns:
        df["Risk Score"] = df["RiskScore_OutOfFold"]

all_df["Group"]   = "All"
white_df["Group"] = "White"
black_df["Group"] = "Black"

combined = pd.concat(
    [all_df[["Group", "Patient ID", "Risk Score"]],
     white_df[["Group", "Patient ID", "Risk Score"]],
     black_df[["Group", "Patient ID", "Risk Score"]]],
    ignore_index=True
)

# Ensure strictly positive (partial hazards should be >0; if any zeros, bump by tiny eps)
min_pos = combined.loc[combined["Risk Score"] > 0, "Risk Score"].min()
eps = min_pos * 1e-6 if pd.notnull(min_pos) else 1e-12
combined["Risk Score"] = combined["Risk Score"].clip(lower=eps)

# -----------------------------
# Robust MAD trimming on log10 scale (per group)
# -----------------------------
def trim_outliers_mad_log(df, value_col="Risk Score", group_col="Group", id_col="Patient ID", threshold=3.5):
    """
    Remove points with robust |z| > threshold, where z computed on log10(value).
    z_robust = 0.6745 * (log10(x) - median_log) / MAD_log
    """
    kept = []
    removed_report = []

    for g, sub in df.groupby(group_col, sort=False):
        vals = sub[value_col].astype(float)
        logv = np.log10(vals)

        med = np.median(logv)
        mad = np.median(np.abs(logv - med))

        if mad == 0 or np.isnan(mad):
            # No dispersion or too few points — keep as is
            kept.append(sub)
            continue

        z = 0.6745 * (logv - med) / mad
        mask_keep = np.abs(z) <= threshold

        kept.append(sub[mask_keep])
        removed = sub[~mask_keep]
        if len(removed) > 0:
            removed_report.append({
                "Group": g,
                "Removed_N": int(len(removed)),
                "IDs": removed[id_col].tolist() if id_col in removed.columns else [],
                "Values": removed[value_col].tolist()
            })

    trimmed = pd.concat(kept, ignore_index=True)

    # Optional: print a brief report of removals
    if removed_report:
        print("\nRemoved outliers (log10-MAD) by group:")
        for r in removed_report:
            print(f"- {r['Group']}: removed {r['Removed_N']} outlier(s)")
    else:
        print("\nNo outliers removed by log10-MAD rule.")

    return trimmed

# Use threshold=3.5 (common robust rule). Tighten to 3.0 if tails remain too fat.
combined_trim = trim_outliers_mad_log(combined, threshold=3.5)

# -----------------------------
# Violin plot (trimmed)
# -----------------------------
plt.figure(figsize=(10, 6))
sns.violinplot(
    x="Group",
    y="Risk Score",
    data=combined_trim,
    inner="quartile",
    palette="Set2",
    cut=0
)
plt.title("PRAD Risk Score Distribution (Outliers removed via log10-MAD)")
plt.xlabel("Group")
plt.ylabel("Risk Score")
plt.grid(True, linestyle="--", alpha=0.6)

out_png = os.path.join(output_folder, "PRAD_RiskScore_Violin_NoOutliers.png")
plt.tight_layout()
plt.savefig(out_png, dpi=300)
plt.close()
print(f"\nSaved violin (trimmed) to: {out_png}")

# -----------------------------
# (Optional) Also show the same plot on log10 scale for diagnostics
# -----------------------------
combined_trim = combined_trim.assign(Log10Risk=np.log10(combined_trim["Risk Score"]))
plt.figure(figsize=(10, 6))
sns.violinplot(
    x="Group",
    y="Log10Risk",
    data=combined_trim,
    inner="quartile",
    palette="Set2",
    cut=0
)
plt.title("PRAD log10(Risk Score) Distribution (Outliers removed via log10-MAD)")
plt.xlabel("Group")
plt.ylabel("log10(Risk Score)")
plt.grid(True, linestyle="--", alpha=0.6)

out_png_log = os.path.join(output_folder, "PRAD_RiskScore_Violin_Log10_NoOutliers_logMAD.png")
plt.tight_layout()
plt.savefig(out_png_log, dpi=300)
plt.close()
print(f"Saved violin (log10, trimmed) to: {out_png_log}")

# -----------------------------
# Summary stats post-trim
# -----------------------------
print("\nSummary after trimming (original scale):")
print(combined_trim.groupby("Group")["Risk Score"].describe())

print("\nSummary after trimming (log10 scale):")
print(combined_trim.groupby("Group")["Log10Risk"].describe())
