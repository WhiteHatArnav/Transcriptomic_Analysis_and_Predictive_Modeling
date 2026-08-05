import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# =====================================================
# PATHS
# =====================================================

output_dir = "/Users/arnavjoshi/Desktop/KFoldLassoCox/BRCAKFoldLassoCoxRecurrence/"
os.makedirs(output_dir, exist_ok=True)

risk_files = {
    "White": "/Users/arnavjoshi/Desktop/KFoldLassoCox/BRCAKFoldLassoCoxRecurrence/RiskScores_BRCA_White_DFI_KFold10.csv",
    "All":   "/Users/arnavjoshi/Desktop/KFoldLassoCox/BRCAKFoldLassoCoxRecurrence/RiskScores_BRCA_All_DFI_KFold10.csv",
    "Black": "/Users/arnavjoshi/Desktop/KFoldLassoCox/BRCAKFoldLassoCoxRecurrence/RiskScores_BRCA_Black_DFI_KFold10.csv",
}

# =====================================================
# LOAD + COMBINE
# =====================================================

dfs = []

for group, path in risk_files.items():
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()
    df["Group"] = group
    dfs.append(df)

combined_df = pd.concat(dfs, ignore_index=True)

# =====================================================
# VIOLIN PLOT
# =====================================================

plt.figure(figsize=(12, 6))
sns.violinplot(
    x="Group",
    y="Risk Score",
    data=combined_df,
    inner="quartile",
    cut=0,
    palette=["#1f77b4", "#2ca02c", "#d62728"]
)

plt.title("BRCA Recurrence(DFI) Risk Score Distribution by Race(10-Fold Cross Validation)")
plt.xlabel("Race Group")
plt.ylabel("Out-of-Fold Risk Score")
plt.grid(True, linestyle="--", alpha=0.6)

plt.tight_layout()
plt.savefig(
    os.path.join(output_dir, "BRCA_Recurrence_RiskScore_Violin_KFold.png"),
    dpi=300
)
plt.close()

# =====================================================
# DESCRIPTIVE STATS (OPTIONAL BUT USEFUL)
# =====================================================

print(
    combined_df
    .groupby("Group")["Risk Score"]
    .describe()
)
