import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# =====================================================
# INPUT FILES (K-FOLD RECURRENCE)
# =====================================================

base_dir = "/Users/arnavjoshi/Desktop/KFoldLassoCox/LUADKFoldLassoCoxRecurrence"

white_file = os.path.join(
    base_dir,
    "RiskScores_LUAD_White_DFI_KFold10.csv"
)

all_file = os.path.join(
    base_dir,
    "RiskScores_LUAD_All_DFI_KFold10.csv"
)

# =====================================================
# LOAD DATA
# =====================================================

white_df = pd.read_csv(white_file)
all_df   = pd.read_csv(all_file)

white_df["Group"] = "White"
all_df["Group"]   = "All"

combined_df = pd.concat(
    [white_df, all_df],
    ignore_index=True
)

# =====================================================
# VIOLIN PLOT
# =====================================================

plt.figure(figsize=(10, 6))

sns.violinplot(
    x="Group",
    y="Risk Score",
    data=combined_df,
    inner="quartile",
    cut=0
)

plt.title("LUAD Recurrence Risk Score Distribution (K-Fold)")
plt.xlabel("Race Group")
plt.ylabel("Risk Score")

plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

# =====================================================
# SAVE
# =====================================================

out_path = os.path.join(
    base_dir,
    "LUAD_Recurrence_KFold_RiskScore_Violin_White_vs_All.png"
)

plt.savefig(out_path, dpi=300)
plt.close()

# =====================================================
# DESCRIPTIVE STATISTICS
# =====================================================

print(
    combined_df
    .groupby("Group")["Risk Score"]
    .describe()
)
