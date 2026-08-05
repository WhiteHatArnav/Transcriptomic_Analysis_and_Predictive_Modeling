import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# -----------------------------
# Paths
# -----------------------------
all_file = "/Users/arnavjoshi/Desktop/KFoldLassoCox/BRCAKFoldLassoCox/BRCA_All_KFold10_RiskScores.csv"
white_file = "/Users/arnavjoshi/Desktop/KFoldLassoCox/BRCAKFoldLassoCox/BRCA_White_KFold10_RiskScores.csv"
black_file = "/Users/arnavjoshi/Desktop/KFoldLassoCox/BRCAKFoldLassoCox/BRCA_Black_KFold10_RiskScores_WithoutOutlier.csv"

output_folder = "/Users/arnavjoshi/Desktop/KFoldLassoCox/BRCAKFoldLassoCox"
os.makedirs(output_folder, exist_ok=True)

# -----------------------------
# Load data
# -----------------------------
all_df = pd.read_csv(all_file)
white_df = pd.read_csv(white_file)
black_df = pd.read_csv(black_file)

# -----------------------------
# Prepare data
# -----------------------------
# Copy RiskScore_OutOfFold into a common column called "Risk Score"
all_df = all_df.copy()
white_df = white_df.copy()
black_df = black_df.copy()

all_df["Risk Score"] = all_df["RiskScore_OutOfFold"]
white_df["Risk Score"] = white_df["RiskScore_OutOfFold"]
black_df["Risk Score"] = black_df["RiskScore_OutOfFold"]

# Add group labels
all_df["Group"] = "All"
white_df["Group"] = "White"
black_df["Group"] = "Black"

# Combine into one dataframe
combined_df = pd.concat([all_df[["Group", "Risk Score"]],
                          white_df[["Group", "Risk Score"]],
                          black_df[["Group", "Risk Score"]]],
                         ignore_index=True)

# -----------------------------
# Violin plot
# -----------------------------
plt.figure(figsize=(10, 6))
sns.violinplot(x="Group", y="Risk Score", data=combined_df, inner="quartile")

plt.title("BRCA Risk Score Distribution: All vs White vs Black")
plt.xlabel("Group")        
plt.ylabel("Risk Score")   
plt.grid(True, linestyle="--", alpha=0.6)

plt.tight_layout()

plot_path = os.path.join(output_folder, "BRCA_RiskScore_ViolinPlot_All_White_Black_KFold_WithoutOutlier.png")
plt.savefig(plot_path, dpi=300)
plt.close()

# -----------------------------
# Summary stats
# -----------------------------
print(combined_df.groupby("Group")["Risk Score"].describe())
print(f"Plot saved to: {plot_path}")
