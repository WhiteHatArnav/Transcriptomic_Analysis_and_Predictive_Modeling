import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# File paths for All, White, and Black
base = "/Users/arnavjoshi/Desktop/BRCAPredictiveAnalytics"
all_file   = f"{base}/RiskScores_All_CV.csv"
white_file = f"{base}/RiskScores_White_CV.csv"
black_file = f"{base}/RiskScores_Black_CV.csv"

# Load data
all_df = pd.read_csv(all_file)
white_df = pd.read_csv(white_file)
black_df = pd.read_csv(black_file)

# Add race group labels
all_df["Group"] = "All"
white_df["Group"] = "White"
black_df["Group"] = "Black"

# Combine for plotting
combined_df = pd.concat([all_df, white_df, black_df], ignore_index=True)

# Create violin plot
plt.figure(figsize=(10, 6))
sns.violinplot(
    x="Group",
    y="Risk Score",
    data=combined_df,
    inner="quartile",
    order=["All", "White", "Black"]
)
plt.title("BRCA Risk Score Distribution(Seventy/Thirty Split Validation): All vs White vs Black")
plt.xlabel("Group")
plt.ylabel("Risk Score")
plt.grid(True, linestyle="--", alpha=0.6)

# Save instead of showing
plt.tight_layout()
plt.savefig(f"{base}/RiskScore_Comparison_All_vs_White_vs_Black.png", dpi=300)
plt.close()

# Summary stats by group
print(combined_df.groupby("Group")["Risk Score"].describe())
