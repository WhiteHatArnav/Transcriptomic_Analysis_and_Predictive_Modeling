import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# File paths for All Race and White Race
all_file = "/Users/arnavjoshi/Desktop/LIHCPredictiveAnalytics/RiskScores_All_CV.csv"
white_file = "/Users/arnavjoshi/Desktop/LIHCPredictiveAnalytics/RiskScores_White_CV.csv"

# Load data
all_df = pd.read_csv(all_file)
white_df = pd.read_csv(white_file)

# Add race group labels
all_df["Group"] = "All"
white_df["Group"] = "White"

# Combine for plotting
combined_df = pd.concat([all_df, white_df])

# Create violin plot
plt.figure(figsize=(10, 6))
sns.violinplot(x="Group", y="Risk Score", data=combined_df, inner="quartile")
plt.title("LIHC Risk Score Distribution(Seventy/Thirty Split Validation): All vs White")
plt.xlabel("Group")
plt.ylabel("Risk Score")
plt.grid(True, linestyle="--", alpha=0.6)

# Save instead of showing
plt.tight_layout()
plt.savefig("/Users/arnavjoshi/Desktop/LIHCPredictiveAnalytics/RiskScore_Comparison_All_vs_White.png", dpi=300)
plt.close()

print(combined_df.groupby("Group")["Risk Score"].describe())
