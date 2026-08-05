import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# File paths for All, White, and Black Race
all_file = "/Users/arnavjoshi/Desktop/LUADPredictiveAnalytics/RiskScores_All_CV.csv"
white_file = "/Users/arnavjoshi/Desktop/LUADPredictiveAnalytics/RiskScores_White_CV.csv"
black_file = "/Users/arnavjoshi/Desktop/LUADPredictiveAnalytics/RiskScores_Black_CV.csv"

# Load data
all_df = pd.read_csv(all_file)
white_df = pd.read_csv(white_file)
black_df = pd.read_csv(black_file)

# Add race group labels
all_df["Group"] = "All"
white_df["Group"] = "White"
black_df["Group"] = "Black"

# Combine for plotting
combined_df = pd.concat([all_df, white_df, black_df])

# Create violin plot
plt.figure(figsize=(12, 6))
sns.violinplot(x="Group", y="Risk Score", data=combined_df, inner="quartile", palette=["#1f77b4", "#ff7f0e", "#2ca02c"])
plt.title("LUAD Risk Score Distribution(Seventy/Thirty Split Validation): All vs White vs Black")
plt.xlabel("Group")
plt.ylabel("Risk Score")
plt.grid(True, linestyle="--", alpha=0.6)

# Save instead of showing
plt.tight_layout()
plt.savefig("/Users/arnavjoshi/Desktop/LUADPredictiveAnalytics/RiskScore_Comparison_By_Race.png", dpi=300)
plt.close()

# Print descriptive statistics
print(combined_df.groupby("Group")["Risk Score"].describe())
