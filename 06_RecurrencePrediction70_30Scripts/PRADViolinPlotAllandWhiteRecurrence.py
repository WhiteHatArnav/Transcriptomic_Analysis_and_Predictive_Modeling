import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ----------------------------------------
# File paths (recurrence-based 70/30 scores)
# ----------------------------------------
all_file = "/Users/arnavjoshi/Desktop/PRADPredictiveAnalyticsRecurrence/RiskScores_All_Recurrence_CV.csv"
white_file = "/Users/arnavjoshi/Desktop/PRADPredictiveAnalyticsRecurrence/RiskScores_White_Recurrence_CV.csv"

# Load data
all_df = pd.read_csv(all_file)
white_df = pd.read_csv(white_file)

# Add race labels
all_df["Group"] = "All"
white_df["Group"] = "White"

# Combine for plotting
combined_df = pd.concat([all_df, white_df], ignore_index=True)

# ----------------------------------------
# Violin plot of recurrence risk scores
# ----------------------------------------
plt.figure(figsize=(10, 6))
sns.violinplot(
    x="Group",
    y="Risk Score",
    data=combined_df,
    inner="quartile"
)

plt.title("PRAD Recurrence(BCR) Risk Score Distribution: All vs White(Seventy/Thirty Split Validation)")
plt.xlabel("Group")
plt.ylabel("Recurrence Risk Score")

plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

# Save output
output_plot = "/Users/arnavjoshi/Desktop/PRADPredictiveAnalyticsRecurrence/Recurrence_RiskScore_Comparison_All_vs_White.png"
plt.savefig(output_plot, dpi=300)
plt.close()

print(f"Plot saved to:\n{output_plot}\n")

# ----------------------------------------
# Print statistical summary
# ----------------------------------------
print("Statistical Summary:")
print(combined_df.groupby("Group")["Risk Score"].describe())
