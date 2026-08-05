import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# =====================================================
# File paths for LIHC recurrence risk scores
# =====================================================

white_file = (
    "/Users/arnavjoshi/Desktop/LIHCPredictiveAnalyticsRecurrence/"
    "RiskScores_LIHC_White_DFI_CV.csv"
)

all_file = (
    "/Users/arnavjoshi/Desktop/LIHCPredictiveAnalyticsRecurrence/"
    "RiskScores_LIHC_All_DFI_CV.csv"
)

# =====================================================
# Load data
# =====================================================

white_df = pd.read_csv(white_file)
all_df = pd.read_csv(all_file)

# Add group labels
white_df["Group"] = "White"
all_df["Group"] = "All"

# Combine for plotting
combined_df = pd.concat([white_df, all_df], ignore_index=True)

# =====================================================
# Violin plot
# =====================================================

plt.figure(figsize=(10, 6))
sns.violinplot(
    x="Group",
    y="Risk Score",
    data=combined_df,
    inner="quartile",
    palette=["#1f77b4", "#2ca02c"]
)

plt.title("LIHC Recurrence(DFI) Risk Score Distribution(Seventy/Thirty Split Validation)")
plt.xlabel("Group")
plt.ylabel("Risk Score")
plt.grid(True, linestyle="--", alpha=0.6)

# Save instead of showing
plt.tight_layout()
plt.savefig(
    "/Users/arnavjoshi/Desktop/LIHCPredictiveAnalyticsRecurrence/"
    "RiskScore_Comparison_LIHC_DFI.png",
    dpi=300
)
plt.close()

# =====================================================
# Descriptive statistics
# =====================================================

print(
    combined_df
    .groupby("Group")["Risk Score"]
    .describe()
)
