import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# =====================================================
# File paths for BRCA recurrence risk scores
# =====================================================

all_file = (
    "/Users/arnavjoshi/Desktop/BRCAPredictiveAnalyticsRecurrence/"
    "RiskScores_BRCA_All_DFI_CV.csv"
)

white_file = (
    "/Users/arnavjoshi/Desktop/BRCAPredictiveAnalyticsRecurrence/"
    "RiskScores_BRCA_White_DFI_CV.csv"
)

black_file = (
    "/Users/arnavjoshi/Desktop/BRCAPredictiveAnalyticsRecurrence/"
    "RiskScores_BRCA_Black_DFI_CV.csv"
)

# =====================================================
# Load data
# =====================================================

all_df = pd.read_csv(all_file)
white_df = pd.read_csv(white_file)
black_df = pd.read_csv(black_file)

# Add race group labels
all_df["Group"] = "All"
white_df["Group"] = "White"
black_df["Group"] = "Black"

# Combine for plotting
combined_df = pd.concat([all_df, white_df, black_df], ignore_index=True)

# =====================================================
# Violin plot
# =====================================================

plt.figure(figsize=(12, 6))
sns.violinplot(
    x="Group",
    y="Risk Score",
    data=combined_df,
    inner="quartile",
    palette=["#1f77b4", "#ff7f0e", "#2ca02c"]
)

plt.title("BRCA Recurrence(DFI) Risk Score Distribution by Race(Seventy/Thirty Split Validation)")
plt.xlabel("Group")
plt.ylabel("Risk Score")
plt.grid(True, linestyle="--", alpha=0.6)

# Save instead of showing
plt.tight_layout()
plt.savefig(
    "/Users/arnavjoshi/Desktop/BRCAPredictiveAnalyticsRecurrence/"
    "RiskScore_Comparison_By_Race_BRCA_DFI.png",
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
