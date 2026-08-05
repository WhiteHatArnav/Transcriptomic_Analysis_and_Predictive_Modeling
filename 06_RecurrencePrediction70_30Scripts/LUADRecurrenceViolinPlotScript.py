import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# =====================================================
# File paths for LUAD recurrence risk scores (70:30)
# =====================================================

white_file = (
    "/Users/arnavjoshi/Desktop/LUADPredictiveAnalyticsRecurrence/"
    "RiskScores_LUAD_White_DFI_CV_7030.csv"
)

black_file = (
    "/Users/arnavjoshi/Desktop/LUADPredictiveAnalyticsRecurrence/"
    "RiskScores_LUAD_Black_DFI_CV_7030.csv"
)

all_file = (
    "/Users/arnavjoshi/Desktop/LUADPredictiveAnalyticsRecurrence/"
    "RiskScores_LUAD_All_DFI_CV_7030.csv"
)

# =====================================================
# Load data
# =====================================================

white_df = pd.read_csv(white_file)
black_df = pd.read_csv(black_file)
all_df = pd.read_csv(all_file)

# Add group labels
white_df["Group"] = "White"
black_df["Group"] = "Black"
all_df["Group"] = "All"

# Combine for plotting
combined_df = pd.concat(
    [white_df, black_df, all_df],
    ignore_index=True
)

# =====================================================
# Violin plot
# =====================================================

plt.figure(figsize=(12, 6))
sns.violinplot(
    x="Group",
    y="Risk Score",
    data=combined_df,
    inner="quartile",
    palette=["#1f77b4", "#d62728", "#2ca02c"]
)

plt.title("LUAD Recurrence(DFI) Risk Score Distribution by Race(Seventy/Thirty Split Cross Validation)")
plt.xlabel("Race Group")
plt.ylabel("Risk Score")
plt.grid(True, linestyle="--", alpha=0.6)

# Save instead of showing
plt.tight_layout()
plt.savefig(
    "/Users/arnavjoshi/Desktop/LUADPredictiveAnalyticsRecurrence/"
    "RiskScore_Comparison_By_Race_LUAD_DFI_7030.png",
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
