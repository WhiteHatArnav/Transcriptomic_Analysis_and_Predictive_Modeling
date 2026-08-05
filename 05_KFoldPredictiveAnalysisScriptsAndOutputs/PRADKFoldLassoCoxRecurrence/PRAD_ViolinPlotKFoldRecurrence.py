import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ----------------------------------------
# File paths for PRAD K-Fold Recurrence Risk Scores
# ----------------------------------------
base_folder = "/Users/arnavjoshi/Desktop/KFoldLassoCox/PRADKFoldLassoCoxRecurrence"

# All races: NON-censoring recurrence K-Fold file
all_file = f"{base_folder}/PRAD_All_KFold10_RiskScores_Recurrence.csv"

# White race: recurrence K-Fold file (no censoring)
white_file = f"{base_folder}/PRAD_White_KFold10_RiskScores_Recurrence.csv"

# Black race: recurrence K-Fold file WITH OS.time censoring
black_file = f"{base_folder}/PRAD_Black_KFold10_RiskScores_Recurrence_OSCensoring.csv"

# Load data
all_df = pd.read_csv(all_file)
white_df = pd.read_csv(white_file)
black_df = pd.read_csv(black_file)

# Add group labels
all_df["Group"] = "All"
white_df["Group"] = "White"
black_df["Group"] = "Black"

# Keep only relevant columns and standardize column name for plotting
all_df = all_df[["Patient ID", "RiskScore_OutOfFold", "Group"]].rename(
    columns={"RiskScore_OutOfFold": "Risk Score"}
)
white_df = white_df[["Patient ID", "RiskScore_OutOfFold", "Group"]].rename(
    columns={"RiskScore_OutOfFold": "Risk Score"}
)
black_df = black_df[["Patient ID", "RiskScore_OutOfFold", "Group"]].rename(
    columns={"RiskScore_OutOfFold": "Risk Score"}
)

# Combine for plotting
combined_df = pd.concat([all_df, white_df, black_df], ignore_index=True)

# ----------------------------------------
# Violin plot of K-Fold recurrence risk scores
# ----------------------------------------
plt.figure(figsize=(10, 6))
sns.violinplot(
    x="Group",
    y="Risk Score",
    data=combined_df,
    inner="quartile"
)

plt.title("PRAD Recurrence(BCR) Risk Score Distribution by Race(10-Fold Cross Validation)")
plt.xlabel("Group")
plt.ylabel("Recurrence Risk Score (Out-of-Fold)")

plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

# Save plot
output_plot = f"{base_folder}/KFold_Recurrence_RiskScore_Comparison_All_White_Black.png"
plt.savefig(output_plot, dpi=300)
plt.close()

print(f"Plot saved to:\n{output_plot}\n")

# ----------------------------------------
# Statistical summary
# ----------------------------------------
print("Statistical Summary (K-Fold Recurrence Risk Scores):")
print(combined_df.groupby("Group")["Risk Score"].describe())
