import pandas as pd

# -----------------------------
# Paths
# -----------------------------
input_path = "/Users/arnavjoshi/Desktop/KFoldLassoCox/BRCAKFoldLassoCox/BRCA_White_KFold10_AllSelectedFeatures.csv"
output_path = "/Users/arnavjoshi/Desktop/KFoldLassoCox/BRCAKFoldLassoCox/BRCA_White_KFold10_StableFeatures.csv"

# -----------------------------
# Parameters
# -----------------------------
freq_threshold = 2  # feature must appear in ≥2 folds to be considered stable

# -----------------------------
# Load and validate input
# -----------------------------
df = pd.read_csv(input_path)

if "Feature" not in df.columns or "Times_Selected_Across_Folds" not in df.columns:
    raise ValueError("Input file must contain 'Feature' and 'Times_Selected_Across_Folds' columns.")

# -----------------------------
# Filter stable features
# -----------------------------
stable_df = df[df["Times_Selected_Across_Folds"] >= freq_threshold].copy()
stable_df = stable_df.sort_values(by="Times_Selected_Across_Folds", ascending=False)

# -----------------------------
# Save output
# -----------------------------
stable_df.to_csv(output_path, index=False)

# -----------------------------
# Summary
# -----------------------------
print(f"Loaded {len(df)} total features.")
print(f"Selected {len(stable_df)} stable features with frequency ≥ {freq_threshold}.")
print(f"Stable feature list saved to:\n{output_path}")
