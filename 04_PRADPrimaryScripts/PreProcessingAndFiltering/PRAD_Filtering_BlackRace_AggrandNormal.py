import pandas as pd

# === File paths (PRAD) ===
input_file = "/Users/arnavjoshi/Desktop/PRADTranscAnalysis/prad_expression_matrix_with_race_sampletype_and_subtype.csv"
output_expression_file = "/Users/arnavjoshi/Desktop/PRADDESeq/prad_black_expression_filtered.csv"
output_metadata_file = "/Users/arnavjoshi/Desktop/PRADDESeq/Black_PRAD_Metadata.csv"

# === Load expression matrix with metadata rows ===
df = pd.read_csv(input_file, index_col=0)

# === Extract metadata rows ===
sample_id_row = df.loc["Sample ID"]
race_row = df.loc["Race"]
sample_type_row = df.loc["Sample Type"]
sample_subtype_row = df.loc["Sample Subtype"]  # Gleason scores stored here

# === Normalize helpers ===
race_lower = race_row.astype(str).str.lower().fillna("")

# === Identify columns by condition ===
# Tumor = Gleason 9 or 10 AND race is Black or African American
# Coerce subtype to numeric to safely catch "9"/"10" or 9/10
subtype_numeric = pd.to_numeric(sample_subtype_row, errors="coerce")
is_black = race_lower == "black or african american"

tumor_cols = df.columns[subtype_numeric.isin([9, 10]) & is_black]
normal_cols = df.columns[(sample_type_row == "Normal") & is_black]

# === Print counts ===
print(f"Number of Gleason 9/10 (Tumor) Black samples: {len(tumor_cols)}")
print(f"Number of Normal Black samples: {len(normal_cols)}")

# === Combine selected columns and preserve original column order ===
selected_cols = [c for c in df.columns if (c in tumor_cols) or (c in normal_cols)]
df_filtered = df[selected_cols]

# === Drop metadata rows (keep gene expression only) ===
rows_to_drop = ["Race", "Sample Subtype", "Sample Type"]
df_filtered = df_filtered.drop(index=rows_to_drop)

# === Save filtered expression matrix ===
df_filtered.to_csv(output_expression_file)

# === Create and save metadata ===
sample_ids = sample_id_row[selected_cols]
condition_labels = ["Tumor" if col in tumor_cols else "Normal" for col in selected_cols]

metadata_df = pd.DataFrame({
    "Sample": sample_ids.values,
    "Condition": condition_labels
})

metadata_df.to_csv(output_metadata_file, index=False)

print(f"\nFiltered expression matrix saved to: {output_expression_file}")
print(f"Metadata file saved to: {output_metadata_file}")
