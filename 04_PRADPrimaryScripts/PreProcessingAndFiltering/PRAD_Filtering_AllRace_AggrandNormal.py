import pandas as pd

# === File paths (PRAD All races) ===
input_file = "/Users/arnavjoshi/Desktop/PRADTranscAnalysis/prad_expression_matrix_with_race_sampletype_and_subtype.csv"
output_expression_file = "/Users/arnavjoshi/Desktop/PRADDESeq/prad_all_expression_filtered.csv"
output_metadata_file = "/Users/arnavjoshi/Desktop/PRADDESeq/All_PRAD_Metadata.csv"

# === Load expression matrix with metadata rows ===
df = pd.read_csv(input_file, index_col=0)

# === Extract metadata rows ===
sample_id_row = df.loc["Sample ID"]
race_row = df.loc["Race"]
sample_type_row = df.loc["Sample Type"]
sample_subtype_row = df.loc["Sample Subtype"]  # Gleason scores stored here

# === Identify columns by condition ===
# Tumor = Gleason 9 or 10 (any race)
subtype_numeric = pd.to_numeric(sample_subtype_row, errors="coerce")

tumor_cols = df.columns[subtype_numeric.isin([9, 10])]
normal_cols = df.columns[sample_type_row == "Normal"]

# === Print counts ===
print(f"Number of Gleason 9/10 (Tumor) samples (all races): {len(tumor_cols)}")
print(f"Number of Normal samples (all races): {len(normal_cols)}")

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
