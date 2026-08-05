import pandas as pd

# === File paths (LIHC/lihc -> LUAD/luad) ===
input_file = "/Users/arnavjoshi/Desktop/LUADTranscAnalysis/luad_expression_matrix_with_race_sampletype_and_subtypes.csv"
output_expression_file = "/Users/arnavjoshi/Desktop/LUADDESeq/luad_white_expression_filtered.csv"
output_metadata_file = "/Users/arnavjoshi/Desktop/LUADDESeq/White_LUAD_Metadata.csv"

# === Load expression matrix with metadata rows ===
df = pd.read_csv(input_file, index_col=0)

# === Extract metadata rows ===
sample_id_row = df.loc["Sample ID"]
race_row = df.loc["Race"]
sample_type_row = df.loc["Sample Type"]
immune_subtype_row = df.loc["Immune Subtype"]

# === Identify columns by condition ===
immune_tumor_values = ["Wound Healing (Immune C1)", "IFN-gamma Dominant (Immune C2)"]

tumor_cols = df.columns[
    immune_subtype_row.isin(immune_tumor_values) &
    (race_row.str.lower() == "white")
]

normal_cols = df.columns[
    (sample_type_row == "Normal") &
    (race_row.str.lower() == "white")
]

# === Print counts ===
print(f"Number of Immune C1/C2 (Tumor) White samples: {len(tumor_cols)}")
print(f"Number of Normal White samples: {len(normal_cols)}")

# === Combine selected columns and filter DataFrame ===
selected_cols = tumor_cols.union(normal_cols)
df_filtered = df[selected_cols]

# === Drop metadata rows (keep gene expression only) ===
rows_to_drop = ["Race", "Immune Subtype", "Sample Type"]
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
