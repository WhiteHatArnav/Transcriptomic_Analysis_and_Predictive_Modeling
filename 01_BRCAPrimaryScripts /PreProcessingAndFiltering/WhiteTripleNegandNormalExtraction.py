import pandas as pd

# File paths (corrected input)
input_file = "/Users/arnavjoshi/Desktop/BRCATranscAnalysis/brca_expression_matrix_with_triple_negative.csv"
output_expression_file = "/Users/arnavjoshi/Desktop/BRCADESeq/brca_white_expression_filtered.csv"
output_metadata_file = "/Users/arnavjoshi/Desktop/BRCADESeq/White_BRCA_Metadata.csv"

# Load expression matrix with gene IDs as index
df = pd.read_csv(input_file, index_col=0)

# Extract metadata rows
sample_id_row = df.loc["Sample ID"]
race_row = df.loc["Race"]
sample_type_row = df.loc["Sample Type"]
triple_negative_row = df.loc["Triple Negative"]

# Find columns meeting each condition
tumor_cols = df.columns[(triple_negative_row == "Yes") & (race_row.str.lower() == "white")]
normal_cols = df.columns[(sample_type_row == "Solid Tissue Normal") & (race_row.str.lower() == "white")]

# Combine selected columns
selected_cols = tumor_cols.union(normal_cols)

# Filter the DataFrame to keep only selected columns
df_filtered = df[selected_cols]

# Drop metadata rows except "Sample ID"
rows_to_drop = ["Race", "Sample Subtype", "Sample Type", "Triple Negative"]
df_filtered = df_filtered.drop(index=rows_to_drop)

# Save the filtered expression matrix
df_filtered.to_csv(output_expression_file)

# Build the metadata DataFrame
sample_ids = sample_id_row[selected_cols]
condition_labels = ["Tumor" if col in tumor_cols else "Normal" for col in selected_cols]

metadata_df = pd.DataFrame({
    "Sample": sample_ids.values,
    "Condition": condition_labels
})

# Save the metadata file
metadata_df.to_csv(output_metadata_file, index=False)

print(f"Filtered expression matrix saved to: {output_expression_file}")
print(f"Metadata file saved to: {output_metadata_file}")
