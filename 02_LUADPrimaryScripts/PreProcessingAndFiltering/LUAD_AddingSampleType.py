import pandas as pd
from collections import Counter

# === File paths ===
input_expression_csv = "/Users/arnavjoshi/Desktop/LUADTranscAnalysis/luad_expression_matrix_with_race.csv"
sample_sheet_tsv = "/Users/arnavjoshi/Desktop/LUADOtherData/gdc_sample_sheet.2025-09-23.tsv"
output_csv = "/Users/arnavjoshi/Desktop/LUADTranscAnalysis/luad_expression_matrix_with_race_and_sampletype.csv"

# === Load expression matrix with 3-level MultiIndex (Patient ID, Sample ID, Race) ===
expr_df = pd.read_csv(input_expression_csv, header=[0, 1, 2], index_col=0)

# === Load sample type data ===
sample_df = pd.read_csv(sample_sheet_tsv, sep='\t', dtype=str)
sample_to_type = sample_df.set_index("Sample ID")["Tissue Type"].to_dict()

# === Extract sample IDs from level 1 of MultiIndex ===
sample_ids = [col[1] for col in expr_df.columns]

# === Build Sample Type row aligned with expression matrix columns ===
sample_type_values = [sample_to_type.get(sid, "NA") for sid in sample_ids]

# === Display sample type counts ===
sample_type_counts = Counter(sample_type_values)
print("\nSample Type distribution among samples:")
for stype, count in sample_type_counts.items():
    print(f"  {stype}: {count}")

# === Create new MultiIndex with 4 header rows:
# Patient ID, Sample ID, Race, Sample Type
new_columns = pd.MultiIndex.from_arrays([
    [col[0] for col in expr_df.columns],  # Patient ID
    [col[1] for col in expr_df.columns],  # Sample ID
    [col[2] for col in expr_df.columns],  # Race
    sample_type_values                    # Sample Type (new)
], names=["Patient ID", "Sample ID", "Race", "Sample Type"])

# === Apply updated MultiIndex to DataFrame ===
expr_df.columns = new_columns

# === Save updated CSV ===
expr_df.to_csv(output_csv)
print(f"\nSample Type row added below Race row and saved to:\n{output_csv}")
