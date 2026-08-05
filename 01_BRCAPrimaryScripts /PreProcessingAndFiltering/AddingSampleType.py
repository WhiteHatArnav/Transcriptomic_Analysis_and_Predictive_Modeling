import pandas as pd
from collections import Counter

# === File paths ===
input_expression_csv = "/Users/arnavjoshi/Desktop/BRCATranscAnalysis/brca_expression_matrix_with_race_and_subtype.csv"
sample_type_tsv = "/Users/arnavjoshi/Desktop/BRCAOtherData/xenaDownload_sample_type.tsv"
output_csv = "/Users/arnavjoshi/Desktop/BRCATranscAnalysis/brca_expression_matrix_with_race_subtype_and_sampletype.csv"

# === Load expression matrix with multi-index header (4 levels) ===
expr_df = pd.read_csv(input_expression_csv, header=[0, 1, 2, 3], index_col=0)

# === Load sample type data ===
sample_type_df = pd.read_csv(sample_type_tsv, sep='\t', dtype=str)
sample_to_type = sample_type_df.set_index("sample")["sample_type"].to_dict()

# === Trim last character of Sample IDs from the expression matrix (level 1) ===
trimmed_sample_ids = [col[1][:-1] if col[1] else "" for col in expr_df.columns]

# === Build sample type row aligned with expression matrix columns ===
sample_type_values = [sample_to_type.get(sid, "NA") for sid in trimmed_sample_ids]

# === Display sample type counts ===
sample_type_counts = Counter(sample_type_values)
print("\nSample Type distribution among samples:")
for stype, count in sample_type_counts.items():
    print(f"  {stype}: {count}")

# === Create new MultiIndex with 5 header rows:
# Patient ID, Sample ID, Race, Sample Subtype, Sample Type
new_columns = pd.MultiIndex.from_arrays([
    [col[0] for col in expr_df.columns],  # Patient ID
    [col[1] for col in expr_df.columns],  # Sample ID (original, untrimmed)
    [col[2] for col in expr_df.columns],  # Race
    [col[3] for col in expr_df.columns],  # Sample Subtype
    sample_type_values                     # Sample Type (new)
], names=["Patient ID", "Sample ID", "Race", "Sample Subtype", "Sample Type"])

# === Apply updated MultiIndex to DataFrame ===
expr_df.columns = new_columns

# === Save updated CSV ===
expr_df.to_csv(output_csv)
print(f"\nSample type row added below Sample Subtype row and saved to:\n{output_csv}")
