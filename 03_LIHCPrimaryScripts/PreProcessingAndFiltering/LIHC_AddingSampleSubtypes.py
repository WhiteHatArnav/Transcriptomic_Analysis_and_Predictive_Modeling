import pandas as pd
from collections import Counter

# === File paths ===
input_expression_csv = "/Users/arnavjoshi/Desktop/LIHCTranscAnalysis/lihc_expression_matrix_with_race_and_sampletype.csv"
subtype_tsv = "/Users/arnavjoshi/Desktop/LIHCOtherData/PANCAN_iClusterData.tsv"
output_csv = "/Users/arnavjoshi/Desktop/LIHCTranscAnalysis/lihc_expression_matrix_with_race_sampletype_and_subtype.csv"

# === Load expression matrix with multi-index header (4 levels: Patient ID, Sample ID, Race, Sample Type) ===
expr_df = pd.read_csv(input_expression_csv, header=[0, 1, 2, 3], index_col=0)

# === Load subtype data ===
subtype_df = pd.read_csv(subtype_tsv, sep='\t', dtype=str)
sample_to_subtype = subtype_df.set_index("sample")["Subtype_Integrative"].to_dict()

# === Trim last character of Sample IDs from the expression matrix (level 1) ===
trimmed_sample_ids = [col[1][:-1] if col[1] else "" for col in expr_df.columns]

# === Map subtype values aligned with trimmed sample IDs ===
subtype_values = [sample_to_subtype.get(sid, "NA") for sid in trimmed_sample_ids]

# === Count and display subtype distribution ===
subtype_counts = Counter(subtype_values)
print("\nSample Subtype distribution among samples:")
for subtype, count in subtype_counts.items():
    print(f"  {subtype}: {count}")

# === Create new MultiIndex with 5 header rows:
# Patient ID, Sample ID, Race, Sample Type, Sample Subtype
new_columns = pd.MultiIndex.from_arrays([
    [col[0] for col in expr_df.columns],  # Patient ID
    [col[1] for col in expr_df.columns],  # Original Sample ID
    [col[2] for col in expr_df.columns],  # Race
    [col[3] for col in expr_df.columns],  # Sample Type
    subtype_values                        # Sample Subtype (new)
], names=["Patient ID", "Sample ID", "Race", "Sample Type", "Sample Subtype"])

# === Apply updated MultiIndex to DataFrame ===
expr_df.columns = new_columns

# === Save updated CSV ===
expr_df.to_csv(output_csv)
print(f"\nSample subtype row added below Sample Type row and saved to:\n{output_csv}")
