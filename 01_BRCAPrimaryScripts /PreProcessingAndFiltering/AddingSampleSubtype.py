import pandas as pd
from collections import Counter

# === File paths ===
input_expression_csv = "/Users/arnavjoshi/Desktop/BRCATranscAnalysis/brca_expression_matrix_with_race.csv"
subtype_tsv = "/Users/arnavjoshi/Desktop/BRCAOtherData/xenaDownload_sample_subtype.tsv"
output_csv = "/Users/arnavjoshi/Desktop/BRCATranscAnalysis/brca_expression_matrix_with_race_and_subtype.csv"

# === Load expression matrix with multi-index header ===
expr_df = pd.read_csv(input_expression_csv, header=[0, 1, 2], index_col=0)

# === Load subtype data ===
subtype_df = pd.read_csv(subtype_tsv, sep='\t', dtype=str)
sample_to_subtype = subtype_df.set_index("sample")["PAM50_mRNA_nature2012"].to_dict()

# === Trim last character of Sample IDs from the expression matrix (level 1) ===
trimmed_sample_ids = [col[1][:-1] if col[1] else "" for col in expr_df.columns]
subtype_values = [sample_to_subtype.get(sid, "NA") for sid in trimmed_sample_ids]

# === Display subtype counts ===
subtype_counts = Counter(subtype_values)
print("\nSample Subtype distribution among samples:")
for subtype, count in subtype_counts.items():
    print(f"  {subtype}: {count}")

# === Create new MultiIndex with 4 header rows (Patient ID, Sample ID, Race, Subtype) ===
new_columns = pd.MultiIndex.from_arrays([
    [col[0] for col in expr_df.columns],  # Patient ID
    [col[1] for col in expr_df.columns],  # Original Sample ID
    [col[2] for col in expr_df.columns],  # Race
    subtype_values                        # Sample Subtype
], names=["Patient ID", "Sample ID", "Race", "Sample Subtype"])

# === Apply updated MultiIndex to DataFrame ===
expr_df.columns = new_columns

# === Save updated CSV ===
expr_df.to_csv(output_csv)
print(f"\n Sample subtype row added below Race row and saved to:\n{output_csv}")
