import pandas as pd
from collections import Counter

# === File paths ===
input_expression_csv = "/Users/arnavjoshi/Desktop/LIHCTranscAnalysis/lihc_expression_matrix.csv"
clinical_tsv = "/Users/arnavjoshi/Desktop/LIHCOtherData/clinical.cart.2025-07-18/clinical.tsv"
output_expression_csv = "/Users/arnavjoshi/Desktop/LIHCTranscAnalysis/lihc_expression_matrix_with_race.csv"

# === Load expression matrix with multi-header (Patient ID, Sample ID) ===
expr_df = pd.read_csv(input_expression_csv, header=[0, 1], index_col=0)

# === Load clinical data ===
clinical_df = pd.read_csv(clinical_tsv, sep='\t', dtype=str)

# === Map Case ID (submitter_id) → Race ===
case_to_race = clinical_df.set_index("cases.submitter_id")["demographic.race"].to_dict()

# === Build Race row (using Case ID / Patient ID from level 0 of column MultiIndex) ===
case_ids = [col[0] for col in expr_df.columns]
race_values = [case_to_race.get(cid, "NA") for cid in case_ids]

# === Display race counts ===
race_counts = Counter(race_values)
print("\nRace distribution among samples:")
for race, count in race_counts.items():
    print(f"  {race}: {count}")

# === Insert Race row as a third header row ===
# First, create a new MultiIndex for 3 header rows
new_columns = pd.MultiIndex.from_arrays([
    [expr_df.columns[i][0] for i in range(len(expr_df.columns))],   # Patient ID
    [expr_df.columns[i][1] for i in range(len(expr_df.columns))],   # Sample ID
    race_values                                                     # Race
], names=["Patient ID", "Sample ID", "Race"])

# Assign the new 3-level header
expr_df.columns = new_columns

# === Save updated CSV with 3 header rows ===
expr_df.to_csv(output_expression_csv)

print(f"\nRace row added below Sample ID row and saved to:\n{output_expression_csv}")
