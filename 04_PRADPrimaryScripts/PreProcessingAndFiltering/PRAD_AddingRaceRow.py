import pandas as pd
from collections import Counter

# === File paths ===
input_expression_csv = "/Users/arnavjoshi/Desktop/PRADTranscAnalysis/prad_expression_matrix_with_race_and_sampletype.csv"
subtype_tsv = "/Users/arnavjoshi/Desktop/PRADOtherData/clinical.cart.2025-10-28/clinical.tsv"
output_csv = "/Users/arnavjoshi/Desktop/PRADTranscAnalysis/prad_expression_matrix_with_race_sampletype_and_subtype.csv"

# === Load expression matrix with multi-index header (4 levels: Patient ID, Sample ID, Race, Sample Type) ===
expr_df = pd.read_csv(input_expression_csv, header=[0, 1, 2, 3], index_col=0)

# === Load clinical (subtype) data: map cases.submitter_id -> diagnoses.gleason_score ===
clinical_df = pd.read_csv(subtype_tsv, sep='\t', dtype=str)

# Keep only the needed columns, drop full-NA rows, and fill missing scores with "NA"
needed_cols = ["cases.submitter_id", "diagnoses.gleason_score"]
clinical_df = clinical_df[needed_cols].dropna(how="all")
clinical_df["diagnoses.gleason_score"] = clinical_df["diagnoses.gleason_score"].fillna("NA")

# If there are duplicate submitter_ids, keep the first non-null score
clinical_df = clinical_df.drop_duplicates(subset=["cases.submitter_id"], keep="first")

# Build mapping dict (PATIENT-level)
case_to_gleason = clinical_df.set_index("cases.submitter_id")["diagnoses.gleason_score"].to_dict()

# === Align Gleason scores to the expression matrix columns using Patient ID (level 0) ===
patient_ids = [col[0] for col in expr_df.columns]  # level-0 header = Patient ID
subtype_values = [case_to_gleason.get(pid, "NA") for pid in patient_ids]

# === Count and display subtype (Gleason) distribution across columns ===
subtype_counts = Counter(subtype_values)
print("\nGleason score distribution among expression-matrix columns (by Patient ID):")
for subtype, count in sorted(subtype_counts.items(), key=lambda x: (str(x[0]), x[1])):
    print(f"  {subtype}: {count}")

# === Create new MultiIndex with 5 header rows:
# Patient ID, Sample ID, Race, Sample Type, Sample Subtype (Gleason score)
new_columns = pd.MultiIndex.from_arrays([
    [col[0] for col in expr_df.columns],  # Patient ID
    [col[1] for col in expr_df.columns],  # Sample ID
    [col[2] for col in expr_df.columns],  # Race
    [col[3] for col in expr_df.columns],  # Sample Type
    subtype_values                        # Sample Subtype (Gleason score)
], names=["Patient ID", "Sample ID", "Race", "Sample Type", "Sample Subtype"])

# === Apply updated MultiIndex to DataFrame ===
expr_df.columns = new_columns

# === Save updated CSV ===
expr_df.to_csv(output_csv)
print(f"\nSample subtype (Gleason score) row added below Sample Type row and saved to:\n{output_csv}")
