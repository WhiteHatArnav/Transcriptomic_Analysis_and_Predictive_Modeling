import pandas as pd
from collections import Counter

# === File paths (LUAD) ===
input_expression_csv = "/Users/arnavjoshi/Desktop/LUADTranscAnalysis/luad_expression_matrix_with_race_and_sampletype.csv"
subtype_tsv = "/Users/arnavjoshi/Desktop/LUADOtherData/LUADSubtypeData.tsv"
output_csv = "/Users/arnavjoshi/Desktop/LUADTranscAnalysis/luad_expression_matrix_with_race_sampletype_and_subtypes.csv"

# === Load expression matrix with multi-index header (4 levels: Patient ID, Sample ID, Race, Sample Type) ===
expr_df = pd.read_csv(input_expression_csv, header=[0, 1, 2, 3], index_col=0)

# === Load subtype data ===
subtype_df = pd.read_csv(subtype_tsv, sep='\t', dtype=str)

# Robustly pick the sample column name if needed
sample_col = "sample"
if sample_col not in subtype_df.columns:
    for cand in ["Sample", "sampleID", "SampleID"]:
        if cand in subtype_df.columns:
            sample_col = cand
            break

# Build lookup dicts
sample_to_integrative = subtype_df.set_index(sample_col)["Subtype_Integrative"].to_dict()
sample_to_immune = subtype_df.set_index(sample_col)["Subtype_Immune_Model_Based"].to_dict()

# === Trim last character of Sample IDs from the expression matrix (level 1) ===

trimmed_sample_ids = [col[1][:-1] if col[1] else "" for col in expr_df.columns]

# === Map subtype values aligned with trimmed sample IDs ===
integrative_values = [sample_to_integrative.get(sid, "NA") for sid in trimmed_sample_ids]
immune_values = [sample_to_immune.get(sid, "NA") for sid in trimmed_sample_ids]

# === Count and display subtype distributions ===
print("\nIntegrated Subtype distribution among samples:")
for subtype, count in Counter(integrative_values).items():
    print(f"  {subtype}: {count}")

print("\nImmune Subtype distribution among samples:")
for subtype, count in Counter(immune_values).items():
    print(f"  {subtype}: {count}")

# === Create new MultiIndex with 6 header rows:
# Patient ID, Sample ID, Race, Sample Type, Integrative Subtype, Immune Subtype
new_columns = pd.MultiIndex.from_arrays([
    [col[0] for col in expr_df.columns],  # Patient ID
    [col[1] for col in expr_df.columns],  # Sample ID (original)
    [col[2] for col in expr_df.columns],  # Race
    [col[3] for col in expr_df.columns],  # Sample Type
    integrative_values,                   # Integrative Subtype (new)
    immune_values                         # Immune Subtype (new)
], names=["Patient ID", "Sample ID", "Race", "Sample Type", "Integrative Subtype", "Immune Subtype"])

# === Apply updated MultiIndex to DataFrame ===
expr_df.columns = new_columns

# === Save updated CSV ===
expr_df.to_csv(output_csv)
print(f"\nIntegrative and Immune subtype rows added below Sample Type row and saved to:\n{output_csv}")
