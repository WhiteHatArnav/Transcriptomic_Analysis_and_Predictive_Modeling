import pandas as pd
from collections import Counter

# === Input file path ===
input_csv = "/Users/arnavjoshi/Desktop/PRADTranscAnalysis/prad_expression_matrix_with_race_and_sampletype.csv"

# === Load expression matrix with 4-level MultiIndex header ===
expr_df = pd.read_csv(input_csv, header=[0, 1, 2, 3], index_col=0)

# === Extract Race and Sample Type from MultiIndex ===
race_list = [col[2].strip().lower() if isinstance(col[2], str) else "na" for col in expr_df.columns]
sample_type_list = [col[3].strip().lower() if isinstance(col[3], str) else "na" for col in expr_df.columns]

# === Count combinations of Race and Sample Type ===
combo_counts = Counter(zip(race_list, sample_type_list))

# === Define helper to print neatly ===
def get_count(race_key, type_key):
    return combo_counts.get((race_key.lower(), type_key.lower()), 0)

# === Display results ===
print("\nRace × Sample Type Distribution:")
print(f"  Black or African American & Tumor: {get_count('black or african american', 'tumor')}")
print(f"  Black or African American & Normal: {get_count('black or african american', 'normal')}")
print(f"  White & Tumor: {get_count('white', 'tumor')}")
print(f"  White & Normal: {get_count('white', 'normal')}")

# === Optional: print all unique combinations for QC ===
print("\nAll combinations found:")
for (race, stype), count in combo_counts.items():
    print(f"  {race} & {stype}: {count}")
