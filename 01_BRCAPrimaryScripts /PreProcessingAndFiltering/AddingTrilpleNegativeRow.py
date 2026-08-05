import pandas as pd

# File paths
expression_file = "/Users/arnavjoshi/Desktop/BRCATranscAnalysis/brca_expression_matrix_with_race_subtype_and_sampletype.csv"
triple_negative_file = "/Users/arnavjoshi/Desktop/BRCAOtherData/TripleNegativeSampleList.xlsx"
output_file = "/Users/arnavjoshi/Desktop/BRCATranscAnalysis/brca_expression_matrix_with_triple_negative.csv"

# Load expression matrix with gene IDs as row index
df = pd.read_csv(expression_file, index_col=0)

# Extract "Sample ID" row (assumed to be the second row in the file)
sample_ids = df.iloc[0]  # index 0 corresponds to second row (after header)

# Remove last character from each Sample ID for matching
processed_sample_ids = sample_ids.apply(lambda x: str(x)[:-1])

# Load list of triple negative sample IDs
tn_samples = pd.read_excel(triple_negative_file)
tn_sample_set = set(tn_samples['sample'].astype(str))

# Create Triple Negative row based on match
triple_negative_row = processed_sample_ids.apply(lambda x: "Yes" if x in tn_sample_set else "No")

# Insert "Triple Negative" row just below "Sample Type"
# First, find the index of the "Sample Type" row
row_labels = list(df.index)
if "Sample Type" not in row_labels:
    raise ValueError('"Sample Type" row not found in the file.')

sample_type_index = row_labels.index("Sample Type")

# Create a new DataFrame with the row inserted in the correct position
df_upper = df.iloc[:sample_type_index + 1]  # includes "Sample Type"
df_lower = df.iloc[sample_type_index + 1:]  # everything below it

# Convert the Triple Negative row to a DataFrame
triple_negative_df = pd.DataFrame([triple_negative_row], index=["Triple Negative"])

# Concatenate all parts
df_updated = pd.concat([df_upper, triple_negative_df, df_lower])

# Save the updated file
df_updated.to_csv(output_file)
print(f"Updated file saved to: {output_file}")
