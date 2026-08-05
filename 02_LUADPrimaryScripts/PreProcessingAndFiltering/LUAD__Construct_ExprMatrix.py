import os
import pandas as pd

# Paths for LUAD
data_dir = "/Users/arnavjoshi/Desktop/LUAD"
sample_sheet_path = "/Users/arnavjoshi/Desktop/LUADOtherData/gdc_sample_sheet.2025-09-23.tsv"
output_csv_path = "/Users/arnavjoshi/Desktop/LUADTranscAnalysis/luad_expression_matrix.csv"

# Step 1: Load sample sheet and build mappings
sample_df = pd.read_csv(sample_sheet_path, sep='\t', dtype=str)
file_to_case = dict(zip(sample_df['File Name'], sample_df['Case ID']))
file_to_sample = dict(zip(sample_df['File Name'], sample_df['Sample ID']))

# Step 2: Initialize expression matrix
expression_matrix = {}
sample_id_row = {}
gene_ids = None

# Step 3: Walk through LUAD directory and process each TSV
for folder_name in os.listdir(data_dir):
    folder_path = os.path.join(data_dir, folder_name)
    if not os.path.isdir(folder_path):
        continue

    tsv_files = [f for f in os.listdir(folder_path) if f.endswith('.tsv')]
    if not tsv_files:
        continue

    tsv_file = tsv_files[0]
    if tsv_file not in file_to_case or tsv_file not in file_to_sample:
        print(f"File Name not found in sample sheet: {tsv_file}")
        continue

    case_id = file_to_case[tsv_file]
    sample_id = file_to_sample[tsv_file]
    tsv_path = os.path.join(folder_path, tsv_file)

    try:
        df = pd.read_csv(tsv_path, sep='\t', comment='#', usecols=['gene_id', 'unstranded'])
        df['gene_id'] = df['gene_id'].str.split('|').str[0]  # Strip version suffix

        if gene_ids is None:
            gene_ids = df['gene_id'].tolist()

        expression_matrix[case_id] = df['unstranded'].values
        sample_id_row[case_id] = sample_id

    except Exception as e:
        print(f"Error reading {tsv_path}: {e}")

# Step 4: Create DataFrame and insert second header for Sample IDs
if expression_matrix:
    matrix_df = pd.DataFrame(expression_matrix, index=gene_ids)
    matrix_df.index.name = "Gene ID"

    # Insert a second row for Sample IDs
    sample_id_series = pd.Series(sample_id_row)
    sample_id_df = pd.DataFrame([["Sample ID"] + [sample_id_series.get(col, "") for col in matrix_df.columns]],
                                columns=["Gene ID"] + list(matrix_df.columns))

    # Reset index and prepend Sample ID row
    matrix_df_reset = matrix_df.reset_index()
    combined_df = pd.concat([sample_id_df, matrix_df_reset], ignore_index=True)

    # Save final CSV
    combined_df.to_csv(output_csv_path, index=False)
    print(f"Expression matrix with Sample IDs saved to: {output_csv_path}")
else:
    print("No expression data was processed.")
