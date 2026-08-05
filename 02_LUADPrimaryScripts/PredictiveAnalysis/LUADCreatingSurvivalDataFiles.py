import pandas as pd
import numpy as np
import os

# File paths (LUAD)
expression_files = {
    "luad_white_survival_data.csv": "/Users/arnavjoshi/Desktop/LUADDESeq/luad_white_expression_filtered.csv",
    "luad_black_survival_data.csv": "/Users/arnavjoshi/Desktop/LUADDESeq/luad_black_expression_filtered.csv",
    "luad_all_survival_data.csv":   "/Users/arnavjoshi/Desktop/LUADDESeq/luad_all_expression_filtered.csv",
}

survival_file = "/Users/arnavjoshi/Desktop/LUADOtherData/LUADSurvivalData.tsv"
output_folder = "/Users/arnavjoshi/Desktop/LUADPredictiveAnalytics"

# Load survival data
surv_df = pd.read_csv(survival_file, sep="\t")
surv_df["sample"] = surv_df["sample"].astype(str)

# Create output directory if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Process each expression file
for output_name, file_path in expression_files.items():
    print(f"\nProcessing file: {file_path}")

    # Load expression data
    expr_df = pd.read_csv(file_path, index_col=0)

    # Transpose the DataFrame (samples become rows)
    expr_df = expr_df.T

    # Log2(count + 1) transformation
    expr_df = np.log2(expr_df + 1)

    # Modify sample IDs by removing the last character (to align with survival 'sample' IDs)
    expr_df.index = expr_df.index.str[:-1]
    expr_df.index.name = "sample"

    # Merge with survival data
    merged_df = expr_df.merge(
        surv_df[["sample", "OS", "OS.time"]],
        on="sample",
        how="left"
    )

    # Print first 5 entries of OS and OS.time
    print("Preview of OS and OS.time columns:")
    print(merged_df[["OS", "OS.time"]].head())

    # Save to output folder
    output_path = os.path.join(output_folder, output_name)
    merged_df.to_csv(output_path, index=False)

    print(f"Processed and saved: {output_path}")
