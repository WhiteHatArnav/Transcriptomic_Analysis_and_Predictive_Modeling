import pandas as pd
import os

# Input file paths for LIHC
input_files = {
    "WhiteRace": "/Users/arnavjoshi/Desktop/LIHCDESeq/LIHCFinalWhiteRace_DESeqTestResults.csv",
    "BlackRace": "/Users/arnavjoshi/Desktop/LIHCDESeq/LIHC_BlackRace_DESeqTestResults.csv",
    "AllRace": "/Users/arnavjoshi/Desktop/LIHCDESeq/LIHC_AllRace_DESeqTestResults.csv"
}

# Output folder for LIHC
output_dir = "/Users/arnavjoshi/Desktop/LIHC_GOandPE"
os.makedirs(output_dir, exist_ok=True)

# Function to process each file
def process_deseq_file(label, filepath):
    # Load CSV
    df = pd.read_csv(filepath)

    # Strip version numbers from the gene ID column (assumed to be the first column)
    df['gene_id_clean'] = df.iloc[:, 0].astype(str).str.replace(r'\..*', '', regex=True)

    # Filter upregulated genes
    up = df[(df['padj'] < 0.05) & (df['log2FoldChange'] > 0)]
    up_genes = up['gene_id_clean'].dropna().unique()

    # Filter downregulated genes
    down = df[(df['padj'] < 0.05) & (df['log2FoldChange'] < 0)]
    down_genes = down['gene_id_clean'].dropna().unique()

    # Save to CSV
    up_output_path = os.path.join(output_dir, f"{label}_Upregulated_Genes.csv")
    down_output_path = os.path.join(output_dir, f"{label}_Downregulated_Genes.csv")

    pd.Series(up_genes, name='gene_id').to_csv(up_output_path, index=False)
    pd.Series(down_genes, name='gene_id').to_csv(down_output_path, index=False)

    print(f"Processed {label}: {len(up_genes)} up, {len(down_genes)} down")

# Run for each LIHC file
for label, filepath in input_files.items():
    process_deseq_file(label, filepath)
