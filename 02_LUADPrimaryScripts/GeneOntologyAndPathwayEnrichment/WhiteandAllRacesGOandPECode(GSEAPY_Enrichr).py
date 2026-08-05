import pandas as pd
import gseapy as gp
from gseapy import enrichr
import os

# === 1️⃣ Define input files and corresponding output folders/names ===
file_mapping = {
    "/Users/arnavjoshi/Desktop/LUAD_GOandPE/WhiteRace_Downregulated_Genes.csv": {
        "folder": "/Users/arnavjoshi/Desktop/LUAD_GOandPE/WhiteDownReg",
        "prefix": "WhiteRace_DownRegulated"
    },
    "/Users/arnavjoshi/Desktop/LUAD_GOandPE/WhiteRace_Upregulated_Genes.csv": {
        "folder": "/Users/arnavjoshi/Desktop/LUAD_GOandPE/WhiteUpReg",
        "prefix": "WhiteRace_UpRegulated"
    },
    "/Users/arnavjoshi/Desktop/LUAD_GOandPE/AllRace_Downregulated_Genes.csv": {
        "folder": "/Users/arnavjoshi/Desktop/LUAD_GOandPE/AllDownReg",
        "prefix": "AllRaces_DownRegulated"
    },
    "/Users/arnavjoshi/Desktop/LUAD_GOandPE/AllRace_Upregulated_Genes.csv": {
        "folder": "/Users/arnavjoshi/Desktop/LUAD_GOandPE/AllUpReg",
        "prefix": "AllRaces_UpRegulated"
    }
}

# === 2️⃣ Enrichr categories mapping ===
enrichr_categories = {
    "GO_Biological_Process_2023": "BPDirect",
    "GO_Cellular_Component_2023": "CCDirect",
    "GO_Molecular_Function_2023": "MFDirect",
    "KEGG_2021_Human": "KEGGPathway"
}

# === 3️⃣ Run enrichment and save CSVs ===
def run_enrichment(file_path, out_folder, prefix):
    os.makedirs(out_folder, exist_ok=True)
    
    df = pd.read_csv(file_path)
    gene_list = df['gene_id'].dropna().astype(str).tolist()
    
    print(f"Processing {prefix} ({len(gene_list)} genes)...")
    
    for category, cat_label in enrichr_categories.items():
        enr = gp.enrichr(
        gene_list=gene_list,
        gene_sets=category,
        outdir=None,   # no automatic output
        cutoff=0.05
        )
        
        if enr.res2d is not None and not enr.res2d.empty:
            res = enr.res2d.copy()
            
            df_out = pd.DataFrame()
            df_out['Category'] = cat_label
            df_out['Term'] = res['Term']
            df_out['Genes'] = res['Genes']
            df_out['Count'] = res['Overlap'].str.split('/').str[0].astype(int)
            df_out['List Total'] = len(gene_list)
            df_out['Pop Hits'] = df_out['Count']
            df_out['P-Value'] = res['P-value']
            df_out['Benjamini'] = res['Adjusted P-value']
            df_out['Fold Enrichment'] = res['Combined Score']
            df_out['Bonferroni'] = res['Old P-value']
            df_out['FDR'] = res['Adjusted P-value']
            df_out['Fisher Exact'] = res['P-value']
            
            # Output file name as requested
            out_file = os.path.join(out_folder, f"{prefix}_{cat_label}_Dotplot.csv")
            df_out.to_csv(out_file, index=False)
            print(f"  Saved {out_file}")
        else:
            print(f"  No enrichment results for {cat_label}")

# === 4️⃣ Loop through all input files ===
for f, info in file_mapping.items():
    run_enrichment(f, info['folder'], info['prefix'])
