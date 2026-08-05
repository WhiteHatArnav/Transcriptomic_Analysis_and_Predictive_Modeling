
import pandas as pd
from gprofiler import GProfiler
import os

# === 1) Input files and output folders/names for PRAD ===
file_mapping = {
    # White Race
    "/Users/arnavjoshi/Desktop/PRAD_GOandPE/WhiteRace_Downregulated_Genes.csv": {
        "folder": "/Users/arnavjoshi/Desktop/PRAD_GOandPE/WhiteDownReg",
        "prefix": "WhiteRace_DownRegulated"
    },
    "/Users/arnavjoshi/Desktop/PRAD_GOandPE/WhiteRace_Upregulated_Genes.csv": {
        "folder": "/Users/arnavjoshi/Desktop/PRAD_GOandPE/WhiteUpReg",
        "prefix": "WhiteRace_UpRegulated"
    },

    # Black Race
    "/Users/arnavjoshi/Desktop/PRAD_GOandPE/BlackRace_Downregulated_Genes.csv": {
        "folder": "/Users/arnavjoshi/Desktop/PRAD_GOandPE/BlackDownReg",
        "prefix": "BlackRace_DownRegulated"
    },
    "/Users/arnavjoshi/Desktop/PRAD_GOandPE/BlackRace_Upregulated_Genes.csv": {
        "folder": "/Users/arnavjoshi/Desktop/PRAD_GOandPE/BlackUpReg",
        "prefix": "BlackRace_UpRegulated"
    },

    # All Races
    "/Users/arnavjoshi/Desktop/PRAD_GOandPE/AllRace_Downregulated_Genes.csv": {
        "folder": "/Users/arnavjoshi/Desktop/PRAD_GOandPE/AllDownReg",
        "prefix": "AllRaces_DownRegulated"
    },
    "/Users/arnavjoshi/Desktop/PRAD_GOandPE/AllRace_Upregulated_Genes.csv": {
        "folder": "/Users/arnavjoshi/Desktop/PRAD_GOandPE/AllUpReg",
        "prefix": "AllRaces_UpRegulated"
    }
}

# === 2) Enrichment categories mapping for g:Profiler ===
gp_categories = {
    "GO:BP": "BPDirect",
    "GO:CC": "CCDirect",
    "GO:MF": "MFDirect",
    "KEGG": "KEGGPathway"
}

# Initialize g:Profiler
gp = GProfiler(return_dataframe=True)

# === 3) Function to run enrichment and save CSV ===
def run_gprofiler(file_path, out_folder, prefix):
    os.makedirs(out_folder, exist_ok=True)

    df = pd.read_csv(file_path)
    gene_list = df["gene_id"].dropna().astype(str).tolist()

    print(f"\nProcessing {prefix} ({len(gene_list)} genes)...")

    # Run g:Profiler
    res = gp.profile(
        organism="hsapiens",
        query=gene_list,
        sources=list(gp_categories.keys())
    )

    if res.empty:
        print("  No enrichment results found.")
        return

    # Map g:Profiler source to category label
    res["Category"] = res["source"].map(gp_categories)

    # Robust GO / pathway ID handling
    if "term_id" in res.columns:
        term_id_col = res["term_id"]
    elif "native" in res.columns:
        term_id_col = res["native"]
    else:
        raise RuntimeError("No GO or pathway identifier column found in g:Profiler output.")

    # Handle gene list output columns
    if "intersection_names" in res.columns:
        genes_col = res["intersection_names"].apply(
            lambda x: ",".join(x) if isinstance(x, list) else str(x)
        )
    elif "intersections" in res.columns:
        genes_col = res["intersections"].apply(str)
    elif "intersection" in res.columns:
        genes_col = res["intersection"].apply(
            lambda x: ",".join(x) if isinstance(x, list) else str(x)
        )
    else:
        genes_col = pd.Series([""] * len(res))

    # Compute Fold Enrichment
    df_out = pd.DataFrame({
        "Category": res["Category"],
        "Term": term_id_col.astype(str) + "~" + res["name"].astype(str),
        "Genes": genes_col,
        "Count": res["intersection_size"],
        "P-Value": res["p_value"],
        "Fold Enrichment": (
            res["intersection_size"] / len(gene_list)
        ) / (
            res["term_size"] / res["effective_domain_size"]
        )
    })

    # Save per-category results
    for cat_label in df_out["Category"].unique():
        out_file = os.path.join(
            out_folder,
            f"{prefix}_{cat_label}_Dotplot.csv"
        )
        df_out[df_out["Category"] == cat_label].to_csv(out_file, index=False)
        print(f"  Saved {out_file} ({len(df_out[df_out['Category'] == cat_label])} terms)")

# === Run enrichment for all PRAD race groups ===
for f, info in file_mapping.items():
    run_gprofiler(f, info["folder"], info["prefix"])



##############################################################################
# Code without GO Code number
##############################################################################


# import pandas as pd
# from gprofiler import GProfiler
# import os

# # === Input files and output folders/names for PRAD ===
# file_mapping = {
#     # --- White Race ---
#     "/Users/arnavjoshi/Desktop/PRAD_GOandPE/WhiteRace_Downregulated_Genes.csv": {
#         "folder": "/Users/arnavjoshi/Desktop/PRAD_GOandPE/WhiteDownReg",
#         "prefix": "WhiteRace_DownRegulated"
#     },
#     "/Users/arnavjoshi/Desktop/PRAD_GOandPE/WhiteRace_Upregulated_Genes.csv": {
#         "folder": "/Users/arnavjoshi/Desktop/PRAD_GOandPE/WhiteUpReg",
#         "prefix": "WhiteRace_UpRegulated"
#     },
    
#     # --- Black Race ---
#     "/Users/arnavjoshi/Desktop/PRAD_GOandPE/BlackRace_Downregulated_Genes.csv": {
#         "folder": "/Users/arnavjoshi/Desktop/PRAD_GOandPE/BlackDownReg",
#         "prefix": "BlackRace_DownRegulated"
#     },
#     "/Users/arnavjoshi/Desktop/PRAD_GOandPE/BlackRace_Upregulated_Genes.csv": {
#         "folder": "/Users/arnavjoshi/Desktop/PRAD_GOandPE/BlackUpReg",
#         "prefix": "BlackRace_UpRegulated"
#     },
    
#     # --- All Races ---
#     "/Users/arnavjoshi/Desktop/PRAD_GOandPE/AllRace_Downregulated_Genes.csv": {
#         "folder": "/Users/arnavjoshi/Desktop/PRAD_GOandPE/AllDownReg",
#         "prefix": "AllRaces_DownRegulated"
#     },
#     "/Users/arnavjoshi/Desktop/PRAD_GOandPE/AllRace_Upregulated_Genes.csv": {
#         "folder": "/Users/arnavjoshi/Desktop/PRAD_GOandPE/AllUpReg",
#         "prefix": "AllRaces_UpRegulated"
#     }
# }

# # ===  Enrichment categories mapping for g:Profiler ===
# gp_categories = {
#     "GO:BP": "BPDirect",
#     "GO:CC": "CCDirect",
#     "GO:MF": "MFDirect",
#     "KEGG": "KEGGPathway"
# }

# # Initialize g:Profiler
# gp = GProfiler(return_dataframe=True)

# # ===  Function to run enrichment and save CSV ===
# def run_gprofiler(file_path, out_folder, prefix):
#     os.makedirs(out_folder, exist_ok=True)
    
#     df = pd.read_csv(file_path)
#     gene_list = df['gene_id'].dropna().astype(str).tolist()
    
#     print(f"\nProcessing {prefix} ({len(gene_list)} genes)...")
    
#     # Run g:Profiler
#     res = gp.profile(
#         organism='hsapiens',
#         query=gene_list,
#         sources=list(gp_categories.keys())
#     )
    
#     if res.empty:
#         print("  No enrichment results found.")
#         return
    
#     # Map g:Profiler source to category label
#     res['Category'] = res['source'].map(gp_categories)
    
#     # Handle gene list output columns
#     if 'intersection_names' in res.columns:
#         genes_col = res['intersection_names'].apply(lambda x: ','.join(x) if isinstance(x, list) else str(x))
#     elif 'intersections' in res.columns:
#         genes_col = res['intersections'].apply(str)
#     elif 'intersection' in res.columns:
#         genes_col = res['intersection'].apply(lambda x: ','.join(x) if isinstance(x, list) else str(x))
#     else:
#         genes_col = pd.Series([''] * len(res))
    
#     # Compute Fold Enrichment
#     df_out = pd.DataFrame({
#         'Category': res['Category'],
#         'Term': res['name'],
#         'Genes': genes_col,
#         'Count': res['intersection_size'],
#         'P-Value': res['p_value'],
#         'Fold Enrichment': res['intersection_size'] / len(gene_list) / (res['term_size'] / res['effective_domain_size'])
#     })
    
#     # Save per-category results
#     for cat_label in df_out['Category'].unique():
#         out_file = os.path.join(out_folder, f"{prefix}_{cat_label}_Dotplot.csv")
#         df_out[df_out['Category'] == cat_label].to_csv(out_file, index=False)
#         print(f"  Saved {out_file} ({len(df_out[df_out['Category'] == cat_label])} terms)")

# # === Run enrichment for all PRAD race groups ===
# for f, info in file_mapping.items():
#     run_gprofiler(f, info['folder'], info['prefix'])
