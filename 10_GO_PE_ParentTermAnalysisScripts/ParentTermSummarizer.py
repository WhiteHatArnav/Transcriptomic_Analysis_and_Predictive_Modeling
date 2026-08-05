import os
import pandas as pd

BASE_DIR = "/Users/arnavjoshi/Desktop/GOParentTermAnalysis"
OUT_DIR = "/Users/arnavjoshi/Desktop/GOParentTermAnalysis/SummaryCSVs"

CANCERS = ["BRCA", "LUAD", "LIHC", "PRAD"]
REGULATIONS = ["UpReg", "DownReg"]
ONTOLOGIES = ["BP", "MF", "CC"]

CATEGORY_FILES = {
    "Unique_Black": "Unique_Black_ParentTermSummary.csv",
    "Unique_White": "Unique_White_ParentTermSummary.csv",
    "Unique_All": "Unique_All_ParentTermSummary.csv",
    "Common_White_Black": "Shared_White_Black_ParentTermSummary.csv",
    "Common_Black_All": "Shared_Black_All_ParentTermSummary.csv",
    "Common_White_All": "Shared_White_All_ParentTermSummary.csv",
    "Common_All_Three": "Common_AllThree_ParentTermSummary.csv"
}

def extract_parent_terms(filepath):
    if not os.path.exists(filepath):
        return []

    df = pd.read_csv(filepath)

    if "parentTerm" not in df.columns:
        return []

    terms = (
        df["parentTerm"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    return sorted(terms)

os.makedirs(OUT_DIR, exist_ok=True)

for cancer in CANCERS:
    for reg in REGULATIONS:
        for ont in ONTOLOGIES:

            column_lists = {}

            max_len = 0
            for header, filename in CATEGORY_FILES.items():
                in_path = os.path.join(
                    BASE_DIR,
                    cancer,
                    reg,
                    ont,
                    filename
                )

                terms = extract_parent_terms(in_path)
                column_lists[header] = terms
                max_len = max(max_len, len(terms))

            padded_columns = {}
            for header, terms in column_lists.items():
                padded_columns[header] = terms + [""] * (max_len - len(terms))

            out_df = pd.DataFrame(padded_columns)

            out_name = f"{cancer}_{reg}_{ont}_GO_ParentTermSummary.csv"
            out_path = os.path.join(OUT_DIR, out_name)

            out_df.to_csv(out_path, index=False)

            print(f"Rewritten: {out_path}")
