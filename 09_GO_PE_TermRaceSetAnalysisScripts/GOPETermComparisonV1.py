import os
import pandas as pd

# =========================
# CONFIGURABLE PARAMETERS
# =========================

CANCERS = ["BRCA", "LUAD", "PRAD", "LIHC"]
RACES = ["White", "Black", "All"]
REGULATIONS = {
    "Up": "UpReg",
    "Down": "DownReg"
}
CHART_KEYS = ["BP", "CC", "MF", "KEGG"]  # BPDirect, CCDirect, MFDirect, KEGGPathway
TOP_N = 30  # number of top rows per file (matching dotplot code)
BASE_TEMPLATE = "/Users/arnavjoshi/Desktop/{cancer}_GOandPE"
OUTPUT_DIR = "/Users/arnavjoshi/Desktop/GOPETermsComparison"


# =========================
# HELPER FUNCTIONS
# =========================

def read_term_column(file_path, top_n=TOP_N):
    """
    Read the enrichment result file and return a set of Term values
    (top_n rows), stripped of leading/trailing whitespace.
    Handles .txt (tab-delimited) and .csv (comma-delimited).
    """
    ext = os.path.splitext(file_path)[1].lower()
    if ext == ".txt":
        df = pd.read_csv(file_path, sep="\t")
    elif ext == ".csv":
        df = pd.read_csv(file_path)
    else:
        print(f"Warning: Skipping unsupported file type: {file_path}")
        return set()

    if "Term" not in df.columns:
        print(f"Warning: 'Term' column not found in {file_path}, skipping.")
        return set()

    df_top = df.iloc[:top_n]
    terms = (
        df_top["Term"]
        .dropna()
        .astype(str)
        .str.strip()
    )
    return set(terms)


def get_chart_file_for_key(directory, chart_key):
    """
    In the given directory, find the first file that contains `chart_key`
    (e.g., 'BP', 'CC', 'MF', 'KEGG') in its filename.
    Returns full path or None if not found.
    """
    if not os.path.isdir(directory):
        print(f"Warning: Directory does not exist: {directory}")
        return None

    candidates = [
        f for f in os.listdir(directory)
        if chart_key in f and (f.endswith(".txt") or f.endswith(".csv"))
    ]

    if not candidates:
        print(f"Warning: No file found with key '{chart_key}' in {directory}")
        return None

    # If multiple, just take the first (you can refine if needed)
    chosen = candidates[0]
    return os.path.join(directory, chosen)


def collect_terms_for_race_dir(directory):
    """
    For a specific race/regulation directory, return a dict:
    { chart_key: set_of_terms }
    where chart_key in CHART_KEYS.
    """
    terms_by_chart = {}
    for key in CHART_KEYS:
        file_path = get_chart_file_for_key(directory, key)
        if file_path is None:
            terms_by_chart[key] = set()
        else:
            terms_by_chart[key] = read_term_column(file_path)
    return terms_by_chart


def classify_term(in_white, in_black, in_all):
    """
    Classify term into:
    - 'Common_AllThree'
    - 'Shared_White_Black'
    - 'Shared_White_All'
    - 'Shared_Black_All'
    - 'Unique_White'
    - 'Unique_Black'
    - 'Unique_All'
    Returns None if term is in none (shouldn't happen if union used).
    """
    total = int(in_white) + int(in_black) + int(in_all)

    if total == 3:
        return "Common_AllThree"
    elif total == 2:
        if in_white and in_black:
            return "Shared_White_Black"
        if in_white and in_all:
            return "Shared_White_All"
        if in_black and in_all:
            return "Shared_Black_All"
    elif total == 1:
        if in_white:
            return "Unique_White"
        if in_black:
            return "Unique_Black"
        if in_all:
            return "Unique_All"

    return None


# =========================
# MAIN LOGIC
# =========================

def main():
    # Ensure output directory exists
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for cancer in CANCERS:
        base_dir = BASE_TEMPLATE.format(cancer=cancer)

        for reg_label, reg_suffix in REGULATIONS.items():
            # Collect terms for each race and each chart
            race_terms = {}  # { race: { chart_key: set_of_terms } }

            for race in RACES:
                dir_path = os.path.join(base_dir, f"{race}{reg_suffix}")
                race_terms[race] = collect_terms_for_race_dir(dir_path)

            # Build comparison rows
            rows = []

            for chart_key in CHART_KEYS:
                white_set = race_terms.get("White", {}).get(chart_key, set())
                black_set = race_terms.get("Black", {}).get(chart_key, set())
                all_set = race_terms.get("All", {}).get(chart_key, set())

                union_terms = sorted(white_set | black_set | all_set)

                for term in union_terms:
                    in_white = term in white_set
                    in_black = term in black_set
                    in_all = term in all_set

                    category = classify_term(in_white, in_black, in_all)
                    if category is None:
                        continue  # Shouldn't happen, but for safety

                    row = {
                        "Cancer": cancer,
                        "Regulation": reg_label,        # Up / Down
                        "ChartType": chart_key,         # BP / CC / MF / KEGG
                        "Term": term,
                        "In_White": int(in_white),
                        "In_Black": int(in_black),
                        "In_All": int(in_all),
                        "GroupCategory": category       # common/shared/unique pattern
                    }
                    rows.append(row)

            df_out = pd.DataFrame(rows)

            # Output filename, e.g. BRCA_UpReg_GO_PE_TermComparison.csv
            out_name = f"{cancer}_{reg_label}Reg_GO_PE_TermComparison.csv"
            out_path = os.path.join(OUTPUT_DIR, out_name)

            df_out.to_csv(out_path, index=False)
            print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
