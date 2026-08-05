import os
import pandas as pd

# =========================
# CONFIG
# =========================

CANCERS = ["BRCA", "LUAD", "PRAD", "LIHC"]
RACES = ["White", "Black", "All"]
REGULATIONS = {
    "Up": "UpReg",
    "Down": "DownReg"
}
CHART_KEYS = ["BP", "CC", "MF", "KEGG"]  # BPDirect, CCDirect, MFDirect, KEGGPathway
TOP_N = 30  # top rows per file, matching your R code

INPUT_BASE_TEMPLATE = "/Users/arnavjoshi/Desktop/{cancer}_GOandPE"
OUTPUT_BASE_DIR = "/Users/arnavjoshi/Desktop/GOPETermsComparison"


# =========================
# HELPERS
# =========================

def read_term_column(file_path, top_n=TOP_N):
    """
    Read top_n terms from a GO/Pathway enrichment file (txt or csv).
    Returns a set of strings.
    """
    ext = os.path.splitext(file_path)[1].lower()
    try:
        if ext == ".txt":
            df = pd.read_csv(file_path, sep="\t")
        elif ext == ".csv":
            df = pd.read_csv(file_path)
        else:
            print(f"Warning: Skipping unsupported file type: {file_path}")
            return set()
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
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
    Find the first file in directory whose name contains chart_key
    and ends with .txt or .csv.
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

    return os.path.join(directory, candidates[0])


def collect_terms_for_race_chart(base_dir, race, reg_suffix, chart_key):
    """
    For a given cancer base dir, race, regulation suffix, and chart key,
    return the set of terms.
    """
    directory = os.path.join(base_dir, f"{race}{reg_suffix}")
    file_path = get_chart_file_for_key(directory, chart_key)
    if file_path is None:
        return set()
    return read_term_column(file_path)


def build_category_lists(white_terms, black_terms, all_terms):
    """
    Given term sets for White, Black, All, return dict of lists for:
    In_White, In_Black, In_All,
    Common_AllThree, Shared_White_Black, Shared_White_All, Shared_Black_All,
    Unique_White, Unique_Black, Unique_All
    Each value is a list of terms (no booleans).
    """

    # Base sets
    in_white = sorted(white_terms)
    in_black = sorted(black_terms)
    in_all = sorted(all_terms)

    # Category sets
    common_all_three = white_terms & black_terms & all_terms

    shared_white_black = (white_terms & black_terms) - all_terms
    shared_white_all = (white_terms & all_terms) - black_terms
    shared_black_all = (black_terms & all_terms) - white_terms

    unique_white = white_terms - (black_terms | all_terms)
    unique_black = black_terms - (white_terms | all_terms)
    unique_all = all_terms - (white_terms | black_terms)

    category_lists = {
        "In_White": sorted(in_white),
        "In_Black": sorted(in_black),
        "In_All": sorted(in_all),
        "Common_AllThree": sorted(common_all_three),
        "Shared_White_Black": sorted(shared_white_black),
        "Shared_White_All": sorted(shared_white_all),
        "Shared_Black_All": sorted(shared_black_all),
        "Unique_White": sorted(unique_white),
        "Unique_Black": sorted(unique_black),
        "Unique_All": sorted(unique_all),
    }

    return category_lists


def lists_to_dataframe(category_lists, column_order):
    """
    Turn a dict {col: list_of_terms} into a DataFrame where each column
    is that list of terms, padded with empty strings so all columns have
    the same length.
    """
    max_len = max(len(category_lists[col]) for col in column_order)
    data = {}
    for col in column_order:
        lst = category_lists[col]
        padded = lst + [""] * (max_len - len(lst))
        data[col] = padded
    return pd.DataFrame(data)


# =========================
# MAIN
# =========================

def main():
    # Column order exactly as requested
    column_order = [
        "In_White",
        "In_Black",
        "In_All",
        "Common_AllThree",
        "Shared_White_Black",
        "Shared_White_All",
        "Shared_Black_All",
        "Unique_White",
        "Unique_Black",
        "Unique_All",
    ]

    for cancer in CANCERS:
        input_base_dir = INPUT_BASE_TEMPLATE.format(cancer=cancer)

        for reg_label, reg_suffix in REGULATIONS.items():
            # Output dir: /GOPETermsComparison/BRCA/UpReg etc.
            reg_output_dir = os.path.join(OUTPUT_BASE_DIR, cancer, reg_suffix)
            os.makedirs(reg_output_dir, exist_ok=True)

            stats_rows = []

            for chart_key in CHART_KEYS:
                # Collect terms for each race
                white_terms = collect_terms_for_race_chart(
                    input_base_dir, "White", reg_suffix, chart_key
                )
                black_terms = collect_terms_for_race_chart(
                    input_base_dir, "Black", reg_suffix, chart_key
                )
                all_terms = collect_terms_for_race_chart(
                    input_base_dir, "All", reg_suffix, chart_key
                )

                # Build category lists
                category_lists = build_category_lists(white_terms, black_terms, all_terms)

                # Convert to DataFrame where each column is a list of terms
                df_chart = lists_to_dataframe(category_lists, column_order)

                # Save per-chart CSV (same naming as before)
                out_chart_name = f"{cancer}_{reg_suffix}_{chart_key}_TermComparison.csv"
                out_chart_path = os.path.join(reg_output_dir, out_chart_name)
                df_chart.to_csv(out_chart_path, index=False)
                print(f"Saved chart comparison: {out_chart_path}")

                # Stats row for this chart
                stats_rows.append({
                    "Cancer": cancer,
                    "Regulation": reg_label,       # Up / Down (human label)
                    "RegulationFolder": reg_suffix,  # UpReg / DownReg (folder label)
                    "ChartType": chart_key,
                    "N_In_White": len(category_lists["In_White"]),
                    "N_In_Black": len(category_lists["In_Black"]),
                    "N_In_All": len(category_lists["In_All"]),
                    "N_Common_AllThree": len(category_lists["Common_AllThree"]),
                    "N_Shared_White_Black": len(category_lists["Shared_White_Black"]),
                    "N_Shared_White_All": len(category_lists["Shared_White_All"]),
                    "N_Shared_Black_All": len(category_lists["Shared_Black_All"]),
                    "N_Unique_White": len(category_lists["Unique_White"]),
                    "N_Unique_Black": len(category_lists["Unique_Black"]),
                    "N_Unique_All": len(category_lists["Unique_All"]),
                })

            # One stats CSV per cancer & regulation
            df_stats = pd.DataFrame(stats_rows)
            stats_name = f"{cancer}_{reg_suffix}_GO_PE_Stats.csv"
            stats_path = os.path.join(reg_output_dir, stats_name)
            df_stats.to_csv(stats_path, index=False)
            print(f"Saved stats: {stats_path}")


if __name__ == "__main__":
    main()
