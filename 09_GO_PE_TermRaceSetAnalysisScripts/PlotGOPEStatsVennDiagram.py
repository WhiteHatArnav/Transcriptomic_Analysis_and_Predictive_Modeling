import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3  

# =========================
# CONFIG
# =========================

CANCERS = ["BRCA", "LUAD", "PRAD", "LIHC"]
REG_SUFFIXES = ["UpReg", "DownReg"]  # must match your folder names
BASE_DIR = "/Users/arnavjoshi/Desktop/GOPETermsComparison"

# Columns in stats CSV corresponding to comparison categories
CATEGORY_COLUMNS = [
    "N_Common_AllThree",
    "N_Shared_White_Black",
    "N_Shared_White_All",
    "N_Shared_Black_All",
    "N_Unique_White",
    "N_Unique_Black",
    "N_Unique_All",
]

CHARTTYPE_ORDER = ["BP", "CC", "MF", "KEGG"]


# =========================
# PLOTTING: WHITE vs BLACK vs ALL VENN
# =========================

def plot_white_black_all_venn(df, cancer, reg_suffix, out_dir):
    """
    For each ChartType (BP, CC, MF, KEGG), create a 3-set Venn diagram
    of White vs Black vs All terms.

    A = White, B = Black, C = All

    Venn3 subset order (for a 7-tuple) is:
        (100, 010, 110, 001, 101, 011, 111)
      = (A only, B only, A∩B only, C only, A∩C only, B∩C only, A∩B∩C)

    Mapping from your columns:
        A only   (100) = N_Unique_White
        B only   (010) = N_Unique_Black
        A∩B only (110) = N_Shared_White_Black
        C only   (001) = N_Unique_All
        A∩C only (101) = N_Shared_White_All
        B∩C only (011) = N_Shared_Black_All
        A∩B∩C   (111) = N_Common_AllThree
    """
    df = df.copy()
    df["ChartType"] = pd.Categorical(
        df["ChartType"], categories=CHARTTYPE_ORDER, ordered=True
    )
    df = df.sort_values("ChartType")

    for _, row in df.iterrows():
        chart_type = row["ChartType"]

        # Pull counts; treat missing/NaN as 0
        n_common_allthree   = int(row.get("N_Common_AllThree", 0) or 0)
        n_shared_wb         = int(row.get("N_Shared_White_Black", 0) or 0)
        n_shared_wa         = int(row.get("N_Shared_White_All", 0) or 0)
        n_shared_ba         = int(row.get("N_Shared_Black_All", 0) or 0)
        n_unique_white      = int(row.get("N_Unique_White", 0) or 0)
        n_unique_black      = int(row.get("N_Unique_Black", 0) or 0)
        n_unique_all        = int(row.get("N_Unique_All", 0) or 0)

        subsets = (
            max(n_unique_white, 0),   # 100: White only
            max(n_unique_black, 0),   # 010: Black only
            max(n_shared_wb, 0),      # 110: White ∩ Black only
            max(n_unique_all, 0),     # 001: All only
            max(n_shared_wa, 0),      # 101: White ∩ All only
            max(n_shared_ba, 0),      # 011: Black ∩ All only
            max(n_common_allthree, 0) # 111: White ∩ Black ∩ All
        )

        # Skip if everything is zero
        if sum(subsets) == 0:
            print(f"All zeros for {cancer} {reg_suffix} {chart_type}, skipping Venn.")
            continue

        fig, ax = plt.subplots(figsize=(6, 6))

        venn3(
            subsets=subsets,
            set_labels=("White", "Black", "All"),
            ax=ax,
        )

        ax.set_title(
            f"{cancer} {reg_suffix} {chart_type} – White / Black / All term overlap"
        )

        plt.tight_layout()

        out_path = os.path.join(
            out_dir,
            f"{cancer}_{reg_suffix}_{chart_type}_WhiteBlackAll_Venn.png"
        )
        plt.savefig(out_path, dpi=300)
        plt.close(fig)

        print(f"Saved 3-set Venn diagram: {out_path}")


# =========================
# MAIN
# =========================

def main():
    for cancer in CANCERS:
        for reg_suffix in REG_SUFFIXES:
            stats_dir = os.path.join(BASE_DIR, cancer, reg_suffix)
            stats_file = f"{cancer}_{reg_suffix}_GO_PE_Stats.csv"
            stats_path = os.path.join(stats_dir, stats_file)

            if not os.path.isfile(stats_path):
                print(f"Stats file not found, skipping: {stats_path}")
                continue

            df_stats = pd.read_csv(stats_path)

            # Keep only expected ChartTypes
            df_stats = df_stats[df_stats["ChartType"].isin(CHARTTYPE_ORDER)]

            if df_stats.empty:
                print(f"No valid ChartType rows in {stats_path}, skipping.")
                continue

            # Generate White–Black–All Venn diagrams
            plot_white_black_all_venn(df_stats, cancer, reg_suffix, stats_dir)


if __name__ == "__main__":
    main()
