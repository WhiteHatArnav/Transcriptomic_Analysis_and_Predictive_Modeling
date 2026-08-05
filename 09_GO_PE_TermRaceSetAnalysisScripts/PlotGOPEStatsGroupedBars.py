import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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

# Human-readable labels for legend (optional but nicer)
CATEGORY_LABELS = {
    "N_Common_AllThree": "Common (All 3)",
    "N_Shared_White_Black": "Shared (White–Black)",
    "N_Shared_White_All": "Shared (White–All)",
    "N_Shared_Black_All": "Shared (Black–All)",
    "N_Unique_White": "Unique (White)",
    "N_Unique_Black": "Unique (Black)",
    "N_Unique_All": "Unique (All)",
}

CHARTTYPE_ORDER = ["BP", "CC", "MF", "KEGG"]


# =========================
# PLOTTING HELPERS
# =========================

def plot_matplotlib_grouped_bar(df, cancer, reg_suffix, out_dir):
    """
    Grouped bar chart using matplotlib.
    One group per ChartType, 7 bars per group (the categories).
    """
    # Ensure ChartType ordering
    df = df.copy()
    df["ChartType"] = pd.Categorical(df["ChartType"],
                                     categories=CHARTTYPE_ORDER,
                                     ordered=True)
    df = df.sort_values("ChartType")

    x = range(len(df))  # positions for each ChartType
    width = 0.1  # bar width (small to fit 7 bars)

    fig, ax = plt.subplots(figsize=(10, 6))

    # color cycle (7 categories)
    colors = plt.cm.tab10.colors  # at least 10 distinct colors

    for i, col in enumerate(CATEGORY_COLUMNS):
        # Offsetting each category around the group center
        offsets = [pos + (i - 3) * width for pos in x]  # center around each x
        ax.bar(
            offsets,
            df[col],
            width,
            label=CATEGORY_LABELS.get(col, col),
            color=colors[i % len(colors)],
        )

    ax.set_xticks(list(x))
    ax.set_xticklabels(df["ChartType"])
    ax.set_ylabel("Number of Terms")
    ax.set_title(f"{cancer} {reg_suffix} – GO/Pathway Term Overlap (matplotlib)")

    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.)

    plt.tight_layout()

    out_path = os.path.join(
        out_dir, f"{cancer}_{reg_suffix}_GO_PE_Stats_groupedbar_matplotlib.png"
    )
    plt.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved matplotlib plot: {out_path}")


def plot_seaborn_grouped_bar(df, cancer, reg_suffix, out_dir):
    """
    Grouped bar chart using seaborn (catplot / barplot on melted data).
    """
    df = df.copy()
    df["ChartType"] = pd.Categorical(df["ChartType"],
                                     categories=CHARTTYPE_ORDER,
                                     ordered=True)

    # Melt into long form: one row per ChartType-category combination
    df_melt = df.melt(
        id_vars=["ChartType"],
        value_vars=CATEGORY_COLUMNS,
        var_name="Category",
        value_name="Count",
    )

    # Replace category codes with nicer labels for legend / color
    df_melt["CategoryLabel"] = df_melt["Category"].map(
        lambda c: CATEGORY_LABELS.get(c, c)
    )

    sns.set(style="whitegrid", context="talk")

    fig, ax = plt.subplots(figsize=(10, 6))

    sns.barplot(
        data=df_melt,
        x="ChartType",
        y="Count",
        hue="CategoryLabel",
        ax=ax,
    )

    ax.set_ylabel("Number of Terms")
    ax.set_title(f"{cancer} {reg_suffix} – GO/Pathway Term Overlap (seaborn)")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.)

    plt.tight_layout()

    out_path = os.path.join(
        out_dir, f"{cancer}_{reg_suffix}_GO_PE_Stats_groupedbar_seaborn.png"
    )
    plt.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved seaborn plot: {out_path}")


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

            # Sanity check: filter only expected ChartTypes
            df_stats = df_stats[df_stats["ChartType"].isin(CHARTTYPE_ORDER)]

            if df_stats.empty:
                print(f"No valid ChartType rows in {stats_path}, skipping.")
                continue

            # Matplotlib version
            plot_matplotlib_grouped_bar(df_stats, cancer, reg_suffix, stats_dir)

            # Seaborn version
            plot_seaborn_grouped_bar(df_stats, cancer, reg_suffix, stats_dir)


if __name__ == "__main__":
    main()
