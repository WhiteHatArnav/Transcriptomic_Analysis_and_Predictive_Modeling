import os
import pandas as pd
import matplotlib.pyplot as plt

BASE_SUMMARY_DIR = "/Users/arnavjoshi/Desktop/GOParentTermAnalysis/SummaryCSVs"
OUTPUT_BASE_DIR = "/Users/arnavjoshi/Desktop/GOParentTermAnalysis"

CANCERS = ["BRCA", "LUAD", "LIHC", "PRAD"]
REGULATIONS = ["UpReg", "DownReg"]
ONTOLOGIES = ["BP", "MF", "CC"]

CATEGORIES = [
    "Unique_Black",
    "Unique_White",
    "Unique_All",
    "Common_White_Black",
    "Common_Black_All",
    "Common_White_All",
    "Common_All_Three"
]

ONTOLOGY_COLORS = {
    "BP": "#1f77b4",
    "MF": "#ff7f0e",
    "CC": "#2ca02c"
}

def count_terms(csv_path):
    if not os.path.exists(csv_path):
        return {cat: 0 for cat in CATEGORIES}

    df = pd.read_csv(csv_path)

    counts = {}
    for cat in CATEGORIES:
        if cat not in df.columns:
            counts[cat] = 0
            continue

        counts[cat] = (
            df[cat]
            .dropna()
            .astype(str)
            .apply(lambda x: x.strip() != "")
            .sum()
        )

    return counts


for cancer in CANCERS:
    for reg in REGULATIONS:

        ontology_counts = {}

        for ont in ONTOLOGIES:
            csv_path = os.path.join(
                BASE_SUMMARY_DIR,
                f"{cancer}_{reg}_{ont}_GO_ParentTermSummary.csv"
            )
            ontology_counts[ont] = count_terms(csv_path)

        data = pd.DataFrame(ontology_counts).T
        data = data[CATEGORIES]

        if data.values.sum() == 0:
            continue

        fig, ax = plt.subplots(figsize=(12, 5))

        bottom = [0] * len(CATEGORIES)

        for ont in ONTOLOGIES:
            values = data.loc[ont].values
            ax.bar(
                CATEGORIES,
                values,
                bottom=bottom,
                label=ont,
                color=ONTOLOGY_COLORS[ont]
            )
            bottom = [b + v for b, v in zip(bottom, values)]

        ax.set_title(f"{cancer} {reg}: Parent GO Term Category Summary", fontsize=12)
        ax.set_ylabel("Number of Parent GO Terms")
        ax.set_xlabel("Race-Sharing Category")

        ax.set_xticklabels(CATEGORIES, rotation=30, ha="right")

        ax.legend(title="GO Ontology", frameon=False)

        plt.tight_layout()

        out_dir = os.path.join(OUTPUT_BASE_DIR, cancer)
        os.makedirs(out_dir, exist_ok=True)

        out_file = os.path.join(
            out_dir,
            f"{cancer}_{reg}_GO_CategoryStackedBar.png"
        )

        plt.savefig(out_file, dpi=300)
        plt.close()

        print(f"Saved: {out_file}")
