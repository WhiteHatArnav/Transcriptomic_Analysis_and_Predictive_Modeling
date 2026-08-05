import os
import pandas as pd
import matplotlib.pyplot as plt
import textwrap

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

MAX_TERMS_SHOWN = 5
WRAP_WIDTH = 22
FONT_SIZE = 7

HEIGHT_SCALE = 1.6
MIN_BAR_HEIGHT = 1.2


def extract_terms(csv_path):
    if not os.path.exists(csv_path):
        return {cat: [] for cat in CATEGORIES}

    df = pd.read_csv(csv_path)

    term_dict = {}
    for cat in CATEGORIES:
        if cat not in df.columns:
            term_dict[cat] = []
            continue

        terms = (
            df[cat]
            .dropna()
            .astype(str)
            .apply(str.strip)
            .loc[lambda x: x != ""]
            .unique()
            .tolist()
        )

        term_dict[cat] = sorted(terms)

    return term_dict


for cancer in CANCERS:
    for reg in REGULATIONS:

        ontology_terms = {}

        for ont in ONTOLOGIES:
            csv_path = os.path.join(
                BASE_SUMMARY_DIR,
                f"{cancer}_{reg}_{ont}_GO_ParentTermSummary.csv"
            )
            ontology_terms[ont] = extract_terms(csv_path)

        fig, ax = plt.subplots(figsize=(16, 7))

        bottom = [0] * len(CATEGORIES)

        for ont in ONTOLOGIES:
            heights = []
            labels = []

            for cat in CATEGORIES:
                terms = ontology_terms[ont][cat]
                raw_count = len(terms)

                if raw_count == 0:
                    heights.append(0)
                    labels.append("")
                else:
                    scaled_height = max(
                        raw_count * HEIGHT_SCALE,
                        MIN_BAR_HEIGHT
                    )
                    heights.append(scaled_height)

                    shown = terms[:MAX_TERMS_SHOWN]
                    text = "\n".join(
                        textwrap.fill(t, WRAP_WIDTH) for t in shown
                    )

                    if raw_count > MAX_TERMS_SHOWN:
                        text += f"\n(+{raw_count - MAX_TERMS_SHOWN} more)"

                    labels.append(text)

            bars = ax.bar(
                CATEGORIES,
                heights,
                bottom=bottom,
                color=ONTOLOGY_COLORS[ont],
                label=ont
            )

            for bar, label, btm, h in zip(bars, labels, bottom, heights):
                if h > 0:
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        btm + h / 2,
                        label,
                        ha="center",
                        va="center",
                        fontsize=FONT_SIZE
                    )

            bottom = [b + h for b, h in zip(bottom, heights)]

        ax.set_title(
            f"{cancer} {reg}: Parent GO Term Categories",
            fontsize=14,
            pad=12
        )

        ax.set_ylabel("Relative Parent GO Term Representation")
        ax.set_xlabel("Race-Sharing Category")

        ax.set_xticklabels(CATEGORIES, rotation=30, ha="right")

        ax.legend(title="GO Ontology", frameon=False)

        plt.tight_layout()

        out_dir = os.path.join(OUTPUT_BASE_DIR, cancer)
        os.makedirs(out_dir, exist_ok=True)

        out_file = os.path.join(
            out_dir,
            f"{cancer}_{reg}_GO_CategoryStackedBar_Annotated.png"
        )

        plt.savefig(out_file, dpi=300)
        plt.close()

        print(f"Saved and overwritten: {out_file}")
