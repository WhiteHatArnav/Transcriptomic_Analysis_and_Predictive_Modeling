import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import plotly.graph_objs as go
from plotly.offline import plot

RACE_PCA_DIRS = {
    "BRCA": {
        "black": "/Users/arnavjoshi/Desktop/BlackSamplesPCA/BRCA_PCA",
        "white": "/Users/arnavjoshi/Desktop/WhiteSamplesPCA/BRCA_PCA",
        "outfile_dir": "/Users/arnavjoshi/Desktop/PCAWhiteBlackCombined/BRCA"
    },
    "LIHC": {
        "black": "/Users/arnavjoshi/Desktop/BlackSamplesPCA/LIHC_PCA",
        "white": "/Users/arnavjoshi/Desktop/WhiteSamplesPCA/LIHC_PCA",
        "outfile_dir": "/Users/arnavjoshi/Desktop/PCAWhiteBlackCombined/LIHC"
    },
    "LUAD": {
        "black": "/Users/arnavjoshi/Desktop/BlackSamplesPCA/LUAD_PCA",
        "white": "/Users/arnavjoshi/Desktop/WhiteSamplesPCA/LUAD_PCA",
        "outfile_dir": "/Users/arnavjoshi/Desktop/PCAWhiteBlackCombined/LUAD"
    },
}

FILENAME = "{cancer}_per_sample_PCA_scores.csv"

def safe_load_scores(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing file: {path}")
    df = pd.read_csv(path)
    for col in ["PC1", "PC2", "PC3", "OS"]:
        if col not in df.columns:
            raise ValueError(f"{col} missing in {path}")
    df["OS"] = pd.to_numeric(df["OS"], errors="coerce").fillna(0).astype(int)
    return df

def load_expression_matrix(path):
    df = pd.read_csv(path)
    df.rename(columns={df.columns[0]: "Patient ID"}, inplace=True)
    df.columns = df.columns.str.strip()
    if "OS" not in df.columns or "OS.time" not in df.columns:
        raise ValueError(f"Missing OS/OS.time in {path}")
    df = df.dropna(subset=["OS", "OS.time"])
    df["OS"] = pd.to_numeric(df["OS"], errors="coerce").fillna(0).astype(int)
    return df

def plot_interactive_3D(df_dict, cancer, title, outpath):
    symbol_map = {
        "White, OS=0": "circle",
        "White, OS=1": "circle-open",
        "Black, OS=0": "square",
        "Black, OS=1": "diamond",
    }
    color_map = {
        "White, OS=0": "#1f77b4",
        "White, OS=1": "#1f77b4",
        "Black, OS=0": "#d62728",
        "Black, OS=1": "#d62728",
    }

    traces = []
    for label, df in df_dict.items():
        if df.empty:
            continue
        traces.append(
            go.Scatter3d(
                x=df["PC1"], y=df["PC2"], z=df["PC3"],
                mode="markers",
                name=label,
                marker=dict(
                    size=5,
                    color=color_map[label],
                    symbol=symbol_map[label],
                    opacity=0.9,
                    line=dict(width=0.5, color="black")
                )
            )
        )

    layout = go.Layout(
        title=title,
        scene=dict(xaxis_title="PC1", yaxis_title="PC2", zaxis_title="PC3"),
        margin=dict(l=0, r=0, b=0, t=30),
        legend=dict(x=0.02, y=0.98)
    )
    fig = go.Figure(data=traces, layout=layout)
    plot(fig, filename=outpath, auto_open=False)
    print(f"Saved interactive 3D plot: {outpath}")

def combined_single_pca_3D(black_path, white_path, cancer, outfile_dir):
    os.makedirs(outfile_dir, exist_ok=True)
    black_df = load_expression_matrix(black_path)
    white_df = load_expression_matrix(white_path)

    gene_cols = [c for c in black_df.columns if c not in ["Patient ID", "OS", "OS.time"]]
    gene_cols = [c for c in gene_cols if c in white_df.columns]

    black_df["Race"] = "Black"
    white_df["Race"] = "White"
    combined_df = pd.concat([black_df, white_df], ignore_index=True)

    X = combined_df[gene_cols].values
    X_scaled = StandardScaler().fit_transform(X)
    pca = PCA(n_components=3, random_state=42)
    pcs = pca.fit_transform(X_scaled)
    combined_df["PC1"], combined_df["PC2"], combined_df["PC3"] = pcs[:, 0], pcs[:, 1], pcs[:, 2]

    df_groups = {
        "White, OS=0": combined_df[(combined_df["Race"] == "White") & (combined_df["OS"] == 0)],
        "White, OS=1": combined_df[(combined_df["Race"] == "White") & (combined_df["OS"] == 1)],
        "Black, OS=0": combined_df[(combined_df["Race"] == "Black") & (combined_df["OS"] == 0)],
        "Black, OS=1": combined_df[(combined_df["Race"] == "Black") & (combined_df["OS"] == 1)],
    }

    outpath = os.path.join(outfile_dir, f"{cancer}_WhiteBlack_singlePCA_3D.html")
    plot_interactive_3D(df_groups, cancer, f"{cancer}: Single PCA (Combined)", outpath)

def overlay_separate_pca_3D(white_df, black_df, cancer, outfile_dir):
    os.makedirs(outfile_dir, exist_ok=True)
    df_groups = {
        "White, OS=0": white_df[white_df["OS"] == 0],
        "White, OS=1": white_df[white_df["OS"] == 1],
        "Black, OS=0": black_df[black_df["OS"] == 0],
        "Black, OS=1": black_df[black_df["OS"] == 1],
    }
    outpath = os.path.join(outfile_dir, f"{cancer}_WhiteBlack_SeparatePCA_3D.html")
    plot_interactive_3D(df_groups, cancer, f"{cancer}: Overlay of Separate PCAs", outpath)

def separate_race_plot_3D(df, cancer, race, outdir):
    os.makedirs(outdir, exist_ok=True)
    df_groups = {
        f"{race}, OS=0": df[df["OS"] == 0],
        f"{race}, OS=1": df[df["OS"] == 1],
    }
    # Remap to the 4-label scheme used by plot_interactive_3D for colors/symbols
    remapped = {}
    for label, sub in df_groups.items():
        if "White" in label:
            key = "White, OS=0" if "OS=0" in label else "White, OS=1"
        else:
            key = "Black, OS=0" if "OS=0" in label else "Black, OS=1"
        remapped[key] = sub

    outpath = os.path.join(outdir, f"{cancer}_{race}_3D_PCA.html")
    plot_interactive_3D(remapped, cancer, f"{cancer} {race} PCA (3D)", outpath)

def main():
    for cancer, d in RACE_PCA_DIRS.items():
        black_path = os.path.join(d["black"], FILENAME.format(cancer=cancer))
        white_path = os.path.join(d["white"], FILENAME.format(cancer=cancer))

        black_df = safe_load_scores(black_path)
        white_df = safe_load_scores(white_path)

        separate_race_plot_3D(black_df, cancer, "Black", d["black"])
        separate_race_plot_3D(white_df, cancer, "White", d["white"])
        overlay_separate_pca_3D(white_df, black_df, cancer, d["outfile_dir"])
        combined_single_pca_3D(black_path, white_path, cancer, d["outfile_dir"])

if __name__ == "__main__":
    main()
