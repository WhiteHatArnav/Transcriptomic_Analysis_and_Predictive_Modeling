import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import plotly.graph_objs as go
from plotly.offline import plot

BLACK_DIR = "/Users/arnavjoshi/Desktop/BlackSamplesPCA/PRAD_PCA"
WHITE_DIR = "/Users/arnavjoshi/Desktop/WhiteSamplesPCA/PRAD_PCA"
COMB_DIR  = "/Users/arnavjoshi/Desktop/PCAWhiteBlackCombined/PRAD"
SCORES_FILE = "PRAD_per_sample_PCA_scores.csv"
BLACK_SCORES = os.path.join(BLACK_DIR, SCORES_FILE)
WHITE_SCORES = os.path.join(WHITE_DIR, SCORES_FILE)

BLACK_INPUT = "/Users/arnavjoshi/Desktop/PRADPredictiveAnalytics/prad_black_survival_data.csv"
WHITE_INPUT = "/Users/arnavjoshi/Desktop/PRADPredictiveAnalytics/prad_white_survival_data.csv"

def safe_load_scores(path):
    df = pd.read_csv(path)
    for col in ["PC1","PC2","PC3","OS"]:
        if col not in df.columns:
            raise ValueError(f"{col} missing in {path}")
    df["OS"] = pd.to_numeric(df["OS"], errors="coerce").fillna(0).astype(int)
    return df

def load_expr(path):
    df = pd.read_csv(path)
    if df.columns[0] != "sample" and df.columns[0] != "Patient ID":
        df = df.rename(columns={df.columns[0]: "sample"})
    df.columns = df.columns.str.strip()
    if "OS" not in df.columns or "OS.time" not in df.columns:
        raise ValueError(f"Missing OS/OS.time in {path}")
    df = df.dropna(subset=["OS","OS.time"]).copy()
    df["OS"] = pd.to_numeric(df["OS"], errors="coerce").fillna(0).astype(int)
    return df

def plot_interactive(df_map, title, outpath):
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
    for label, df in df_map.items():
        if df.empty:
            continue
        traces.append(go.Scatter3d(
            x=df["PC1"], y=df["PC2"], z=df["PC3"],
            mode="markers",
            name=label,
            marker=dict(
                size=5, color=color_map[label], symbol=symbol_map[label],
                opacity=0.9, line=dict(width=0.5, color="black")
            )
        ))
    layout = go.Layout(
        title=title,
        scene=dict(xaxis_title="PC1", yaxis_title="PC2", zaxis_title="PC3"),
        margin=dict(l=0, r=0, b=0, t=30),
        legend=dict(x=0.02, y=0.98)
    )
    fig = go.Figure(data=traces, layout=layout)
    plot(fig, filename=outpath, auto_open=False)
    print(f"Saved: {outpath}")

def separate_per_race():
    os.makedirs(BLACK_DIR, exist_ok=True)
    os.makedirs(WHITE_DIR, exist_ok=True)
    bd = safe_load_scores(BLACK_SCORES)
    wd = safe_load_scores(WHITE_SCORES)

    plot_interactive({
        "Black, OS=0": bd[bd["OS"]==0],
        "Black, OS=1": bd[bd["OS"]==1],
    }, "PRAD Black PCA (3D)", os.path.join(BLACK_DIR, "PRAD_Black_3D_PCA.html"))

    plot_interactive({
        "White, OS=0": wd[wd["OS"]==0],
        "White, OS=1": wd[wd["OS"]==1],
    }, "PRAD White PCA (3D)", os.path.join(WHITE_DIR, "PRAD_White_3D_PCA.html"))

def overlay_separate_pca():
    os.makedirs(COMB_DIR, exist_ok=True)
    bd = safe_load_scores(BLACK_SCORES)
    wd = safe_load_scores(WHITE_SCORES)
    plot_interactive({
        "White, OS=0": wd[wd["OS"]==0],
        "White, OS=1": wd[wd["OS"]==1],
        "Black, OS=0": bd[bd["OS"]==0],
        "Black, OS=1": bd[bd["OS"]==1],
    }, "PRAD: Overlay of Separate PCAs (3D)",
       os.path.join(COMB_DIR, "PRAD_WhiteBlack_SeparatePCA_3D.html"))

def single_pca_combined():
    os.makedirs(COMB_DIR, exist_ok=True)
    bd = load_expr(BLACK_INPUT)
    wd = load_expr(WHITE_INPUT)

    gene_cols = [c for c in bd.columns if c not in ["sample","Patient ID","OS","OS.time"]]
    gene_cols = [c for c in gene_cols if c in wd.columns]

    bd["Race"]="Black"; wd["Race"]="White"
    comb = pd.concat([bd,wd], ignore_index=True)
    X = comb[gene_cols].values
    Xs = StandardScaler().fit_transform(X)
    pca = PCA(n_components=3, random_state=42)
    pcs = pca.fit_transform(Xs)
    comb["PC1"], comb["PC2"], comb["PC3"] = pcs[:,0], pcs[:,1], pcs[:,2]

    plot_interactive({
        "White, OS=0": comb[(comb["Race"]=="White") & (comb["OS"]==0)],
        "White, OS=1": comb[(comb["Race"]=="White") & (comb["OS"]==1)],
        "Black, OS=0": comb[(comb["Race"]=="Black") & (comb["OS"]==0)],
        "Black, OS=1": comb[(comb["Race"]=="Black") & (comb["OS"]==1)],
    }, "PRAD: Single PCA (Combined, 3D)",
       os.path.join(COMB_DIR, "PRAD_WhiteBlack_singlePCA_3D.html"))

def main():
    separate_per_race()
    overlay_separate_pca()
    single_pca_combined()

if __name__ == "__main__":
    main()
