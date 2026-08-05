import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test  # added for p-value

# -----------------------------
# Configuration
# -----------------------------
OUT_DIRS = {
    "KFold": "/Users/arnavjoshi/Desktop/AnalyticalKMCurves/KFold",
    "SeventyThirty": "/Users/arnavjoshi/Desktop/AnalyticalKMCurves/SeventyThirtySplit",
}

# Survival data roots
SURVIVAL_ROOTS = {
    "BRCA": "/Users/arnavjoshi/Desktop/BRCAPredictiveAnalytics",
    "LIHC": "/Users/arnavjoshi/Desktop/LIHCPredictiveAnalytics",
    "LUAD": "/Users/arnavjoshi/Desktop/LUADPredictiveAnalytics",
}

RACES = ["Black", "White", "All"]

KFOLD_RS_ROOTS = {
    "BRCA": "/Users/arnavjoshi/Desktop/KFoldLassoCox/BRCAKFoldLassoCox",
    "LIHC": "/Users/arnavjoshi/Desktop/KFoldLassoCox/LIHCKFoldLassoCox",
    "LUAD": "/Users/arnavjoshi/Desktop/KFoldLassoCox/LUADKFoldLassoCox",
}
KFOLD_RS_PATTERN = "{cancer}_{race}_KFold10_RiskScores.csv"

SPLIT70_RS_ROOTS = {
    "BRCA": "/Users/arnavjoshi/Desktop/BRCAPredictiveAnalytics",
    "LIHC": "/Users/arnavjoshi/Desktop/LIHCPredictiveAnalytics",
    "LUAD": "/Users/arnavjoshi/Desktop/LUADPredictiveAnalytics",
}
SPLIT70_RS_PATTERN = "RiskScores_{race}_CV.csv"

PATIENT_COL = "Patient ID"
EVENT_COL = "OS"
TIME_COL = "OS.time"

# -----------------------------
# Helpers
# -----------------------------
def ensure_dir(p):
    os.makedirs(p, exist_ok=True)

def load_survival(cancer, race):
    root = SURVIVAL_ROOTS[cancer]
    fname = f"{cancer.lower()}_{race.lower()}_survival_data.csv"
    path = os.path.join(root, fname)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Survival file not found: {path}")
    df = pd.read_csv(path)
    if df.columns[0] != PATIENT_COL:
        df = df.rename(columns={df.columns[0]: PATIENT_COL})
    df.columns = [c.strip() for c in df.columns]
    df = df.dropna(subset=[EVENT_COL, TIME_COL]).copy()
    return df

def load_kfold_risk(cancer, race):
    root = KFOLD_RS_ROOTS[cancer]
    fname = KFOLD_RS_PATTERN.format(cancer=cancer, race=race)
    path = os.path.join(root, fname)
    if not os.path.exists(path):
        raise FileNotFoundError(f"KFold risk file not found: {path}")
    df = pd.read_csv(path)
    if df.columns[0] != PATIENT_COL:
        df = df.rename(columns={df.columns[0]: PATIENT_COL})
    df.columns = [c.strip() for c in df.columns]
    if "RiskScore_OutOfFold" not in df.columns:
        raise ValueError(f"'RiskScore_OutOfFold' column missing in {path}")
    return df[[PATIENT_COL, "RiskScore_OutOfFold"]].copy()

def load_split70_risk(cancer, race):
    root = SPLIT70_RS_ROOTS[cancer]
    fname = SPLIT70_RS_PATTERN.format(race=race)
    path = os.path.join(root, fname)
    if not os.path.exists(path):
        raise FileNotFoundError(f"70-30 risk file not found: {path}")
    df = pd.read_csv(path)
    if df.columns[0] != PATIENT_COL:
        df = df.rename(columns={df.columns[0]: PATIENT_COL})
    df.columns = [c.strip() for c in df.columns]
    if "Risk Score" not in df.columns:
        raise ValueError(f"'Risk Score' column missing in {path}")
    return df[[PATIENT_COL, "Risk Score"]].copy()

def median_split(series):
    med = np.median(series)
    return series >= med, med

def km_plot(df, group_col, out_path, title):
    """Plot KM curves for two groups defined in df[group_col], with log-rank p-value shown bottom-left."""
    groups = sorted(df[group_col].dropna().unique())
    if len(groups) != 2:
        print(f"Skipping KM plot (need 2 groups) at {out_path}; found groups: {groups}")
        return

    kmf = KaplanMeierFitter()
    plt.figure(figsize=(7.5, 6.0))

    for g in groups:
        sub = df[df[group_col] == g]
        if sub.empty:
            continue
        kmf.fit(durations=sub[TIME_COL].values, event_observed=sub[EVENT_COL].values, label=str(g))
        kmf.plot(ci_show=True)

    # Compute log-rank p-value
    g1, g2 = groups[0], groups[1]
    sub1 = df[df[group_col] == g1]
    sub2 = df[df[group_col] == g2]
    try:
        lr = logrank_test(
            sub1[TIME_COL].values, sub2[TIME_COL].values,
            event_observed_A=sub1[EVENT_COL].values,
            event_observed_B=sub2[EVENT_COL].values
        )
        pval = lr.p_value
        ax = plt.gca()
        xmin = ax.get_xlim()[0]
        ymin = ax.get_ylim()[0]
        ax.text(
            xmin + 0.02*(ax.get_xlim()[1]-ax.get_xlim()[0]),
            ymin + 0.05*(ax.get_ylim()[1]-ax.get_ylim()[0]),
            f"log-rank p = {pval:.3e}",
            ha="left", va="bottom", fontsize=10,
            bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.8)
        )
    except Exception as e:
        print(f"[WARN] Failed to compute log-rank p-value for {out_path}: {e}")

    plt.title(title)
    plt.xlabel("Time")
    plt.ylabel("Survival probability")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved KM: {out_path}")

def process_one(cancer, race, mode):
    out_root = OUT_DIRS["KFold"] if mode == "KFold" else OUT_DIRS["SeventyThirty"]
    ensure_dir(out_root)
    surv_df = load_survival(cancer, race)
    if mode == "KFold":
        rs_df = load_kfold_risk(cancer, race)
        rs_col = "RiskScore_OutOfFold"
    else:
        rs_df = load_split70_risk(cancer, race)
        rs_col = "Risk Score"

    merged = surv_df.merge(rs_df, on=PATIENT_COL, how="inner").copy()
    if merged.shape[0] < 2:
        print(f"Too few samples after merge for {cancer}-{race}-{mode}; n={merged.shape[0]}")
        return

    high_mask, risk_med = median_split(merged[rs_col].values)
    merged["RiskGroup"] = np.where(high_mask, "High Risk", "Low Risk")

    short_mask, time_med = median_split(merged[TIME_COL].values)
    merged["TimeGroup"] = np.where(short_mask, f"Short (≤ median {time_med:.1f})", f"Long (> median {time_med:.1f})")

    risk_out = os.path.join(out_root, f"{cancer}_{race}_KM_byRisk_{mode}.png")
    time_out = os.path.join(out_root, f"{cancer}_{race}_KM_bySurvivalTime_{mode}.png")

    km_plot(
        merged,
        group_col="RiskGroup",
        out_path=risk_out,
        title=f"{cancer} {race} – KM by Risk (median split at {risk_med:.4g}) [{mode}]"
    )
    km_plot(
        merged,
        group_col="TimeGroup",
        out_path=time_out,
        title=f"{cancer} {race} – KM by Survival Time (median {time_med:.1f}) [{mode}]"
    )

def main():
    for mode in ["KFold", "SeventyThirty"]:
        ensure_dir(OUT_DIRS[mode])

    cancers = ["BRCA", "LIHC", "LUAD"]
    for cancer in cancers:
        for race in RACES:
            try:
                process_one(cancer, race, mode="KFold")
            except Exception as e:
                print(f"[WARN] KFold {cancer}-{race}: {e}")
            try:
                process_one(cancer, race, mode="SeventyThirty")
            except Exception as e:
                print(f"[WARN] SeventyThirty {cancer}-{race}: {e}")

if __name__ == "__main__":
    main()
