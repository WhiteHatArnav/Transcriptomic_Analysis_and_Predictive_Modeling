import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

# ======================================================
# OUTPUT DIRECTORIES
# ======================================================

OUT_DIRS = {
    "KFold": "/Users/arnavjoshi/Desktop/AnalyticalKMCurves/LUADRecurrence/KFold",
    "SeventyThirty": "/Users/arnavjoshi/Desktop/AnalyticalKMCurves/LUADRecurrence/SeventyThirty",
}

# ======================================================
# INPUT ROOT
# ======================================================

SURVIVAL_ROOT = "/Users/arnavjoshi/Desktop/LUADPredictiveAnalyticsRecurrence"

# ---------------- KFold Risk Score Paths ----------------
KFOLD_RS_FILES = {
    "White": "/Users/arnavjoshi/Desktop/KFoldLassoCox/LUADKFoldLassoCoxRecurrence/RiskScores_LUAD_White_DFI_KFold10.csv",
    "All":   "/Users/arnavjoshi/Desktop/KFoldLassoCox/LUADKFoldLassoCoxRecurrence/RiskScores_LUAD_All_DFI_KFold10.csv",
    "Black": "/Users/arnavjoshi/Desktop/KFoldLassoCox/LUADKFoldLassoCoxRecurrence/RiskScores_LUAD_Black_DFI_KFold10.csv",
}

# ---------------- 70:30 Risk Score Paths (FIXED) ----------------
SPLIT70_RS_FILES = {
    "White": "/Users/arnavjoshi/Desktop/LUADPredictiveAnalyticsRecurrence/RiskScores_LUAD_White_DFI_CV_7030.csv",
    "All":   "/Users/arnavjoshi/Desktop/LUADPredictiveAnalyticsRecurrence/RiskScores_LUAD_All_DFI_CV_7030.csv",
    "Black": "/Users/arnavjoshi/Desktop/LUADPredictiveAnalyticsRecurrence/RiskScores_LUAD_Black_DFI_CV_7030.csv",
}

# ======================================================
# COLUMN NAMES
# ======================================================

PATIENT_COL = "Patient ID"
EVENT_COL = "DFI"
TIME_COL = "DFI.time"
RISK_COL = "Risk Score"

# ======================================================
# HELPERS
# ======================================================

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def load_survival(race):
    path = os.path.join(
        SURVIVAL_ROOT,
        f"luad_{race.lower()}_DFI_survival_data.csv"
    )
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    df = pd.read_csv(path)
    df.rename(columns={df.columns[0]: PATIENT_COL}, inplace=True)
    df.columns = df.columns.str.strip()
    return df[[PATIENT_COL, EVENT_COL, TIME_COL]].copy()

def median_split(values):
    med = np.median(values)
    return values >= med, med

def km_plot(df, out_path, title):
    groups = df["RiskGroup"].unique()
    if len(groups) != 2:
        print(f"[SKIP] Invalid groups for {out_path}")
        return

    kmf = KaplanMeierFitter()
    plt.figure(figsize=(8, 6))
    ax = plt.gca()

    for g in ["Low Risk", "High Risk"]:
        sub = df[df["RiskGroup"] == g]
        kmf.fit(
            sub[TIME_COL],
            event_observed=sub[EVENT_COL],
            label=g
        )
        kmf.plot(ax=ax, ci_show=True)

    low = df[df["RiskGroup"] == "Low Risk"]
    high = df[df["RiskGroup"] == "High Risk"]

    lr = logrank_test(
        low[TIME_COL], high[TIME_COL],
        event_observed_A=low[EVENT_COL],
        event_observed_B=high[EVENT_COL]
    )

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    ax.text(
        xmin + 0.03 * (xmax - xmin),
        ymin + 0.06 * (ymax - ymin),
        f"log-rank p = {lr.p_value:.3e}",
        bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.9)
    )

    ax.set_title(title)
    ax.set_xlabel("Disease-Free Interval (days)")
    ax.set_ylabel("Disease-Free Probability")

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved: {out_path}")

# ======================================================
# CORE PROCESSOR
# ======================================================

def process_one(race, mode):
    ensure_dir(OUT_DIRS[mode])

    surv = load_survival(race)
    rs_path = KFOLD_RS_FILES[race] if mode == "KFold" else SPLIT70_RS_FILES[race]

    rs = pd.read_csv(rs_path)
    rs.rename(columns={rs.columns[0]: PATIENT_COL}, inplace=True)
    rs.columns = rs.columns.str.strip()

    merged = surv.merge(
        rs[[PATIENT_COL, RISK_COL]],
        on=PATIENT_COL,
        how="inner"
    )

    if merged.shape[0] < 5:
        print(f"[SKIP] Too few samples: LUAD {race} {mode}")
        return

    high_mask, _ = median_split(merged[RISK_COL].values)
    merged["RiskGroup"] = np.where(high_mask, "High Risk", "Low Risk")

    out_path = os.path.join(
        OUT_DIRS[mode],
        f"LUAD_{race}_KM_byRisk_{mode}_DFI.png"
    )

    km_plot(
        merged,
        out_path,
        title=f"LUAD {race} – DFI by Risk ({mode})"
    )

# ======================================================
# MAIN
# ======================================================

def main():
    for p in OUT_DIRS.values():
        ensure_dir(p)

    # -------- KFold --------
    for race in ["White", "All", "Black"]:
        try:
            process_one(race, "KFold")
        except Exception as e:
            print(f"[WARN] KFold {race}: {e}")

    # -------- 70:30 --------
    for race in ["White", "All", "Black"]:
        try:
            process_one(race, "SeventyThirty")
        except Exception as e:
            print(f"[WARN] 70:30 {race}: {e}")

if __name__ == "__main__":
    main()
