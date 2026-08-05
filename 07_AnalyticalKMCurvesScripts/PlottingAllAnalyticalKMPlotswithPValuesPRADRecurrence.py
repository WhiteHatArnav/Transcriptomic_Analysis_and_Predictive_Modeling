import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

# ======================================================
# Configuration
# ======================================================

OUT_DIRS = {
    "KFold": "/Users/arnavjoshi/Desktop/AnalyticalKMCurves/KFold",
    "SeventyThirty": "/Users/arnavjoshi/Desktop/AnalyticalKMCurves/SeventyThirtySplit",
}

# Expression survival data roots (PRAD original survival files)
SURVIVAL_ROOT_PRAD = "/Users/arnavjoshi/Desktop/PRADPredictiveAnalytics"

# Recurrence table
RECURRENCE_PATH = "/Users/arnavjoshi/Desktop/PRADOtherData/PRADReccurenceData.tsv"

# Races
KFOLD_RACES = ["All", "White", "Black"]
SPLIT70_RACES = ["All", "White"]   # BLACK DOES NOT EXIST

# ---------------- KFold Risk Score Paths ----------------
KFOLD_RS_FILES = {
    "All":   "/Users/arnavjoshi/Desktop/KFoldLassoCox/PRADKFoldLassoCoxRecurrence/PRAD_All_KFold10_RiskScores_Recurrence.csv",
    "White": "/Users/arnavjoshi/Desktop/KFoldLassoCox/PRADKFoldLassoCoxRecurrence/PRAD_White_KFold10_RiskScores_Recurrence.csv",
    "Black": "/Users/arnavjoshi/Desktop/KFoldLassoCox/PRADKFoldLassoCoxRecurrence/PRAD_Black_KFold10_RiskScores_Recurrence_OSCensoring.csv"
}

# ---------------- 70/30 Risk Score Paths ----------------
SPLIT70_RS_FILES = {
    "All":   "/Users/arnavjoshi/Desktop/PRADPredictiveAnalyticsRecurrence/RiskScores_All_Recurrence_CV.csv",
    "White": "/Users/arnavjoshi/Desktop/PRADPredictiveAnalyticsRecurrence/RiskScores_White_Recurrence_CV.csv"
    # BLACK DOES NOT EXIST
}

# Columns
PATIENT_COL = "Patient ID"
REC_EVENT = "biochemical_recurrence"
REC_TIME = "days_to_first_biochemical_recurrence"

# NOTE for Black KFold only: OS fallback
OS_EVENT = "OS"
OS_TIME = "OS.time"


# ======================================================
# Helper functions
# ======================================================

def ensure_dir(p):
    os.makedirs(p, exist_ok=True)

def load_expression_survival(race):
    """
    Loads prad_<race>_survival_data.csv
    """
    fname = f"prad_{race.lower()}_survival_data.csv"
    path = os.path.join(SURVIVAL_ROOT_PRAD, fname)

    if not os.path.exists(path):
        raise FileNotFoundError(f"Survival file not found: {path}")

    df = pd.read_csv(path)
    df.rename(columns={df.columns[0]: PATIENT_COL}, inplace=True)
    df.columns = df.columns.str.strip()
    return df


def load_recurrence_table():
    rec = pd.read_csv(RECURRENCE_PATH, sep="\t")
    rec.rename(columns={rec.columns[0]: PATIENT_COL}, inplace=True)
    rec.columns = rec.columns.str.strip()

    # Normalize YES/NO -> 1/0
    rec[REC_EVENT] = rec[REC_EVENT].astype(str).str.upper().map({"YES": 1, "NO": 0})
    return rec[[PATIENT_COL, REC_EVENT, REC_TIME]].copy()


def merge_survival_and_recurrence(expr_df, recurrence_df, race, mode):
    """
    Returns dataframe with:
       PATIENT_ID, REC_EVENT, REC_TIME
    For Black KFold only → fill recurrence time from OS.time if missing.
    """
    merged = expr_df.merge(recurrence_df, on=PATIENT_COL, how="left")

    # Black KFold special case: fill recurrence time with OS.time
    if race == "Black" and mode == "KFold":
        missing = merged[REC_TIME].isna()

        if OS_TIME in merged.columns:
            merged.loc[missing, REC_TIME] = merged.loc[missing, OS_TIME]
            merged.loc[missing, REC_EVENT] = merged.loc[missing, OS_EVENT]
        else:
            print("[WARN] Black KFold OS fallback requested but OS columns missing.")

    # Drop rows missing recurrence time or event
    merged = merged.dropna(subset=[REC_TIME, REC_EVENT]).copy()
    merged[REC_TIME] = pd.to_numeric(merged[REC_TIME], errors="coerce")
    merged = merged.dropna(subset=[REC_TIME])

    return merged[[PATIENT_COL, REC_EVENT, REC_TIME]].copy()


def km_plot(df, group_col, out_path, title):
    groups = sorted(df[group_col].dropna().unique())

    if len(groups) != 2:
        print(f"Skipping KM plot at {out_path}, groups={groups}")
        return

    kmf = KaplanMeierFitter()
    plt.figure(figsize=(7.8, 6.2))
    ax = plt.gca()

    # Plot curves
    for g in groups:
        sub = df[df[group_col] == g]
        kmf.fit(durations=sub[REC_TIME].values,
                event_observed=sub[REC_EVENT].values,
                label=str(g))
        kmf.plot(ci_show=True, ax=ax)

    # Legend
    leg = ax.legend(loc="upper right")
    if leg is not None and leg.get_frame() is not None:
        leg.get_frame().set_alpha(0.85)

    # Log-rank test
    sub1 = df[df[group_col] == groups[0]]
    sub2 = df[df[group_col] == groups[1]]

    try:
        lr = logrank_test(
            sub1[REC_TIME], sub2[REC_TIME],
            event_observed_A=sub1[REC_EVENT],
            event_observed_B=sub2[REC_EVENT]
        )
        pval = lr.p_value

        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        ax.text(
            xmin + 0.03 * (xmax - xmin),
            ymin + 0.06 * (ymax - ymin),
            f"log-rank p = {pval:.3e}",
            bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.9)
        )
    except Exception as e:
        print(f"[WARN] Log-rank failed for {out_path}: {e}")

    ax.set_title(title)
    ax.set_xlabel("Recurrence-Free Time")
    ax.set_ylabel("Recurrence-Free Probability")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved KM curve: {out_path}")


def median_split(series):
    med = np.median(series)
    return (series >= med), med


# ======================================================
# Processing Function
# ======================================================

def process_one(race, mode):
    """
    mode ∈ {"KFold", "SeventyThirty"}
    """
    out_root = OUT_DIRS[mode]
    ensure_dir(out_root)

    expr = load_expression_survival(race)
    rec = load_recurrence_table()
    merged = merge_survival_and_recurrence(expr, rec, race, mode)

    # Load risk scores
    if mode == "KFold":
        rs_path = KFOLD_RS_FILES[race]
    else:
        if race not in SPLIT70_RACES:
            print(f"Skipping {race} 70/30 (file not available)")
            return
        rs_path = SPLIT70_RS_FILES[race]

    if not os.path.exists(rs_path):
        print(f"[WARN] Risk score file missing: {rs_path}")
        return

    rs_df = pd.read_csv(rs_path)
    rs_df.rename(columns={rs_df.columns[0]: PATIENT_COL}, inplace=True)
    rs_df.columns = rs_df.columns.str.strip()

    # Determine correct risk column
    rcol = "RiskScore_OutOfFold" if mode == "KFold" else "Risk Score"
    if rcol not in rs_df.columns:
        raise ValueError(f"{rcol} missing in {rs_path}")

    merged = merged.merge(rs_df[[PATIENT_COL, rcol]], on=PATIENT_COL, how="inner")

    if merged.shape[0] < 2:
        print(f"Too few samples for {race}-{mode}")
        return

    # Risk-group split
    high_mask, risk_med = median_split(merged[rcol].values)
    merged["RiskGroup"] = np.where(high_mask, "High Risk", "Low Risk")

    # Time split
    short_mask, time_med = median_split(merged[REC_TIME].values)
    merged["TimeGroup"] = np.where(short_mask,
                                   f"Short (≤ median {time_med:.1f})",
                                   f"Long (> median {time_med:.1f})")

    # Output paths
    risk_out = os.path.join(out_root, f"PRAD_{race}_KM_byRisk_{mode}_Recurrence.png")
    time_out = os.path.join(out_root, f"PRAD_{race}_KM_byTime_{mode}_Recurrence.png")

    km_plot(
        merged,
        group_col="RiskGroup",
        out_path=risk_out,
        title=f"PRAD {race} - Recurrence KM by Risk (median={risk_med:.3g}) [{mode}]"
    )

    km_plot(
        merged,
        group_col="TimeGroup",
        out_path=time_out,
        title=f"PRAD {race} - Recurrence KM by Time (median={time_med:.1f}) [{mode}]"
    )


# ======================================================
# Main Runner
# ======================================================

def main():
    ensure_dir(OUT_DIRS["KFold"])
    ensure_dir(OUT_DIRS["SeventyThirty"])

    # ---- KFold all races ----
    for race in KFOLD_RACES:
        try:
            process_one(race, mode="KFold")
        except Exception as e:
            print(f"[WARN] KFold {race}: {e}")

    # ---- 70/30 only All + White ----
    for race in SPLIT70_RACES:
        try:
            process_one(race, mode="SeventyThirty")
        except Exception as e:
            print(f"[WARN] 70/30 {race}: {e}")


if __name__ == "__main__":
    main()
