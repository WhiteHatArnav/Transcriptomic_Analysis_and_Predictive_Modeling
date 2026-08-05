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

OUT_DIR = "/Users/arnavjoshi/Desktop/AnalyticalKMCurves/KFold"

SURVIVAL_ROOT = "/Users/arnavjoshi/Desktop/PRADPredictiveAnalytics"
RECURRENCE_PATH = "/Users/arnavjoshi/Desktop/PRADOtherData/PRADReccurenceData.tsv"

KFOLD_RS_FILES = {
    "All":   "/Users/arnavjoshi/Desktop/KFoldLassoCox/PRADKFoldLassoCoxRecurrence/PRAD_All_KFold10_RiskScores_Recurrence.csv",
    "White": "/Users/arnavjoshi/Desktop/KFoldLassoCox/PRADKFoldLassoCoxRecurrence/PRAD_White_KFold10_RiskScores_Recurrence.csv",
    "Black": "/Users/arnavjoshi/Desktop/KFoldLassoCox/PRADKFoldLassoCoxRecurrence/PRAD_Black_KFold10_RiskScores_Recurrence_OSCensoring.csv"
}

RACES = ["All", "White", "Black"]

PATIENT_COL = "Patient ID"
REC_EVENT = "biochemical_recurrence"
REC_TIME = "days_to_first_biochemical_recurrence"
OS_EVENT = "OS"
OS_TIME = "OS.time"

# ======================================================
# Helpers
# ======================================================

def ensure_dir(p):
    os.makedirs(p, exist_ok=True)

def load_expression_survival(race):
    fname = f"prad_{race.lower()}_survival_data.csv"
    path = os.path.join(SURVIVAL_ROOT, fname)
    df = pd.read_csv(path)
    df.rename(columns={df.columns[0]: PATIENT_COL}, inplace=True)
    df.columns = df.columns.str.strip()
    return df

def load_recurrence():
    rec = pd.read_csv(RECURRENCE_PATH, sep="\t")
    rec.rename(columns={rec.columns[0]: PATIENT_COL}, inplace=True)
    rec.columns = rec.columns.str.strip()
    rec[REC_EVENT] = rec[REC_EVENT].astype(str).str.upper().map({"YES": 1, "NO": 0})
    return rec[[PATIENT_COL, REC_EVENT, REC_TIME]].copy()

def median_split(series):
    med = np.median(series)
    return (series >= med), med

# ======================================================
# KM Plotting
# ======================================================

def km_plot(df, group_col, out_path, title):
    """Full-time KM."""
    groups = sorted(df[group_col].dropna().unique())
    if len(groups) != 2:
        print(f"[WARN] Skip KM: need 2 groups, found {groups}")
        return

    kmf = KaplanMeierFitter()
    plt.figure(figsize=(7.8, 6.2))
    ax = plt.gca()

    for g in groups:
        sub = df[df[group_col] == g]
        kmf.fit(sub[REC_TIME], event_observed=sub[REC_EVENT], label=str(g))
        kmf.plot(ci_show=True, ax=ax)

    leg = ax.legend(loc="upper right")
    if leg is not None and leg.get_frame() is not None:
        leg.get_frame().set_alpha(0.85)

    g1, g2 = groups
    sub1 = df[df[group_col] == g1]
    sub2 = df[df[group_col] == g2]

    # Log-rank test
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
    except:
        pass

    ax.set_title(title)
    ax.set_xlabel("Recurrence-Free Time")
    ax.set_ylabel("Recurrence-Free Probability")

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved KM: {out_path}")


def km_plot_time_constrained(df, group_col, out_path, title):
    """
    TRUE 250–750 KM: FILTER DATA BEFORE FITTING.
    """

    # Filter samples inside 250–750 day interval
    df_tc = df[(df[REC_TIME] >= 250) & (df[REC_TIME] <= 750)].copy()

    if df_tc.shape[0] < 2:
        print(f"[WARN] Too few filtered samples for time-constrained KM: {out_path}")
        return

    groups = sorted(df_tc[group_col].dropna().unique())
    if len(groups) != 2:
        print(f"[WARN] Need 2 groups for time-constrained KM: {groups}")
        return

    kmf = KaplanMeierFitter()
    plt.figure(figsize=(7.8, 6.2))
    ax = plt.gca()

    # Fit using filtered data
    for g in groups:
        sub = df_tc[df_tc[group_col] == g]
        kmf.fit(sub[REC_TIME], event_observed=sub[REC_EVENT], label=str(g))
        kmf.plot(ci_show=True, ax=ax)

    leg = ax.legend(loc="upper right")
    if leg is not None and leg.get_frame() is not None:
        leg.get_frame().set_alpha(0.85)

    g1, g2 = groups
    sub1 = df_tc[df_tc[group_col] == g1]
    sub2 = df_tc[df_tc[group_col] == g2]

    # Log-rank test using only filtered data
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
            250 + 0.03 * (750 - 250),
            ymin + 0.06 * (ymax - ymin),
            f"log-rank p (250–750d) = {pval:.3e}",
            bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.9)
        )
    except:
        pass

    ax.set_title(title + " (250–750 days, filtered)")
    ax.set_xlabel("Recurrence-Free Time")
    ax.set_ylabel("Recurrence-Free Probability")
    ax.set_xlim(250, 750)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved filtered KM: {out_path}")

# ======================================================
# Main Processing
# ======================================================

def process_race(race):
    expr = load_expression_survival(race)
    rec = load_recurrence()

    merged = expr.merge(rec, on=PATIENT_COL, how="left")

    if race == "Black":
        missing = merged[REC_TIME].isna()
        merged.loc[missing, REC_TIME] = merged.loc[missing, OS_TIME]
        merged.loc[missing, REC_EVENT] = merged.loc[missing, OS_EVENT]

    merged = merged.dropna(subset=[REC_TIME, REC_EVENT]).copy()
    merged[REC_TIME] = pd.to_numeric(merged[REC_TIME], errors="coerce")
    merged = merged.dropna(subset=[REC_TIME])

    # Load risk scores
    rs_path = KFOLD_RS_FILES[race]
    rs = pd.read_csv(rs_path)
    rs.rename(columns={rs.columns[0]: PATIENT_COL}, inplace=True)
    rs.columns = rs.columns.str.strip()

    merged = merged.merge(rs[[PATIENT_COL, "RiskScore_OutOfFold"]], on=PATIENT_COL, how="inner")

    # Median split
    high_mask, med = median_split(merged["RiskScore_OutOfFold"].values)
    merged["RiskGroup"] = np.where(high_mask, "High Risk", "Low Risk")

    # Full KM curve
    full_out = os.path.join(OUT_DIR, f"PRAD_{race}_KM_byRisk_KFold_Recurrence.png")
    km_plot(
        merged,
        group_col="RiskGroup",
        out_path=full_out,
        title=f"PRAD {race} - Recurrence KM by Risk (median={med:.3g}) [KFold]"
    )

    # TRUE FILTERED KM curve
    tc_out = os.path.join(OUT_DIR, f"PRAD_{race}_KM_byRisk_KFold_Recurrence_time_constrained.png")
    km_plot_time_constrained(
        merged,
        group_col="RiskGroup",
        out_path=tc_out,
        title=f"PRAD {race} - Recurrence KM by Risk [KFold]"
    )


def main():
    ensure_dir(OUT_DIR)
    for race in RACES:
        try:
            process_race(race)
        except Exception as e:
            print(f"[ERROR] Failed for {race}: {e}")


if __name__ == "__main__":
    main()
