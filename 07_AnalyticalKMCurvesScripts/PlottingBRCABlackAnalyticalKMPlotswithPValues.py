import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

# -----------------------------
# Configuration (BRCA Black – 70/30 ONLY)
# -----------------------------
OUT_DIR = "/Users/arnavjoshi/Desktop/AnalyticalKMCurves/SeventyThirtySplit"
SURVIVAL_ROOT = "/Users/arnavjoshi/Desktop/BRCAPredictiveAnalytics"
RISK_ROOT = "/Users/arnavjoshi/Desktop/BRCAPredictiveAnalytics"

PATIENT_COL = "Patient ID"
EVENT_COL = "OS"
TIME_COL = "OS.time"

os.makedirs(OUT_DIR, exist_ok=True)

# -----------------------------
# Load survival data
# -----------------------------
survival_file = os.path.join(
    SURVIVAL_ROOT,
    "brca_black_survival_data.csv"
)

surv_df = pd.read_csv(survival_file)
if surv_df.columns[0] != PATIENT_COL:
    surv_df = surv_df.rename(columns={surv_df.columns[0]: PATIENT_COL})

surv_df.columns = [c.strip() for c in surv_df.columns]
surv_df = surv_df.dropna(subset=[EVENT_COL, TIME_COL]).copy()

# -----------------------------
# Load 70/30 risk scores
# -----------------------------
risk_file = os.path.join(
    RISK_ROOT,
    "RiskScores_Black_CV.csv"
)

risk_df = pd.read_csv(risk_file)
if risk_df.columns[0] != PATIENT_COL:
    risk_df = risk_df.rename(columns={risk_df.columns[0]: PATIENT_COL})

risk_df.columns = [c.strip() for c in risk_df.columns]
risk_df = risk_df[[PATIENT_COL, "Risk Score"]].copy()

# -----------------------------
# Merge
# -----------------------------
merged = surv_df.merge(risk_df, on=PATIENT_COL, how="inner")
if merged.shape[0] < 2:
    raise ValueError("Too few samples after merge for BRCA Black 70/30 split")

# -----------------------------
# Median splits
# -----------------------------
risk_median = np.median(merged["Risk Score"])
merged["RiskGroup"] = np.where(
    merged["Risk Score"] >= risk_median,
    "High Risk",
    "Low Risk"
)

time_median = np.median(merged[TIME_COL])
merged["TimeGroup"] = np.where(
    merged[TIME_COL] >= time_median,
    f"Short (≤ median {time_median:.1f})",
    f"Long (> median {time_median:.1f})"
)

# -----------------------------
# KM plot helper
# -----------------------------
def km_plot(df, group_col, out_path, title):
    groups = sorted(df[group_col].unique())
    if len(groups) != 2:
        print(f"Skipping plot, found groups: {groups}")
        return

    kmf = KaplanMeierFitter()
    plt.figure(figsize=(7.5, 6.0))

    for g in groups:
        sub = df[df[group_col] == g]
        kmf.fit(
            durations=sub[TIME_COL],
            event_observed=sub[EVENT_COL],
            label=g
        )
        kmf.plot(ci_show=True)

    g1, g2 = groups
    sub1 = df[df[group_col] == g1]
    sub2 = df[df[group_col] == g2]

    lr = logrank_test(
        sub1[TIME_COL], sub2[TIME_COL],
        event_observed_A=sub1[EVENT_COL],
        event_observed_B=sub2[EVENT_COL]
    )

    ax = plt.gca()
    ax.text(
        0.02, 0.05,
        f"log-rank p = {lr.p_value:.3e}",
        transform=ax.transAxes,
        bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.8)
    )

    plt.title(title)
    plt.xlabel("Time")
    plt.ylabel("Survival probability")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved KM: {out_path}")

# -----------------------------
# Plots
# -----------------------------
km_plot(
    merged,
    "RiskGroup",
    os.path.join(OUT_DIR, "BRCA_Black_KM_byRisk_SeventyThirty.png"),
    f"BRCA Black – KM by Risk (median {risk_median:.4g}) [70/30]"
)

km_plot(
    merged,
    "TimeGroup",
    os.path.join(OUT_DIR, "BRCA_Black_KM_bySurvivalTime_SeventyThirty.png"),
    f"BRCA Black – KM by Survival Time (median {time_median:.1f}) [70/30]"
)
