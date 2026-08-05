import os
import numpy as np
import pandas as pd
import warnings
import multiprocessing as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LassoCV
from sklearn.model_selection import KFold
from sklearn.exceptions import ConvergenceWarning
from lifelines import CoxPHFitter, KaplanMeierFitter

# -----------------------------
# Paths (PRAD, Black cohort)
# -----------------------------
input_folder = "/Users/arnavjoshi/Desktop/PRADPredictiveAnalytics"
input_filename = "prad_black_survival_data.csv"
input_path = os.path.join(input_folder, input_filename)

recurrence_path = "/Users/arnavjoshi/Desktop/PRADOtherData/PRADReccurenceData.tsv"

output_folder = "/Users/arnavjoshi/Desktop/KFoldLassoCox/PRADKFoldLassoCoxRecurrence"
os.makedirs(output_folder, exist_ok=True)

# -----------------------------
# Logging helper
# -----------------------------
run_log_lines = []

def logprint(msg):
    print(msg, flush=True)
    run_log_lines.append(str(msg))

# -----------------------------
# Helper: run LASSO in a subprocess and capture ConvergenceWarning
# -----------------------------
def _fit_lasso_in_subprocess(X_train, y_train_event, cv_folds_for_lasso, return_dict):
    """
    Child process target. Fit LassoCV on (X_train, y_train_event).
    Store result in return_dict:
      return_dict["status"] in {"ok","warn","fail"}
      return_dict["coef"]   = numpy array of coefficients or None
    """
    status, coefs = "ok", None
    try:
        with warnings.catch_warnings(record=True) as wlist:
            warnings.simplefilter("always", ConvergenceWarning)
            lasso_cv = LassoCV(
                alphas=np.logspace(-4, -1, 20),
                cv=cv_folds_for_lasso,
                random_state=42,
                max_iter=2000,
                n_jobs=1  # avoid thread thrash inside subprocess
            )
            lasso_cv.fit(X_train, y_train_event)
            coefs = lasso_cv.coef_
            for w in wlist:
                if issubclass(w.category, ConvergenceWarning):
                    status = "warn"
                    break
    except Exception:
        status = "fail"
        coefs = None
    return_dict["status"] = status
    return_dict["coef"] = coefs

def fit_lasso_with_timeout(X_train, y_train_event, cv_folds_for_lasso, timeout_seconds=60):
    """
    Run LASSO in a child process with timeout.
    Returns (status, coefs), where status in {"ok","warn","fail","timeout"}.
    """
    manager = mp.Manager()
    shared = manager.dict()
    p = mp.Process(
        target=_fit_lasso_in_subprocess,
        args=(X_train, y_train_event, cv_folds_for_lasso, shared)
    )
    p.start()
    p.join(timeout_seconds)
    if p.is_alive():
        p.terminate()
        p.join()
        return ("timeout", None)
    return (shared.get("status", "fail"), shared.get("coef", None))

# -----------------------------
# Fallback LASSO (quick) used ONLY if timeout/fail
# -----------------------------
def fallback_fit_lasso_quick(X_train, y_train_event, cv_folds_for_lasso):
    try:
        lasso_cv = LassoCV(
            alphas=np.logspace(-3, -1, 12),
            cv=cv_folds_for_lasso,
            random_state=42,
            max_iter=3000,
            n_jobs=1
        )
        lasso_cv.fit(X_train, y_train_event)
        return lasso_cv.coef_
    except Exception:
        return None

# -----------------------------
# Main
# -----------------------------
def main():
    logprint(f"Loading {input_filename} ...")
    expr = pd.read_csv(input_path)
    expr.rename(columns={expr.columns[0]: "Patient ID"}, inplace=True)
    expr.columns = expr.columns.str.strip()

    patient_id_column = "Patient ID"
    os_time_col = "OS.time"

    # Load recurrence data
    logprint("Loading recurrence file ...")
    rec = pd.read_csv(recurrence_path, sep="\t")
    rec.rename(columns={rec.columns[0]: "Patient ID"}, inplace=True)
    rec.columns = rec.columns.str.strip()

    event_column = "biochemical_recurrence"
    raw_rec_time_col = "days_to_first_biochemical_recurrence"
    rec_time_col = "Recurrence_time_days"

    # Map YES/NO -> 1/0, preserve NaN initially
    rec[event_column] = rec[event_column].astype(str).str.strip().str.upper()
    rec[event_column] = rec[event_column].map({"YES": 1, "NO": 0})

    rec_subset = rec[[patient_id_column, event_column, raw_rec_time_col]]

    # Merge expression + recurrence
    data = expr.merge(rec_subset, on=patient_id_column, how="left")

    # Numeric OS.time and recurrence time
    if os_time_col not in data.columns:
        logprint(f"ERROR: {os_time_col} not found in {input_filename}.")
        raise SystemExit

    data.loc[:, os_time_col] = pd.to_numeric(data[os_time_col], errors="coerce")
    data.loc[:, raw_rec_time_col] = pd.to_numeric(data[raw_rec_time_col], errors="coerce")

    # Rule 3: missing biochemical_recurrence -> 0 (no documented recurrence)
    data.loc[:, event_column] = data[event_column].fillna(0).astype(int)

    # Start unified recurrence time from raw recurrence time
    data[rec_time_col] = data[raw_rec_time_col]

    # For rows where recurrence time is missing, fill from OS.time (censoring)
    mask_missing_rec_time = data[rec_time_col].isna()
    data.loc[mask_missing_rec_time, rec_time_col] = data.loc[mask_missing_rec_time, os_time_col]

    # Rule 4: drop cases with recurrence = 1 but missing raw recurrence time
    invalid_yes_no_time = (data[event_column] == 1) & (data[raw_rec_time_col].isna())
    if invalid_yes_no_time.any():
        logprint("Dropping samples with biochemical_recurrence = YES but missing recurrence time:")
        logprint(data.loc[invalid_yes_no_time, [patient_id_column, event_column, raw_rec_time_col, os_time_col]])
        data = data.loc[~invalid_yes_no_time].copy()

    # Drop rows where unified recurrence time is still missing
    data_clean = data.dropna(subset=[rec_time_col]).reset_index(drop=True)

    # Summary of OS.time censoring usage
    os_filled_mask = data_clean[raw_rec_time_col].isna()
    os_filled_summary = data_clean.loc[
        os_filled_mask,
        [patient_id_column, event_column, raw_rec_time_col, os_time_col, rec_time_col]
    ].copy()
    os_filled_summary["used_OS_time_as_censoring"] = True

    os_censor_summary_path = os.path.join(
        output_folder,
        "OS_Censoring_Summary_Black_Recurrence.csv"
    )
    os_filled_summary.to_csv(os_censor_summary_path, index=False)
    logprint(f"Summary of OS.time-based censoring saved to {os_censor_summary_path}")
    logprint(f"Number of samples censored at OS.time: {os_filled_summary.shape[0]}")

    # Gene columns: assume Ensembl IDs
    gene_cols = [c for c in data_clean.columns if c.startswith("ENSG")]

    n_total = data_clean.shape[0]
    logprint(f"Total usable patients after recurrence construction: {n_total}")
    logprint(f"Total gene features before selection: {len(gene_cols)}")
    logprint("Biochemical recurrence value counts (1 = YES, 0 = NO):")
    logprint(data_clean[event_column].value_counts())

    # KM curve for recurrence-free survival
    kmf = KaplanMeierFitter()
    kmf.fit(data_clean[rec_time_col], event_observed=data_clean[event_column])
    plt.figure(figsize=(8, 6))
    kmf.plot_survival_function()
    plt.title("Kaplan–Meier Curve – PRAD Black Cohort – Biochemical Recurrence-Free (OS.time censoring)")
    plt.xlabel("Days to recurrence or censoring")
    plt.ylabel("Recurrence-Free Probability")
    km_plot_path = os.path.join(output_folder, "KM_PRAD_Black_Recurrence_OSCensoring.png")
    plt.savefig(km_plot_path, dpi=300)
    plt.close()
    logprint(f"Kaplan–Meier recurrence curve saved to {km_plot_path}")

    # 10-fold CV
    k = 10
    kf = KFold(n_splits=k, shuffle=True, random_state=42)

    feature_selection_counts = {}
    per_patient_risk = []
    fold_idx = 1

    for train_idx, test_idx in kf.split(data_clean):
        logprint(f"\n=== Fold {fold_idx}/{k} ===")
        train_meta = data_clean.iloc[train_idx].reset_index(drop=True)
        test_meta  = data_clean.iloc[test_idx].reset_index(drop=True)

        y_tr_event = train_meta[event_column].values
        y_tr_time  = train_meta[rec_time_col].values
        X_tr_raw   = train_meta[gene_cols].values
        X_te_raw   = test_meta[gene_cols].values

        scaler = StandardScaler()
        X_tr = scaler.fit_transform(X_tr_raw)
        X_te = scaler.transform(X_te_raw)

        n_tr = X_tr.shape[0]
        inner_cv = min(5, max(2, n_tr // 2))
        logprint(f"  Training samples: {n_tr}, Inner CV folds: {inner_cv}")

        # Primary LASSO with timeout and ConvergenceWarning capture
        status, coef = fit_lasso_with_timeout(X_tr, y_tr_event, inner_cv, timeout_seconds=60)
        logprint(f"  Primary LASSO status: {status}")

        # Skip fold if ConvergenceWarning
        if status == "warn":
            logprint("  ConvergenceWarning detected. Skipping this fold.")
            fold_idx += 1
            continue

        # If timeout/fail or coef None, attempt fallback quick LASSO
        if (status in ["timeout", "fail"]) or (coef is None):
            logprint("  Attempting fallback LASSO (quick settings)...")
            coef = fallback_fit_lasso_quick(X_tr, y_tr_event, inner_cv)
            if coef is None:
                logprint("  Fallback LASSO failed to produce coefficients. Skipping this fold.")
                fold_idx += 1
                continue
            else:
                logprint("  Fallback LASSO produced coefficients.")

        selected_mask = coef != 0
        fold_features = np.array(gene_cols)[selected_mask]
        logprint(f"  Selected {len(fold_features)} features.")

        if len(fold_features) == 0:
            logprint("  No features selected after LASSO. Skipping Cox for this fold.")
            fold_idx += 1
            continue

        for f in fold_features:
            feature_selection_counts[f] = feature_selection_counts.get(f, 0) + 1

        X_tr_sel = X_tr[:, selected_mask]
        X_te_sel = X_te[:, selected_mask]

        train_df = pd.concat(
            [pd.DataFrame(X_tr_sel, columns=fold_features),
             pd.DataFrame({event_column: y_tr_event, rec_time_col: y_tr_time})],
            axis=1
        )
        test_df = pd.concat(
            [pd.DataFrame(X_te_sel, columns=fold_features),
             test_meta[[event_column, rec_time_col]].reset_index(drop=True)],
            axis=1
        )

        # Variance filter
        vars_ = train_df[fold_features].var()
        keep_features = vars_[vars_ > 1e-5].index.tolist()
        dropped = len(fold_features) - len(keep_features)
        logprint(f"  Dropped {dropped} near-zero variance features.")

        if len(keep_features) == 0:
            logprint("  After variance filter, no features remain. Skipping Cox for this fold.")
            fold_idx += 1
            continue

        train_df_cox = pd.concat(
            [train_df[keep_features], train_df[[event_column, rec_time_col]]],
            axis=1
        )
        test_df_cox = pd.concat(
            [test_df[keep_features], test_df[[event_column, rec_time_col]]],
            axis=1
        )

        cox = CoxPHFitter(penalizer=0.1)
        try:
            cox.fit(train_df_cox, duration_col=rec_time_col, event_col=event_column)
            logprint("  Cox model converged.")
            fold_risk = cox.predict_partial_hazard(test_df_cox[keep_features])
            per_patient_risk.append(pd.DataFrame({
                "Patient ID": test_meta[patient_id_column].values,
                "RiskScore_OutOfFold": fold_risk.values,
                "Fold": fold_idx
            }))
            logprint("  Out-of-fold risk scores computed.")
        except Exception as e:
            logprint(f"  Cox model failed to converge in this fold: {e}")

        fold_idx += 1

    # Aggregate results
    if per_patient_risk:
        risk_scores_all = pd.concat(per_patient_risk, axis=0, ignore_index=True)
    else:
        risk_scores_all = pd.DataFrame(columns=["Patient ID", "RiskScore_OutOfFold", "Fold"])

    merged_results = risk_scores_all.merge(
        data_clean[[patient_id_column, event_column, rec_time_col]],
        on="Patient ID",
        how="left"
    ).sort_values(by="Patient ID").reset_index(drop=True)

    feature_summary_df = (
        pd.DataFrame({
            "Feature": list(feature_selection_counts.keys()),
            "Times_Selected_Across_Folds": list(feature_selection_counts.values())
        }).sort_values(by="Times_Selected_Across_Folds", ascending=False)
        if feature_selection_counts else
        pd.DataFrame(columns=["Feature", "Times_Selected_Across_Folds"])
    )

    # Save outputs
    risk_scores_path = os.path.join(output_folder, "PRAD_Black_KFold10_RiskScores_Recurrence_OSCensoring.csv")
    merged_results.to_csv(risk_scores_path, index=False)
    logprint(f"Risk scores saved to {risk_scores_path}")

    feature_summary_path = os.path.join(output_folder, "PRAD_Black_KFold10_FeatureSelectionSummary_Recurrence_OSCensoring.csv")
    feature_summary_df.to_csv(feature_summary_path, index=False)
    logprint(f"Feature summary saved to {feature_summary_path}")

    all_features_path = os.path.join(output_folder, "PRAD_Black_KFold10_AllSelectedFeatures_Recurrence_OSCensoring.csv")
    feature_summary_df.to_csv(all_features_path, index=False)
    logprint(f"All unique features (with frequencies) saved to {all_features_path}")

    freq_threshold = 2
    stable_features_df = feature_summary_df[
        feature_summary_df["Times_Selected_Across_Folds"] >= freq_threshold
    ].reset_index(drop=True)
    stable_features_path = os.path.join(output_folder, "PRAD_Black_KFold10_StableFeatures_Recurrence_OSCensoring.csv")
    stable_features_df.to_csv(stable_features_path, index=False)
    logprint(f"Stable features saved to {stable_features_path}")

    # Summary logs
    usable_folds = risk_scores_all["Fold"].nunique() if not risk_scores_all.empty else 0
    log_lines = [
        "PRAD Black Cohort – 10-Fold LASSO + Cox (Biochemical Recurrence, OS.time Censoring)",
        "Skip policy: skip folds with ConvergenceWarning; timeouts/failures use quick fallback LASSO.",
        "-------------------------------------------------------------------------------------------",
        f"Total patients with usable recurrence info: {n_total}",
        f"Total genes pre-selection: {len(gene_cols)}",
        f"Usable folds contributing risk scores: {usable_folds}",
        f"Unique features selected across usable folds: {feature_summary_df.shape[0]}",
        f"Stability threshold: {freq_threshold} folds",
        f"Stable features: {stable_features_df.shape[0]}",
        "",
        "Recurrence time construction:",
        "* If days_to_first_biochemical_recurrence present, used as event time.",
        "* If biochemical_recurrence = 0 or missing, recurrence-free time censored at OS.time.",
        "* If biochemical_recurrence = 1 but recurrence time missing, sample dropped."
    ]
    summary_log_path = os.path.join(output_folder, "PRAD_Black_KFold10_Log_Recurrence_OSCensoring.txt")
    with open(summary_log_path, "w") as f:
        f.write("\n".join(log_lines))
    logprint(f"Summary log saved to {summary_log_path}")

    run_output_log_path = os.path.join(output_folder, "PRAD_Black_KFold10_RunOutputLog_Recurrence_OSCensoring.txt")
    with open(run_output_log_path, "w") as f:
        f.write("\n".join(run_log_lines))
    logprint(f"Full run output log saved to {run_output_log_path}")

    logprint("Done.")

if __name__ == "__main__":
    mp.set_start_method("spawn", force=True)
    main()
