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

output_folder = "/Users/arnavjoshi/Desktop/KFoldLassoCox/PRADKFoldLassoCox"
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
# Fallback LASSO (no timeout): used ONLY if we timed out or failed
# NOTE: We IGNORE ConvergenceWarning here; we only skip on explicit warning
#       from the timed/primary run per your instruction.
# -----------------------------
def fallback_fit_lasso_quick(X_train, y_train_event, cv_folds_for_lasso):
    try:
        lasso_cv = LassoCV(
            alphas=np.logspace(-3, -1, 12),  # smaller grid, larger alphas
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
    data = pd.read_csv(input_path)
    data.rename(columns={data.columns[0]: "Patient ID"}, inplace=True)
    data.columns = data.columns.str.strip()

    patient_id_column = "Patient ID"
    event_column = "OS"
    time_column = "OS.time"

    data_clean = data.dropna(subset=[event_column, time_column]).reset_index(drop=True)
    gene_cols = [c for c in data_clean.columns if c not in [patient_id_column, event_column, time_column]]

    n_total = data_clean.shape[0]
    logprint(f"Total usable patients after survival cleanup: {n_total}")
    logprint(f"Total gene features before selection: {len(gene_cols)}")
    logprint("OS value counts:")
    logprint(data_clean[event_column].value_counts())

    # Descriptive KM
    kmf = KaplanMeierFitter()
    kmf.fit(data_clean[time_column], event_observed=data_clean[event_column])
    plt.figure(figsize=(8, 6))
    kmf.plot_survival_function()
    plt.title("Kaplan–Meier Curve – PRAD Black Cohort")
    plt.xlabel("Time")
    plt.ylabel("Survival Probability")
    km_plot_path = os.path.join(output_folder, "KM_PRAD_Black.png")
    plt.savefig(km_plot_path, dpi=300)
    plt.close()
    logprint(f"Kaplan–Meier curve saved to {km_plot_path}")

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
        y_tr_time  = train_meta[time_column].values
        X_tr_raw   = train_meta[gene_cols].values
        X_te_raw   = test_meta[gene_cols].values

        scaler = StandardScaler()
        X_tr = scaler.fit_transform(X_tr_raw)
        X_te = scaler.transform(X_te_raw)

        n_tr = X_tr.shape[0]
        inner_cv = min(5, max(2, n_tr // 2))
        logprint(f"  Training samples: {n_tr}, Inner CV folds: {inner_cv}")

        # Primary LASSO with timeout (captures ConvergenceWarning)
        status, coef = fit_lasso_with_timeout(X_tr, y_tr_event, inner_cv, timeout_seconds=60)
        logprint(f"  Primary LASSO status: {status}")

        # === SKIP CONDITION (as requested): skip ONLY if status == 'warn' ===
        if status == "warn":
            logprint("  ConvergenceWarning detected. Skipping this fold.")
            fold_idx += 1
            continue

        # If timeout/fail/None, try fallback (do NOT skip just because of timeout)
        if (status in ["timeout", "fail"]) or (coef is None):
            logprint("  Attempting fallback LASSO (quick settings)...")
            coef = fallback_fit_lasso_quick(X_tr, y_tr_event, inner_cv)
            if coef is None:
                logprint("  Fallback LASSO failed to produce coefficients. Skipping this fold.")
                fold_idx += 1
                continue
            else:
                logprint("  Fallback LASSO produced coefficients.")

        # Proceed with selected features
        selected_mask = coef != 0
        fold_features = np.array(gene_cols)[selected_mask]
        logprint(f"  Selected {len(fold_features)} features.")

        if len(fold_features) == 0:
            logprint("  No features selected after LASSO. Skipping Cox for this fold.")
            fold_idx += 1
            continue

        for f in fold_features:
            feature_selection_counts[f] = feature_selection_counts.get(f, 0) + 1

        # Prepare Cox data
        X_tr_sel = X_tr[:, selected_mask]
        X_te_sel = X_te[:, selected_mask]

        train_df = pd.concat(
            [pd.DataFrame(X_tr_sel, columns=fold_features),
             pd.DataFrame({event_column: y_tr_event, time_column: y_tr_time})],
            axis=1
        )
        test_df = pd.concat(
            [pd.DataFrame(X_te_sel, columns=fold_features),
             test_meta[[event_column, time_column]].reset_index(drop=True)],
            axis=1
        )

        # Variance filter before Cox
        vars_ = train_df[fold_features].var()
        keep_features = vars_[vars_ > 1e-5].index.tolist()
        dropped = len(fold_features) - len(keep_features)
        logprint(f"  Dropped {dropped} near-zero variance features.")

        if len(keep_features) == 0:
            logprint("  After variance filter, no features remain. Skipping Cox for this fold.")
            fold_idx += 1
            continue

        train_df_cox = pd.concat([train_df[keep_features], train_df[[event_column, time_column]]], axis=1)
        test_df_cox  = pd.concat([test_df[keep_features],  test_df[[event_column, time_column]]],  axis=1)

        cox = CoxPHFitter(penalizer=0.1)
        try:
            cox.fit(train_df_cox, duration_col=time_column, event_col=event_column)
            logprint("  Cox model converged.")
            fold_risk = cox.predict_partial_hazard(test_df_cox[keep_features])
            per_patient_risk.append(pd.DataFrame({
                "Patient ID": test_meta["Patient ID"].values,
                "RiskScore_OutOfFold": fold_risk.values,
                "Fold": fold_idx
            }))
            logprint("  Out-of-fold risk scores computed.")
        except Exception as e:
            logprint(f"  Cox model failed to converge in this fold: {e}")
            # continue to next fold

        fold_idx += 1

    # Aggregate results
    risk_scores_all = (
        pd.concat(per_patient_risk, axis=0, ignore_index=True)
        if per_patient_risk else
        pd.DataFrame(columns=["Patient ID", "RiskScore_OutOfFold", "Fold"])
    )

    merged_results = risk_scores_all.merge(
        data_clean[["Patient ID", event_column, time_column]],
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
    risk_scores_path = os.path.join(output_folder, "PRAD_Black_KFold10_RiskScores.csv")
    merged_results.to_csv(risk_scores_path, index=False)
    logprint(f"Risk scores saved to {risk_scores_path}")

    feature_summary_path = os.path.join(output_folder, "PRAD_Black_KFold10_FeatureSelectionSummary.csv")
    feature_summary_df.to_csv(feature_summary_path, index=False)
    logprint(f"Feature summary saved to {feature_summary_path}")

    all_features_path = os.path.join(output_folder, "PRAD_Black_KFold10_AllSelectedFeatures.csv")
    feature_summary_df.to_csv(all_features_path, index=False)
    logprint(f"All unique features (with frequencies) saved to {all_features_path}")

    # Choose a low threshold since some folds may be skipped
    freq_threshold = 2
    stable_features_df = feature_summary_df[
        feature_summary_df["Times_Selected_Across_Folds"] >= freq_threshold
    ].reset_index(drop=True)
    stable_features_path = os.path.join(output_folder, "PRAD_Black_KFold10_StableFeatures.csv")
    stable_features_df.to_csv(stable_features_path, index=False)
    logprint(f"Stable features saved to {stable_features_path}")

    # Summary logs
    log_lines = [
        "PRAD Black Cohort – 10-Fold LASSO + Cox",
        "Skip policy: ONLY skip if a ConvergenceWarning occurred; timeouts use a quick fallback LASSO.",
        "-------------------------------------------------------------------------",
        f"Total patients: {n_total}",
        f"Total genes pre-selection: {len(gene_cols)}",
        f"Usable folds contributing risk scores: {risk_scores_all['Fold'].nunique() if not risk_scores_all.empty else 0}",
        f"Unique features selected across usable folds: {feature_summary_df.shape[0]}",
        f"Stability threshold: {freq_threshold} folds",
        f"Stable features: {stable_features_df.shape[0]}"
    ]
    summary_log_path = os.path.join(output_folder, "PRAD_Black_KFold10_Log.txt")
    with open(summary_log_path, "w") as f:
        f.write("\n".join(log_lines))
    logprint(f"Summary log saved to {summary_log_path}")

    run_output_log_path = os.path.join(output_folder, "PRAD_Black_KFold10_RunOutputLog.txt")
    with open(run_output_log_path, "w") as f:
        f.write("\n".join(run_log_lines))
    logprint(f"Full run output log saved to {run_output_log_path}")

    logprint("Done.")

if __name__ == "__main__":
    mp.set_start_method("spawn", force=True)
    main()


############################### Code with Time Limit per Fold ################################
# import os
# import numpy as np
# import pandas as pd
# import warnings
# import multiprocessing as mp
# import matplotlib
# matplotlib.use("Agg")
# import matplotlib.pyplot as plt

# from sklearn.preprocessing import StandardScaler
# from sklearn.linear_model import LassoCV
# from sklearn.model_selection import KFold
# from sklearn.exceptions import ConvergenceWarning
# from lifelines import CoxPHFitter, KaplanMeierFitter

# # -----------------------------
# # Paths (PRAD, Black cohort)
# # -----------------------------
# input_folder = "/Users/arnavjoshi/Desktop/PRADPredictiveAnalytics"
# input_filename = "prad_black_survival_data.csv"
# input_path = os.path.join(input_folder, input_filename)

# output_folder = "/Users/arnavjoshi/Desktop/KFoldLassoCox/PRADKFoldLassoCox"
# os.makedirs(output_folder, exist_ok=True)

# # -----------------------------
# # Logging helper
# # -----------------------------
# run_log_lines = []

# def logprint(msg):
#     """Print to console immediately and also record in run_log_lines."""
#     print(msg, flush=True)
#     run_log_lines.append(str(msg))

# # -----------------------------
# # Helper: run LASSO in a subprocess with a time limit
# # -----------------------------
# def _fit_lasso_in_subprocess(X_train, y_train_event, cv_folds_for_lasso, return_dict):
#     """
#     Child process target. Fit LassoCV on (X_train, y_train_event).
#     Store result in return_dict:
#       return_dict["status"] in {"ok", "warn", "fail"}
#       return_dict["coef"]   = numpy array of coefficients or None
#     """
#     status = "ok"
#     coefs = None
#     try:
#         with warnings.catch_warnings(record=True) as wlist:
#             warnings.simplefilter("always", ConvergenceWarning)

#             lasso_cv = LassoCV(
#                 alphas=np.logspace(-4, -1, 20),   # reduced grid for speed
#                 cv=cv_folds_for_lasso,
#                 random_state=42,
#                 max_iter=2000,                    # lower max_iter to avoid grinding
#                 n_jobs=-1                         # parallelize across cores
#             )
#             lasso_cv.fit(X_train, y_train_event)
#             coefs = lasso_cv.coef_

#             # mark as "warn" if we saw ConvergenceWarning
#             for w in wlist:
#                 if issubclass(w.category, ConvergenceWarning):
#                     status = "warn"
#                     break

#     except Exception:
#         status = "fail"
#         coefs = None

#     return_dict["status"] = status
#     return_dict["coef"] = coefs

# def fit_lasso_with_timeout(X_train, y_train_event, cv_folds_for_lasso, timeout_seconds=60):
#     """
#     Launch _fit_lasso_in_subprocess(...) in a child process.
#     If it finishes within timeout_seconds: return (status, coef array or None)
#     If it does not finish: kill it and return ("timeout", None)
#     """
#     manager = mp.Manager()
#     shared_dict = manager.dict()

#     proc = mp.Process(
#         target=_fit_lasso_in_subprocess,
#         args=(X_train, y_train_event, cv_folds_for_lasso, shared_dict)
#     )
#     proc.start()
#     proc.join(timeout_seconds)

#     if proc.is_alive():
#         proc.terminate()
#         proc.join()
#         return ("timeout", None)

#     status = shared_dict.get("status", "fail")
#     coefs = shared_dict.get("coef", None)
#     return (status, coefs)

# # -----------------------------
# # Main
# # -----------------------------
# def main():
#     logprint(f"Loading {input_filename} ...")
#     data = pd.read_csv(input_path)

#     # Standardize column names
#     data.rename(columns={data.columns[0]: "Patient ID"}, inplace=True)
#     data.columns = data.columns.str.strip()

#     patient_id_column = "Patient ID"
#     event_column = "OS"
#     time_column = "OS.time"

#     # Drop rows missing survival info
#     data_clean = data.dropna(subset=[event_column, time_column]).reset_index(drop=True)

#     # Identify gene expression columns (all but ID/OS/OS.time)
#     gene_cols = [
#         col for col in data_clean.columns
#         if col not in [patient_id_column, event_column, time_column]
#     ]

#     n_total_samples = data_clean.shape[0]
#     logprint(f"Total usable patients after survival cleanup: {n_total_samples}")
#     logprint(f"Total gene features before selection: {len(gene_cols)}")

#     logprint("OS value counts:")
#     logprint(data_clean[event_column].value_counts())

#     # -----------------------------
#     # Kaplan–Meier curve (descriptive)
#     # -----------------------------
#     kmf = KaplanMeierFitter()
#     T_all = data_clean[time_column]
#     E_all = data_clean[event_column]
#     kmf.fit(T_all, event_observed=E_all)

#     plt.figure(figsize=(8, 6))
#     ax = kmf.plot_survival_function()
#     ax.set_title("Kaplan–Meier Curve – PRAD Black Cohort")
#     ax.set_xlabel("Time")
#     ax.set_ylabel("Survival Probability")

#     km_plot_path = os.path.join(output_folder, "KM_PRAD_Black.png")
#     plt.savefig(km_plot_path, dpi=300)
#     plt.close()
#     logprint(f"Kaplan–Meier curve saved to {km_plot_path}")

#     # -----------------------------
#     # K-fold setup (10-fold)
#     # -----------------------------
#     k = 10
#     kf = KFold(n_splits=k, shuffle=True, random_state=42)

#     all_selected_features = set()
#     feature_selection_counts = {}
#     per_patient_risk_scores = []

#     fold_idx = 1

#     # -----------------------------
#     # K-fold loop
#     # -----------------------------
#     for train_index, test_index in kf.split(data_clean):
#         logprint(f"\n=== Fold {fold_idx}/{k} ===")

#         # Split data for this fold
#         train_meta = data_clean.iloc[train_index].reset_index(drop=True)
#         test_meta = data_clean.iloc[test_index].reset_index(drop=True)

#         y_train_event = train_meta[event_column].values
#         y_train_time = train_meta[time_column].values

#         X_train_raw = train_meta[gene_cols].values
#         X_test_raw = test_meta[gene_cols].values

#         # Standardize features within this fold
#         scaler = StandardScaler()
#         X_train = scaler.fit_transform(X_train_raw)
#         X_test = scaler.transform(X_test_raw)

#         n_train_samples = X_train.shape[0]
#         cv_folds_for_lasso = min(5, max(2, n_train_samples // 2))
#         logprint(f"  Training samples: {n_train_samples}, Inner CV folds: {cv_folds_for_lasso}")
#         logprint("  Fitting LASSO with 60s timeout...")

#         # Time-constrained LASSO step
#         status, coef = fit_lasso_with_timeout(
#             X_train,
#             y_train_event,
#             cv_folds_for_lasso,
#             timeout_seconds=60
#         )

#         # status can be: "ok", "warn", "fail", or "timeout"
#         logprint(f"  LASSO status for this fold: {status}")

#         # Skip unusable folds
#         if status in ["timeout", "fail"] or coef is None:
#             logprint("  Skipping this fold (no usable LASSO result).")
#             fold_idx += 1
#             continue

#         # Treat ConvergenceWarning as unusable as requested
#         if status == "warn":
#             logprint("  LASSO produced convergence warning. Skipping this fold.")
#             fold_idx += 1
#             continue

#         # Non-zero coefficients mark selected features
#         selected_mask = coef != 0
#         fold_selected_features = np.array(gene_cols)[selected_mask]
#         logprint(f"  Fold {fold_idx}: selected {len(fold_selected_features)} features")

#         if len(fold_selected_features) == 0:
#             logprint("  No features selected after LASSO. Skipping Cox for this fold.")
#             fold_idx += 1
#             continue

#         # Update stability counts for features
#         for f in fold_selected_features:
#             all_selected_features.add(f)
#             feature_selection_counts[f] = feature_selection_counts.get(f, 0) + 1

#         # -----------------------------
#         # Prepare Cox PH data for this fold
#         # -----------------------------
#         X_train_sel = X_train[:, selected_mask]
#         X_test_sel = X_test[:, selected_mask]

#         train_df = pd.concat(
#             [
#                 pd.DataFrame(X_train_sel, columns=fold_selected_features),
#                 pd.DataFrame({
#                     event_column: y_train_event,
#                     time_column: y_train_time
#                 })
#             ],
#             axis=1
#         )

#         test_df = pd.concat(
#             [
#                 pd.DataFrame(X_test_sel, columns=fold_selected_features),
#                 test_meta[[event_column, time_column]].reset_index(drop=True)
#             ],
#             axis=1
#         )

#         # Drop near-zero variance features before Cox
#         variances = train_df[fold_selected_features].var()
#         keep_features = variances[variances > 1e-5].index.tolist()

#         dropped = len(fold_selected_features) - len(keep_features)
#         logprint(f"  Dropped {dropped} near-zero variance features.")

#         if len(keep_features) == 0:
#             logprint("  After variance filter, no features remain. Skipping Cox for this fold.")
#             fold_idx += 1
#             continue

#         train_df_cox = pd.concat(
#             [train_df[keep_features], train_df[[event_column, time_column]]],
#             axis=1
#         )
#         test_df_cox = pd.concat(
#             [test_df[keep_features], test_df[[event_column, time_column]]],
#             axis=1
#         )

#         # -----------------------------
#         # Cox PH model for this fold
#         # -----------------------------
#         cox_model = CoxPHFitter(penalizer=0.1)
#         try:
#             cox_model.fit(
#                 train_df_cox,
#                 duration_col=time_column,
#                 event_col=event_column
#             )
#             logprint("  Cox model converged.")

#             fold_risk = cox_model.predict_partial_hazard(test_df_cox[keep_features])

#             per_patient_risk_scores.append(pd.DataFrame({
#                 "Patient ID": test_meta["Patient ID"].values,
#                 "RiskScore_OutOfFold": fold_risk.values,
#                 "Fold": fold_idx
#             }))

#             logprint("  Out-of-fold risk scores computed.")
#         except Exception as e:
#             logprint(f"  Cox model failed to converge in this fold: {e}")
#             # Skip adding scores for this fold

#         fold_idx += 1

#     # -----------------------------
#     # After K-fold aggregation
#     # -----------------------------
#     if per_patient_risk_scores:
#         risk_scores_all = pd.concat(per_patient_risk_scores, axis=0, ignore_index=True)
#     else:
#         risk_scores_all = pd.DataFrame(columns=["Patient ID", "RiskScore_OutOfFold", "Fold"])

#     merged_results = risk_scores_all.merge(
#         data_clean[["Patient ID", event_column, time_column]],
#         on="Patient ID",
#         how="left"
#     ).sort_values(by="Patient ID").reset_index(drop=True)

#     # Build feature selection frequency table
#     feature_summary_df = pd.DataFrame({
#         "Feature": list(feature_selection_counts.keys()),
#         "Times_Selected_Across_Folds": list(feature_selection_counts.values())
#     }).sort_values(by="Times_Selected_Across_Folds", ascending=False)

#     # Stable features threshold: at least 2 folds (out of up to 10 usable folds)
#     freq_threshold = 2
#     logprint(f"\nFrequency threshold for stable features: {freq_threshold} folds")

#     stable_features_df = feature_summary_df[
#         feature_summary_df["Times_Selected_Across_Folds"] >= freq_threshold
#     ].reset_index(drop=True)

#     # -----------------------------
#     # Save outputs (PRAD Black, KFold10)
#     # -----------------------------
#     risk_scores_path = os.path.join(output_folder, "PRAD_Black_KFold10_RiskScores.csv")
#     merged_results.to_csv(risk_scores_path, index=False)
#     logprint(f"Risk scores saved to {risk_scores_path}")

#     feature_summary_path = os.path.join(output_folder, "PRAD_Black_KFold10_FeatureSelectionSummary.csv")
#     feature_summary_df.to_csv(feature_summary_path, index=False)
#     logprint(f"Feature summary saved to {feature_summary_path}")

#     all_features_path = os.path.join(output_folder, "PRAD_Black_KFold10_AllSelectedFeatures.csv")
#     feature_summary_df.to_csv(all_features_path, index=False)
#     logprint(f"All unique features (with frequencies) saved to {all_features_path}")

#     stable_features_path = os.path.join(output_folder, "PRAD_Black_KFold10_StableFeatures.csv")
#     stable_features_df.to_csv(stable_features_path, index=False)
#     logprint(f"Stable features saved to {stable_features_path}")

#     # Summary log
#     log_lines = [
#         "PRAD Black Cohort – 10-Fold LASSO + Cox (60s timeout per fold)",
#         "----------------------------------------------------------------",
#         f"Total patients: {n_total_samples}",
#         f"Total genes pre-selection: {len(gene_cols)}",
#         f"Unique features selected across usable folds: {len(feature_summary_df)}",
#         f"Frequency threshold for stability: {freq_threshold} folds",
#         f"Stable (high-frequency) features: {stable_features_df.shape[0]}",
#         "",
#         "Interpretation:",
#         "* Each fold attempts LASSO feature selection on 90% of patients, then fits Cox PH on those features.",
#         "* If LASSO timed out, failed, or raised a ConvergenceWarning, that fold was skipped (per your instruction).",
#         "* Remaining folds generated out-of-fold Cox risk scores (RiskScore_OutOfFold) for held-out patients.",
#         "* Stable features are those repeatedly selected across multiple usable folds.",
#     ]
#     summary_log_path = os.path.join(output_folder, "PRAD_Black_KFold10_Log.txt")
#     with open(summary_log_path, "w") as f:
#         f.write("\n".join(log_lines))
#     logprint(f"Summary log saved to {summary_log_path}")

#     # Full run log
#     run_output_log_path = os.path.join(output_folder, "PRAD_Black_KFold10_RunOutputLog.txt")
#     with open(run_output_log_path, "w") as f:
#         f.write("\n".join(run_log_lines))
#     logprint(f"Full run output log saved to {run_output_log_path}")

#     logprint("Done.")

# if __name__ == '__main__':
#     mp.set_start_method('spawn', force=True)  # safer on macOS/Windows
#     main()
