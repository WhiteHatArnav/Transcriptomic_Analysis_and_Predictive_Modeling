import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LassoCV
from sklearn.model_selection import KFold
from lifelines import CoxPHFitter, KaplanMeierFitter

# -----------------------------
# Paths
# -----------------------------
input_folder = "/Users/arnavjoshi/Desktop/BRCAPredictiveAnalytics"
input_filename = "brca_all_survival_data.csv"
input_path = os.path.join(input_folder, input_filename)

output_folder = "/Users/arnavjoshi/Desktop/KFoldLassoCox/BRCAKFoldLassoCox"
os.makedirs(output_folder, exist_ok=True)

# -----------------------------
# Logging helper
# -----------------------------
run_log_lines = []

def logprint(msg):
    """Print to console and also record to run_log_lines."""
    print(msg)
    run_log_lines.append(str(msg))

# -----------------------------
# Load and prepare data
# -----------------------------
logprint(f"Loading {input_filename} ...")
data = pd.read_csv(input_path)

# Standardize column names
data.rename(columns={data.columns[0]: "Patient ID"}, inplace=True)
data.columns = data.columns.str.strip()

patient_id_column = "Patient ID"
event_column = "OS"
time_column = "OS.time"

# Drop rows missing survival info
data_clean = data.dropna(subset=[event_column, time_column]).reset_index(drop=True)

# Identify gene expression columns
gene_cols = [
    col for col in data_clean.columns
    if col not in [patient_id_column, event_column, time_column]
]

n_total_samples = data_clean.shape[0]
logprint(f"Total usable patients after survival cleanup: {n_total_samples}")
logprint(f"Total gene features before selection: {len(gene_cols)}")

logprint("OS value counts:")
logprint(data_clean[event_column].value_counts())

# -----------------------------
# Kaplan–Meier curve
# -----------------------------
kmf = KaplanMeierFitter()
T_all = data_clean[time_column]
E_all = data_clean[event_column]
kmf.fit(T_all, event_observed=E_all)

plt.figure(figsize=(8, 6))
ax = kmf.plot_survival_function()
ax.set_title("Kaplan–Meier Curve – BRCA All Cohort")
ax.set_xlabel("Time")
ax.set_ylabel("Survival Probability")

km_plot_path = os.path.join(output_folder, "KM_BRCA_All.png")
plt.savefig(km_plot_path, dpi=300)
plt.close()
logprint(f"Kaplan–Meier curve saved to {km_plot_path}")

# -----------------------------
# K-fold setup
# -----------------------------
k = 10
kf = KFold(n_splits=k, shuffle=True, random_state=42)

all_selected_features = set()
feature_selection_counts = {}
per_patient_risk_scores = []

fold_idx = 1

# -----------------------------
# K-fold loop
# -----------------------------
for train_index, test_index in kf.split(data_clean):
    logprint(f"\n=== Fold {fold_idx}/{k} ===")

    # Partition cohort
    train_meta = data_clean.iloc[train_index].reset_index(drop=True)
    test_meta = data_clean.iloc[test_index].reset_index(drop=True)

    y_train_event = train_meta[event_column].values
    y_train_time = train_meta[time_column].values

    X_train_raw = train_meta[gene_cols].values
    X_test_raw = test_meta[gene_cols].values

    # Standardize expression within this fold
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train_raw)
    X_test = scaler.transform(X_test_raw)

    # Inner CV folds for LASSO
    n_train_samples = X_train.shape[0]
    cv_folds_for_lasso = min(5, max(2, n_train_samples // 2))
    logprint(f"  Training samples: {n_train_samples}, Inner CV folds: {cv_folds_for_lasso}")

    # LASSO on event indicator (OS) for feature screening
    lasso_cv = LassoCV(
        alphas=np.logspace(-4, -1, 100),
        cv=cv_folds_for_lasso,
        random_state=42,
        max_iter=10000
    )

    logprint("  Fitting LASSO...")
    lasso_cv.fit(X_train, y_train_event)
    coef = lasso_cv.coef_

    selected_mask = coef != 0
    fold_selected_features = np.array(gene_cols)[selected_mask]

    logprint(f"  Fold {fold_idx}: selected {len(fold_selected_features)} features")

    # Track which genes were chosen this fold
    for f in fold_selected_features:
        all_selected_features.add(f)
        feature_selection_counts[f] = feature_selection_counts.get(f, 0) + 1

    # Build fold-specific Cox training frame
    X_train_sel = X_train[:, selected_mask]
    X_test_sel = X_test[:, selected_mask]

    train_df = pd.concat(
        [
            pd.DataFrame(X_train_sel, columns=fold_selected_features),
            pd.DataFrame({
                event_column: y_train_event,
                time_column: y_train_time
            })
        ],
        axis=1
    )

    test_df = pd.concat(
        [
            pd.DataFrame(X_test_sel, columns=fold_selected_features),
            test_meta[[event_column, time_column]].reset_index(drop=True)
        ],
        axis=1
    )

    # Drop near-zero variance features for Cox stability
    if len(fold_selected_features) > 0:
        variances = train_df[fold_selected_features].var()
        keep_features = variances[variances > 1e-5].index.tolist()
    else:
        keep_features = []

    dropped = len(fold_selected_features) - len(keep_features)
    logprint(f"  Dropped {dropped} near-zero variance features.")

    if len(keep_features) == 0:
        logprint("  No stable features remain for Cox; skipping Cox fit.")
        fold_idx += 1
        continue

    train_df_cox = pd.concat(
        [train_df[keep_features], train_df[[event_column, time_column]]],
        axis=1
    )
    test_df_cox = pd.concat(
        [test_df[keep_features], test_df[[event_column, time_column]]],
        axis=1
    )

    # Cox PH fit on the fold-specific features
    cox_model = CoxPHFitter(penalizer=0.1)
    try:
        cox_model.fit(
            train_df_cox,
            duration_col=time_column,
            event_col=event_column
        )
        logprint("  Cox model converged.")

        # Predict risk for held-out patients (out-of-fold risk)
        fold_risk = cox_model.predict_partial_hazard(test_df_cox[keep_features])

        per_patient_risk_scores.append(pd.DataFrame({
            "Patient ID": test_meta["Patient ID"].values,
            "RiskScore_OutOfFold": fold_risk.values,
            "Fold": fold_idx
        }))

        logprint("  Out-of-fold risk scores computed.")
    except Exception as e:
        logprint(f"  Cox model failed to converge: {e}")

    fold_idx += 1

# -----------------------------
# After K-fold: aggregate
# -----------------------------
if per_patient_risk_scores:
    risk_scores_all = pd.concat(per_patient_risk_scores, axis=0, ignore_index=True)
else:
    risk_scores_all = pd.DataFrame(columns=["Patient ID", "RiskScore_OutOfFold", "Fold"])

merged_results = risk_scores_all.merge(
    data_clean[["Patient ID", event_column, time_column]],
    on="Patient ID",
    how="left"
).sort_values(by="Patient ID").reset_index(drop=True)

# Feature frequency summary across folds
feature_summary_df = pd.DataFrame({
    "Feature": list(feature_selection_counts.keys()),
    "Times_Selected_Across_Folds": list(feature_selection_counts.values())
}).sort_values(by="Times_Selected_Across_Folds", ascending=False)

# Frequency threshold for "stable" features
# Threshold rule:
# - round(0.1 * n_total_samples), but at least 2 folds
freq_threshold = 2
logprint(f"\nFrequency threshold for stable features: {freq_threshold} folds")

stable_features_df = feature_summary_df[
    feature_summary_df["Times_Selected_Across_Folds"] >= freq_threshold
].reset_index(drop=True)

# -----------------------------
# Save outputs
# -----------------------------
risk_scores_path = os.path.join(output_folder, "BRCA_All_KFold10_RiskScores.csv")
merged_results.to_csv(risk_scores_path, index=False)
logprint(f"Risk scores saved to {risk_scores_path}")

feature_summary_path = os.path.join(output_folder, "BRCA_All_KFold10_FeatureSelectionSummary.csv")
feature_summary_df.to_csv(feature_summary_path, index=False)
logprint(f"Feature summary saved to {feature_summary_path}")

all_features_path = os.path.join(output_folder, "BRCA_All_KFold10_AllSelectedFeatures.csv")
feature_summary_df.to_csv(all_features_path, index=False)
logprint(f"All unique features (with frequencies) saved to {all_features_path}")

stable_features_path = os.path.join(output_folder, "BRCA_All_KFold10_StableFeatures.csv")
stable_features_df.to_csv(stable_features_path, index=False)
logprint(f"Stable features saved to {stable_features_path}")

# Summary log
log_lines = [
    "BRCA All Cohort – 10-Fold LASSO + Cox",
    "---------------------------------------------",
    f"Total patients: {n_total_samples}",
    f"Total genes pre-selection: {len(gene_cols)}",
    f"Unique features selected across folds: {len(feature_summary_df)}",
    f"Frequency threshold for stability: {freq_threshold}",
    f"Stable (high-frequency) features: {stable_features_df.shape[0]}",
    "",
    "Interpretation:",
    "* RiskScore_OutOfFold = Cox partial hazard for each patient from a model that did not train on that patient.",
    "* Times_Selected_Across_Folds reflects reproducibility of a gene being associated with event status.",
    "* High-frequency features are stronger survival-associated candidates.",
    "* Penalizer=0.1 is used in Cox for numerical stability.",
    "* 10-fold CV increases effective training size per model vs. a single 70/30 split, which is important in the All cohort because patient heterogeneity is high."
]
summary_log_path = os.path.join(output_folder, "BRCA_All_KFold10_Log.txt")
with open(summary_log_path, "w") as f:
    f.write("\n".join(log_lines))
logprint(f"Summary log saved to {summary_log_path}")

# Full run output log (console transcript)
run_output_log_path = os.path.join(output_folder, "BRCA_All_KFold10_RunOutputLog.txt")
with open(run_output_log_path, "w") as f:
    f.write("\n".join(run_log_lines))
logprint(f"Full run output log saved to {run_output_log_path}")

logprint("Done.")
