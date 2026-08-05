import os, time, warnings
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LassoCV, ElasticNetCV
from sklearn.exceptions import ConvergenceWarning
from lifelines import CoxPHFitter, KaplanMeierFitter
import matplotlib
matplotlib.use('Agg')  # no display
import matplotlib.pyplot as plt

# ===== Paths =====
input_folder = "/Users/arnavjoshi/Desktop/BRCAPredictiveAnalytics"
output_folder = input_folder
filename = "brca_black_survival_data.csv"
input_path = os.path.join(input_folder, filename)

print(f"Loading {filename}...")
data = pd.read_csv(input_path)
data.rename(columns={data.columns[0]: "Patient ID"}, inplace=True)
data.columns = data.columns.str.strip()

patient_id_column = "Patient ID"
event_column = "OS"
time_column = "OS.time"

# ===== Basic cleaning =====
data_clean = data.dropna(subset=[event_column, time_column])
gene_expression_columns = [
    c for c in data_clean.columns if c not in [patient_id_column, event_column, time_column]
]
print(f"Total genes (no prefilter): {len(gene_expression_columns)}")

# ===== Quick KM plot =====
kmf = KaplanMeierFitter()
T = data_clean[time_column]
E = data_clean[event_column]
kmf.fit(T, event_observed=E)
plt.figure(figsize=(8, 6))
ax = kmf.plot_survival_function()
ax.set_title("Kaplan-Meier Curve - Black")
ax.set_xlabel("Time")
ax.set_ylabel("Survival Probability")
plt.tight_layout()
plt.savefig(os.path.join(output_folder, "KM_Black.png"), dpi=300)
plt.close()
print("Kaplan-Meier curve saved.")

# ===== Standardize all genes (float32 for speed/memory) =====
expr = data_clean[gene_expression_columns]
scaler = StandardScaler(with_mean=True, with_std=True)
X = scaler.fit_transform(expr).astype(np.float32)
y = data_clean[[event_column, time_column]].reset_index(drop=True)

# ===== Split 70/30 =====
# If class imbalance is strong, stratify helps stability
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.30, random_state=42, stratify=y[event_column]
)

# ===== Primary model: LASSO (fast/stable settings with full p) =====
lasso_cv = LassoCV(
    alphas=np.logspace(-4, -2, 48),  # weaker/narrower to avoid over-shrinkage
    n_alphas=48,
    eps=1e-2,                        # coarser path => fewer steps
    cv=5,
    random_state=42,
    max_iter=200000,                 # give it room since p is huge
    tol=1e-3,
    selection="random",              # faster for very high p
    n_jobs=1                         # avoid BLAS oversubscription
)

print("Starting primary LASSO fit on ALL genes (no prefilter)...")
t0 = time.time()
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=ConvergenceWarning)
    lasso_cv.fit(X_train, y_train[event_column])
elapsed = time.time() - t0
print(f"LASSO fit done in {elapsed:.1f}s. Alpha*: {lasso_cv.alpha_:.3g}")

coef = lasso_cv.coef_.astype(np.float64, copy=True)
coef[np.abs(coef) < 1e-8] = 0.0  # remove numerical fuzz only
selected_mask = coef != 0
selected_features = np.array(gene_expression_columns)[selected_mask]
print(f"Selected {selected_mask.sum()} features with LASSO.")

# ===== Fallback: ElasticNet if LASSO is too sparse =====
used_model = "LASSO"
if selected_mask.sum() < 5:
    print("Few features from LASSO; trying ElasticNetCV fallback (l1_ratio≈0.9-1.0).")
    enet = ElasticNetCV(
        l1_ratio=[0.9, 0.95, 1.0],
        alphas=np.logspace(-4, -2, 48),
        cv=5,
        random_state=42,
        max_iter=200000,
        tol=1e-3,
        selection="random",
        n_jobs=1,
    )
    t1 = time.time()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        enet.fit(X_train, y_train[event_column])
    print(f"ElasticNet fit done in {time.time() - t1:.1f}s. "
          f"Alpha*: {enet.alpha_:.3g}, l1_ratio*: {enet.l1_ratio_}")

    coef = enet.coef_.astype(np.float64, copy=True)
    coef[np.abs(coef) < 1e-8] = 0.0
    selected_mask = coef != 0
    selected_features = np.array(gene_expression_columns)[selected_mask]
    print(f"Selected {selected_mask.sum()} features with ElasticNet.")
    used_model = "ElasticNet"

# ===== Save selected features =====
features_out = os.path.join(output_folder, "Features_Black_CV.csv")
pd.DataFrame({"Feature": selected_features, "Coefficient": coef[selected_mask]}).to_csv(
    features_out, index=False
)
print(f"Selected features saved to {features_out}")

# ===== Build train/test DataFrames for Cox =====
train_df = pd.concat(
    [pd.DataFrame(X_train[:, selected_mask], columns=selected_features),
     y_train.reset_index(drop=True)],
    axis=1
)
test_df = pd.concat(
    [pd.DataFrame(X_test[:, selected_mask], columns=selected_features),
     y_test.reset_index(drop=True)],
    axis=1
)

# Remove near-zero variance among selected features (final polish for Cox)
if len(selected_features) > 0:
    variances = train_df[selected_features].var()
    features_to_keep = variances[variances > 1e-5].index.tolist()
else:
    features_to_keep = []

print(f"Features kept for Cox: {len(features_to_keep)} (dropped "
      f"{len(selected_features) - len(features_to_keep)} near-zero var).")

train_df = pd.concat([train_df[features_to_keep], train_df[[event_column, time_column]]], axis=1)
test_df = pd.concat([test_df[features_to_keep], test_df[[event_column, time_column]]], axis=1)

# ===== Cox model =====
if len(features_to_keep) == 0:
    raise RuntimeError("No features survived selection+variance filter; consider widening alpha grid or ElasticNet.")

cox_model = CoxPHFitter(penalizer=0.1)
try:
    cox_model.fit(train_df, duration_col=time_column, event_col=event_column)
    print(f"Cox model fitted using {used_model}-selected features.")
    cox_model.print_summary()

    risk_scores = cox_model.predict_partial_hazard(test_df[features_to_keep])
    risk_df = pd.DataFrame({
        "Patient ID": data_clean.iloc[y_test.index][patient_id_column].values,
        "Risk Score": risk_scores.values
    })
    out_scores = os.path.join(output_folder, "RiskScores_Black_CV.csv")
    risk_df.to_csv(out_scores, index=False)
    print(f"Risk scores saved to {out_scores}")
except Exception as e:
    print("Cox model failed to converge or another error occurred:", e)
