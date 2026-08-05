import os, time, warnings
import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, LeaveOneOut
from sklearn.linear_model import LassoLarsCV, LassoLars, ElasticNetCV
from sklearn.exceptions import ConvergenceWarning

from lifelines import CoxPHFitter, KaplanMeierFitter
import matplotlib
matplotlib.use('Agg')  # no display
import matplotlib.pyplot as plt

# ===== Paths =====
input_folder = "/Users/arnavjoshi/Desktop/LIHCPredictiveAnalytics"
output_folder = input_folder
filename = "lihc_black_survival_data.csv"
input_path = os.path.join(input_folder, filename)

print(f"Loading {filename}...")
data = pd.read_csv(input_path)
data.rename(columns={data.columns[0]: "Patient ID"}, inplace=True)
data.columns = data.columns.str.strip()

patient_id_column = "Patient ID"
event_column = "OS"
time_column = "OS.time"

# ===== Basic cleaning =====
data_clean = data.dropna(subset=[event_column, time_column]).reset_index(drop=True)
gene_expression_columns = [c for c in data_clean.columns if c not in [patient_id_column, event_column, time_column]]
print(f"Total genes (no prefilter): {len(gene_expression_columns)}")

# ===== Quick KM plot =====
kmf = KaplanMeierFitter()
T = data_clean[time_column]
E = data_clean[event_column]
kmf.fit(T, event_observed=E)
plt.figure(figsize=(8, 6))
ax = kmf.plot_survival_function()
ax.set_title("Kaplan-Meier Curve - LIHC Black")
ax.set_xlabel("Time")
ax.set_ylabel("Survival Probability")
plt.tight_layout()
plt.savefig(os.path.join(output_folder, "KM_LIHC_Black.png"), dpi=300)
plt.close()
print("Kaplan-Meier curve saved.")

# ===== Standardize all genes (float32 for speed/memory) =====
expr = data_clean[gene_expression_columns]
scaler = StandardScaler(with_mean=True, with_std=True)
X = scaler.fit_transform(expr).astype(np.float32)
y = data_clean[[event_column, time_column]].reset_index(drop=True)

n_samples = X.shape[0]
print(f"Samples: {n_samples}")

# ===== Split 70/30 (no stratify for n very small) =====
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=min(0.30, 1.0 - 1.0/n_samples), random_state=42
)

# ===== CV scheme =====
cv_scheme = LeaveOneOut() if n_samples <= 10 else 5

# ===== Helper: pack selected features into DataFrames =====
def pack_frames(Xtr, Xte, cols, mask, ytr, yte):
    selected = np.array(cols)[mask]
    train_df = pd.concat(
        [pd.DataFrame(Xtr[:, mask], columns=selected), ytr.reset_index(drop=True)], axis=1
    )
    test_df = pd.concat(
        [pd.DataFrame(Xte[:, mask], columns=selected), yte.reset_index(drop=True)], axis=1
    )
    return selected, train_df, test_df

# ===== 1) Primary: LASSO via LARS path (fast for p≫n) =====
used_model = "LassoLarsCV"
print("Starting LassoLarsCV (LOO if n<=10)...")
t0 = time.time()
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=ConvergenceWarning)
    lars_cv = LassoLarsCV(cv=cv_scheme, max_iter=1000, n_jobs=1)
    lars_cv.fit(X_train, y_train[event_column].values)

elapsed = time.time() - t0
print(f"LassoLarsCV done in {elapsed:.1f}s.")
coef = lars_cv.coef_.astype(np.float64, copy=True)
coef[np.abs(coef) < 1e-8] = 0.0
mask = coef != 0
print(f"LassoLarsCV selected {mask.sum()} features.")

# ===== If too few features, try a slightly smaller alpha with plain LassoLars =====
TARGET_MIN = 10  # aim for at least ~10 features for Cox
if mask.sum() < TARGET_MIN:
    print("Few features from LassoLarsCV; nudging alpha smaller with LassoLars.")
    # Use a very small alpha to allow more features; warm-start is not needed for LARS
    # We try a small sequence of alphas until we exceed TARGET_MIN or hit floor.
    for alpha in [1e-6, 5e-7, 1e-7]:
        try:
            ll = LassoLars(alpha=alpha, max_iter=2000)
            ll.fit(X_train, y_train[event_column].values)
            coef2 = ll.coef_.astype(np.float64, copy=True)
            coef2[np.abs(coef2) < 1e-8] = 0.0
            mask2 = coef2 != 0
            print(f"LassoLars (alpha={alpha:g}) selected {mask2.sum()} features.")
            if mask2.sum() >= TARGET_MIN:
                coef, mask = coef2, mask2
                used_model = f"LassoLars(alpha={alpha:g})"
                break
        except Exception as e:
            print(f"LassoLars alpha={alpha:g} failed: {e}")

# ===== 2) Fallback: ElasticNetCV if still too sparse =====
if mask.sum() < TARGET_MIN:
    print("Still few features; trying ElasticNetCV fallback (tiny-n friendly).")
    enet = ElasticNetCV(
        l1_ratio=[0.8, 0.9, 0.95, 1.0],
        alphas=np.logspace(-5, -3, 36),
        cv=cv_scheme,
        random_state=42,
        max_iter=200000,
        tol=1e-3,
        selection="random",
        n_jobs=1,
    )
    t1 = time.time()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        enet.fit(X_train, y_train[event_column].values)
    print(f"ElasticNetCV done in {time.time() - t1:.1f}s. "
          f"alpha*={enet.alpha_:.3g}, l1_ratio*={enet.l1_ratio_}")
    coef = enet.coef_.astype(np.float64, copy=True)
    coef[np.abs(coef) < 1e-8] = 0.0
    mask = coef != 0
    used_model = "ElasticNetCV"

print(f"Final selected features: {mask.sum()} using {used_model}")

# ===== Save selected features =====
features = np.array(gene_expression_columns)[mask]
features_out = os.path.join(output_folder, "Features_LIHC_Black_CV.csv")
pd.DataFrame({"Feature": features, "Coefficient": coef[mask]}).to_csv(features_out, index=False)
print(f"Selected features saved to {features_out}")

# ===== Build train/test DataFrames for Cox =====
selected, train_df, test_df = pack_frames(
    X_train, X_test, gene_expression_columns, mask, y_train, y_test
)

# Final guard: if nothing survived, relax Cox penalizer will still allow a fit
if len(selected) == 0:
    raise RuntimeError("No features selected even after fallbacks; consider pooling cohorts or widening alpha range.")

# Remove near-zero variance among selected features (rare after standardization, but safe)
variances = train_df[selected].var()
keep = variances[variances > 1e-5].index.tolist()
dropped = len(selected) - len(keep)
print(f"Cox pre-clean: kept {len(keep)} features (dropped {dropped} near-zero var).")

train_df = pd.concat([train_df[keep], train_df[[event_column, time_column]]], axis=1)
test_df  = pd.concat([test_df[keep],  test_df[[event_column, time_column]]], axis=1)

# ===== Cox model (stronger penalizer for n tiny) =====
cox_model = CoxPHFitter(penalizer= 0.1 )
try:
    cox_model.fit(train_df, duration_col=time_column, event_col=event_column)
    print(f"Cox model fitted using {used_model}-selected features (penalizer=0.1).")
    cox_model.print_summary()

    risk_scores = cox_model.predict_partial_hazard(test_df[keep])
    risk_df = pd.DataFrame({
        "Patient ID": data_clean.iloc[y_test.index][patient_id_column].values,
        "Risk Score": risk_scores.values
    })
    out_scores = os.path.join(output_folder, "RiskScores_LIHC_Black_CV.csv")
    risk_df.to_csv(out_scores, index=False)
    print(f"Risk scores saved to {out_scores}")
except Exception as e:
    print("Cox model failed to converge or another error occurred:", e)
