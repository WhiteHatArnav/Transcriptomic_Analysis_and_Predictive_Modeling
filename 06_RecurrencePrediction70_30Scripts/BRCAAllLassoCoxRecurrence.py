import pandas as pd
import numpy as np
import os
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LassoCV
from lifelines import CoxPHFitter, KaplanMeierFitter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =====================================================
# PATHS
# =====================================================

# Original BRCA expression + OS survival file (ALL races)
expression_path = (
    "/Users/arnavjoshi/Desktop/BRCAPredictiveAnalytics/"
    "brca_all_survival_data.csv"
)

# BRCA DFI clinical file
dfi_path = "/Users/arnavjoshi/Desktop/BRCAOtherData/BRCA_DFI_Data.tsv"

# Output folder (recurrence analysis)
output_folder = "/Users/arnavjoshi/Desktop/BRCAPredictiveAnalyticsRecurrence"
os.makedirs(output_folder, exist_ok=True)

# New DFI survival file (created, not overwriting anything)
dfi_survival_path = os.path.join(
    output_folder,
    "brca_all_DFI_survival_data.csv"
)

# =====================================================
# PART 1 — CREATE DFI SURVIVAL FILE
# =====================================================

print("Loading original BRCA expression + OS file (ALL races)...")
expr = pd.read_csv(expression_path)
expr.rename(columns={expr.columns[0]: "Patient ID"}, inplace=True)
expr.columns = expr.columns.str.strip()

print("Loading BRCA DFI file...")
dfi = pd.read_csv(dfi_path, sep="\t")
dfi.rename(columns={dfi.columns[0]: "Patient ID"}, inplace=True)
dfi.columns = dfi.columns.str.strip()

# Keep only required DFI columns
dfi = dfi[["Patient ID", "DFI", "DFI.time"]]

# Merge
merged = expr.merge(dfi, on="Patient ID", how="left")

# Drop OS columns (non-destructive)
drop_cols = [c for c in ["OS", "OS.time"] if c in merged.columns]
merged = merged.drop(columns=drop_cols)

# Clean DFI fields
merged = merged.dropna(subset=["DFI", "DFI.time"])
merged["DFI.time"] = pd.to_numeric(merged["DFI.time"], errors="coerce")
merged = merged.dropna(subset=["DFI.time"])

# Save new survival file
merged.to_csv(dfi_survival_path, index=False)
print(f"DFI survival file saved:\n{dfi_survival_path}")

# =====================================================
# PART 2 — LASSO–COX RECURRENCE MODEL (70:30)
# =====================================================

print("Loading DFI survival file for modeling...")
data = pd.read_csv(dfi_survival_path)
data.columns = data.columns.str.strip()

patient_id_column = "Patient ID"
event_column = "DFI"
time_column = "DFI.time"

# Identify gene expression columns
non_gene_cols = {patient_id_column, event_column, time_column}
gene_expression_columns = [
    c for c in data.columns if c not in non_gene_cols
]

print(f"Total genes: {len(gene_expression_columns)}")
print("DFI event distribution (1 = recurrence, 0 = censored):")
print(data[event_column].value_counts())

# =====================================================
# Kaplan–Meier Curve (Disease-Free Interval)
# =====================================================

kmf = KaplanMeierFitter()
kmf.fit(
    durations=data[time_column],
    event_observed=data[event_column]
)

plt.figure(figsize=(8, 6))
ax = kmf.plot_survival_function()
ax.set_title("Kaplan–Meier Curve – BRCA All – Disease-Free Interval")
ax.set_xlabel("Days to first disease-free failure / censoring")
ax.set_ylabel("Disease-Free Probability")
plt.tight_layout()
plt.savefig(os.path.join(output_folder, "KM_BRCA_All_DFI.png"))
plt.close()

# =====================================================
# Standardization
# =====================================================

scaler = StandardScaler()
X = scaler.fit_transform(data[gene_expression_columns])
y = data[[event_column, time_column]]

# 70:30 split
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.3, random_state=42
)

# =====================================================
# LASSO Feature Selection (DFI event)
# =====================================================

lasso_cv = LassoCV(
    alphas=np.logspace(-4, -1, 100),
    cv=5,
    random_state=42,
    max_iter=10000
)

print("Running LASSO on DFI event...")
lasso_cv.fit(X_train, y_train[event_column])

selected_mask = lasso_cv.coef_ != 0
selected_features = np.array(gene_expression_columns)[selected_mask]
selected_coefficients = lasso_cv.coef_[selected_mask]

print(f"Selected {len(selected_features)} features.")

pd.DataFrame({
    "Feature": selected_features,
    "Coefficient": selected_coefficients
}).to_csv(
    os.path.join(output_folder, "Features_BRCA_All_DFI_CV.csv"),
    index=False
)

# =====================================================
# Cox Proportional Hazards Model
# =====================================================

train_df = pd.concat(
    [
        pd.DataFrame(X_train[:, selected_mask], columns=selected_features),
        y_train.reset_index(drop=True)
    ],
    axis=1
)

test_df = pd.concat(
    [
        pd.DataFrame(X_test[:, selected_mask], columns=selected_features),
        y_test.reset_index(drop=True)
    ],
    axis=1
)

# Remove near-zero variance genes
variances = train_df[selected_features].var()
features_to_keep = variances[variances > 1e-5].index.tolist()

train_df = train_df[features_to_keep + [event_column, time_column]]
test_df = test_df[features_to_keep + [event_column, time_column]]

cox_model = CoxPHFitter(penalizer=0.1)
cox_model.fit(train_df, duration_col=time_column, event_col=event_column)

# =====================================================
# Risk Scores
# =====================================================

risk_scores = cox_model.predict_partial_hazard(test_df[features_to_keep])
patient_ids_test = data.loc[y_test.index, patient_id_column].values

pd.DataFrame({
    "Patient ID": patient_ids_test,
    "Risk Score": risk_scores.values
}).to_csv(
    os.path.join(output_folder, "RiskScores_BRCA_All_DFI_CV.csv"),
    index=False
)

print("BRCA ALL-race DFI recurrence modeling complete.")
