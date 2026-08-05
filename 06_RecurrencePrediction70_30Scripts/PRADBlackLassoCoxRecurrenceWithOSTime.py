import pandas as pd
import numpy as np
import os
import sys
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LassoCV
from lifelines import CoxPHFitter, KaplanMeierFitter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# === Paths (PRAD BLACK) ===
expression_folder = "/Users/arnavjoshi/Desktop/PRADPredictiveAnalytics"
expression_filename = "prad_black_survival_data.csv"  # Black cohort expression + OS
expression_path = os.path.join(expression_folder, expression_filename)

recurrence_path = "/Users/arnavjoshi/Desktop/PRADOtherData/PRADReccurenceData.tsv"

output_folder = "/Users/arnavjoshi/Desktop/PRADPredictiveAnalyticsRecurrence"
os.makedirs(output_folder, exist_ok=True)

print(f"Loading expression file: {expression_filename}...")
expr = pd.read_csv(expression_path)
expr.rename(columns={expr.columns[0]: "Patient ID"}, inplace=True)
expr.columns = expr.columns.str.strip()

print("Loading recurrence file...")
rec = pd.read_csv(recurrence_path, sep="\t")
rec.rename(columns={rec.columns[0]: "Patient ID"}, inplace=True)
rec.columns = rec.columns.str.strip()

# --- Column names ---
patient_id_column = "Patient ID"
event_column = "biochemical_recurrence"                     # YES / NO / NaN
raw_rec_time_col = "days_to_first_biochemical_recurrence"  # event time when present
os_time_col = "OS.time"                                     # from survival CSV
rec_time_col = "Recurrence_time_days"                       # unified time column

# === Prepare recurrence data ===

# Map YES/NO/NaN -> 1/0 (NaN treated as 0 = no documented recurrence)
rec[event_column] = rec[event_column].astype(str).str.strip().str.upper()
rec[event_column] = rec[event_column].map({"YES": 1, "NO": 0})
# After mapping, remaining NaN = no documented recurrence
rec[event_column] = rec[event_column].fillna(0).astype(int)

# Keep only relevant recurrence columns
rec_subset = rec[[patient_id_column, event_column, raw_rec_time_col]]

# === Merge expression + recurrence on Patient ID ===
data = expr.merge(rec_subset, on=patient_id_column, how="left")

# Ensure OS.time is present and numeric
if os_time_col not in data.columns:
    print(f"ERROR: {os_time_col} not found in {expression_filename}. "
          "Cannot construct recurrence time.")
    sys.exit(1)

data.loc[:, os_time_col] = pd.to_numeric(data[os_time_col], errors="coerce")

# Ensure raw recurrence time is numeric
data.loc[:, raw_rec_time_col] = pd.to_numeric(data[raw_rec_time_col], errors="coerce")

# For any rows where event_column is NaN after merge (should be rare), treat as 0
data.loc[:, event_column] = data[event_column].fillna(0).astype(int)

# Build unified recurrence time column
# Step 1: start with raw recurrence time
data[rec_time_col] = data[raw_rec_time_col]

# Step 2: for rows where recurrence time is missing, fill from OS.time
mask_missing_rec_time = data[rec_time_col].isna()
data.loc[mask_missing_rec_time, rec_time_col] = data.loc[mask_missing_rec_time, os_time_col]

# Step 3: identify invalid rows:
# event = 1 (YES) but raw recurrence time is missing -> cannot assign OS.time as event time
invalid_yes_no_time = (data[event_column] == 1) & (data[raw_rec_time_col].isna())
if invalid_yes_no_time.any():
    print("Dropping samples with biochemical_recurrence = YES but missing recurrence time:")
    print(data.loc[invalid_yes_no_time, [patient_id_column, event_column, raw_rec_time_col, os_time_col]])
    data = data.loc[~invalid_yes_no_time].copy()

# Step 4: drop rows where unified recurrence time or OS.time still missing
data_clean = data.dropna(subset=[rec_time_col, os_time_col])

print(f"Total samples in Black expression file: {expr.shape[0]}")
print(f"Samples with usable recurrence info after processing: {data_clean.shape[0]}")

# === Summary CSV of samples where OS.time was used as recurrence time ===
# OS.time is used whenever raw_rec_time_col was NaN but rec_time_col was filled from OS.time,
# and event_column == 0 (no documented recurrence).
os_filled_mask = data_clean[raw_rec_time_col].isna()
os_filled_summary = data_clean.loc[os_filled_mask, [
    patient_id_column,
    event_column,
    raw_rec_time_col,
    os_time_col,
    rec_time_col
]].copy()

# Add a simple flag for clarity
os_filled_summary["used_OS_time_as_censoring"] = True

summary_output_path = os.path.join(
    output_folder,
    "OS_Censoring_Summary_Black_Recurrence.csv"
)
os_filled_summary.to_csv(summary_output_path, index=False)
print(f"Summary of OS.time-based censoring saved to {summary_output_path}")
print(f"Number of samples censored at OS.time: {os_filled_summary.shape[0]}")

if data_clean.shape[0] < 3:
    print("Warning: very few samples with usable recurrence information "
          f"(n = {data_clean.shape[0]}). Lasso-Cox results may be unstable.")

# === Identify gene expression columns ===
non_gene_cols = {
    patient_id_column,
    event_column,
    rec_time_col,
    raw_rec_time_col,
    os_time_col,
    "OS"  # if OS column exists
}
gene_expression_columns = [col for col in data_clean.columns if col not in non_gene_cols]

print(f"Total genes: {len(gene_expression_columns)}")

# Print recurrence counts
print("Biochemical recurrence (1 = YES, 0 = NO) value counts:")
print(data_clean[event_column].value_counts(dropna=False))

# === Kaplan-Meier curve for recurrence-free survival ===
kmf = KaplanMeierFitter()
T = data_clean[rec_time_col]
E = data_clean[event_column]

kmf.fit(T, event_observed=E)
plt.figure(figsize=(8, 6))
ax = kmf.plot_survival_function()
ax.set_title("Kaplan-Meier Curve - Black (PRAD) - Biochemical Recurrence-Free")
ax.set_xlabel("Days to recurrence or censoring")
ax.set_ylabel("Recurrence-Free Probability")
plt.tight_layout()
km_path = os.path.join(output_folder, "KM_Black_Recurrence.png")
plt.savefig(km_path)
plt.close()
print(f"Kaplan-Meier recurrence curve saved to {km_path}")

# === Standardize gene expression data ===
scaler = StandardScaler()
X = scaler.fit_transform(data_clean[gene_expression_columns])
y = data_clean[[event_column, rec_time_col]]

# Train-test split (70% train, 30% test)
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.3, random_state=42
)

n_train = X_train.shape[0]
print(f"Number of training samples after split: {n_train}")

if n_train < 3:
    print("Not enough training samples for reliable LassoCV and CoxPH in the Black recurrence cohort "
          f"(n_train = {n_train}). Skipping Lasso-Cox modeling.")
    sys.exit(0)

# === LassoCV for feature selection (on recurrence indicator) ===
cv_folds = min(5, n_train)
print(f"Using cv={cv_folds} folds for LassoCV.")

lasso_cv = LassoCV(
    alphas=np.logspace(-4, -1, 100),
    cv=cv_folds,
    random_state=42,
    max_iter=10000
)

print("Starting LassoCV fitting on biochemical recurrence...")
lasso_cv.fit(X_train, y_train[event_column])
print("LassoCV fitting done.")

selected_mask = lasso_cv.coef_ != 0
selected_features = np.array(gene_expression_columns)[selected_mask]
selected_coefficients = lasso_cv.coef_[selected_mask]

print(f"Selected {len(selected_features)} features.")

selected_features_df = pd.DataFrame({
    "Feature": selected_features,
    "Coefficient": selected_coefficients
})

features_output_path = os.path.join(output_folder, "Features_Black_Recurrence_CV.csv")
selected_features_df.to_csv(features_output_path, index=False)
print(f"Selected features saved to {features_output_path}")

# === Prepare train/test data for Cox regression ===
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

# Remove near-zero variance features
if len(selected_features) > 0:
    variances = train_df[selected_features].var()
    features_to_keep = variances[variances > 1e-5].index.tolist()
else:
    features_to_keep = []

print(f"Removing {len(selected_features) - len(features_to_keep)} near-zero variance features before Cox fit.")

train_df = pd.concat(
    [train_df[features_to_keep], train_df[[event_column, rec_time_col]]],
    axis=1
)
test_df = pd.concat(
    [test_df[features_to_keep], test_df[[event_column, rec_time_col]]],
    axis=1
)

# === Cox proportional hazards model ===
cox_model = CoxPHFitter(penalizer=0.1)

try:
    if len(features_to_keep) == 0:
        print("No features were selected by Lasso for the Black recurrence cohort. "
              "Skipping CoxPH fitting.")
        sys.exit(0)

    cox_model.fit(train_df, duration_col=rec_time_col, event_col=event_column)
    print("Cox Model Summary:")
    cox_model.print_summary()

    # Predict risk scores on test set
    risk_scores = cox_model.predict_partial_hazard(test_df[features_to_keep])

    # Align Patient IDs with test indices
    patient_ids_test = data_clean.loc[y_test.index, patient_id_column].values

    risk_scores_df = pd.DataFrame({
        "Patient ID": patient_ids_test,
        "Risk Score": risk_scores.values
    })

    risk_scores_output_path = os.path.join(output_folder, "RiskScores_Black_Recurrence_CV.csv")
    risk_scores_df.to_csv(risk_scores_output_path, index=False)
    print(f"Risk scores saved to {risk_scores_output_path}")

except Exception as e:
    print("Cox model failed to converge or another error occurred:", e)
