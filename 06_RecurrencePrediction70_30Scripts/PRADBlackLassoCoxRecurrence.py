import pandas as pd
import numpy as np
import os
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LassoCV
from lifelines import CoxPHFitter, KaplanMeierFitter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# === Paths (PRAD BLACK) ===
expression_folder = "/Users/arnavjoshi/Desktop/PRADPredictiveAnalytics"
expression_filename = "prad_black_survival_data.csv"  # contains gene expression + Patient IDs (Black race)
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

# --- Recurrence variables ---
patient_id_column = "Patient ID"
event_column = "biochemical_recurrence"
time_column = "days_to_first_biochemical_recurrence"

# Normalize YES/NO to 1/0
rec[event_column] = rec[event_column].astype(str).str.strip().str.upper()
rec[event_column] = rec[event_column].map({"YES": 1, "NO": 0})

# Keep only relevant recurrence columns for merge
rec_subset = rec[[patient_id_column, event_column, time_column]]

# --- Merge expression + recurrence (on Patient ID) ---
data = expr.merge(rec_subset, on=patient_id_column, how="left")

# Drop rows missing recurrence info
data_clean = data.dropna(subset=[event_column, time_column])

# Ensure numeric time (use .loc to avoid SettingWithCopyWarning)
data_clean.loc[:, time_column] = pd.to_numeric(data_clean[time_column], errors="coerce")
data_clean = data_clean.dropna(subset=[time_column])

# Identify gene expression columns (everything except ID + recurrence vars)
non_gene_cols = {patient_id_column, event_column, time_column}
gene_expression_columns = [col for col in data_clean.columns if col not in non_gene_cols]

print(f"Total genes: {len(gene_expression_columns)}")

# Print recurrence counts
print("Biochemical recurrence value counts (1 = YES, 0 = NO):")
print(data_clean[event_column].value_counts(dropna=False))

# === Kaplan-Meier curve for recurrence-free survival ===
kmf = KaplanMeierFitter()
T = data_clean[time_column]
E = data_clean[event_column]

kmf.fit(T, event_observed=E)
plt.figure(figsize=(8, 6))
ax = kmf.plot_survival_function()
ax.set_title("Kaplan-Meier Curve - Black (PRAD) - Biochemical Recurrence-Free")
ax.set_xlabel("Days to first biochemical recurrence / censoring")
ax.set_ylabel("Recurrence-Free Probability")
plt.tight_layout()
plt.savefig(os.path.join(output_folder, "KM_Black_Recurrence.png"))
plt.close()
print("Kaplan-Meier recurrence curve saved.")

# === Standardize gene expression data ===
scaler = StandardScaler()
X = scaler.fit_transform(data_clean[gene_expression_columns])
y = data_clean[[event_column, time_column]]

# Train-test split (70% train, 30% test)
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.3, random_state=42
)

# === LassoCV for feature selection (on recurrence indicator) ===
lasso_cv = LassoCV(
    alphas=np.logspace(-4, -1, 100),
    cv=5,
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

# Remove near-zero variance features in training set to improve Cox convergence
variances = train_df[selected_features].var()
features_to_keep = variances[variances > 1e-5].index.tolist()

print(f"Removing {len(selected_features) - len(features_to_keep)} near-zero variance features before Cox fit.")

train_df = pd.concat(
    [train_df[features_to_keep], train_df[[event_column, time_column]]],
    axis=1
)
test_df = pd.concat(
    [test_df[features_to_keep], test_df[[event_column, time_column]]],
    axis=1
)

# === Fit Cox proportional hazards model with penalizer ===
cox_model = CoxPHFitter(penalizer=0.1)

try:
    cox_model.fit(train_df, duration_col=time_column, event_col=event_column)
    print("Cox Model Summary:")
    cox_model.print_summary()

    # Predict risk scores on test set
    risk_scores = cox_model.predict_partial_hazard(test_df[features_to_keep])

    # Use .loc for correct Patient ID alignment
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
