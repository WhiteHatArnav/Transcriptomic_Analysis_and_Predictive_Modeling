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

# === Paths ===
input_folder = "/Users/arnavjoshi/Desktop/LIHCPredictiveAnalytics"
output_folder = input_folder

filename = "lihc_white_survival_data.csv"
input_path = os.path.join(input_folder, filename)

print(f"Loading {filename}...")
data = pd.read_csv(input_path)
data.rename(columns={data.columns[0]: "Patient ID"}, inplace=True)
data.columns = data.columns.str.strip()

patient_id_column = "Patient ID"
event_column = "OS"
time_column = "OS.time"

# Drop rows missing survival info
data_clean = data.dropna(subset=[event_column, time_column])

# Identify gene expression columns
gene_expression_columns = [
    col for col in data_clean.columns if col not in [patient_id_column, event_column, time_column]
]

print(f"Total genes: {len(gene_expression_columns)}")

# Print OS counts
print("OS value counts:")
print(data_clean[event_column].value_counts())

# Plot and save Kaplan-Meier curve
kmf = KaplanMeierFitter()
T = data_clean[time_column]
E = data_clean[event_column]
kmf.fit(T, event_observed=E)
plt.figure(figsize=(8, 6))
ax = kmf.plot_survival_function()
ax.set_title("Kaplan-Meier Curve - White")
ax.set_xlabel("Time")
ax.set_ylabel("Survival Probability")
plt.savefig(os.path.join(output_folder, "KM_White.png"))
plt.close()
print("Kaplan-Meier curve saved.")

# Standardize gene expression data
scaler = StandardScaler()
X = scaler.fit_transform(data_clean[gene_expression_columns])
y = data_clean[[event_column, time_column]]

# Train-test split (70% train, 30% test)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# LassoCV for feature selection
lasso_cv = LassoCV(
    alphas=np.logspace(-4, -1, 100),
    cv=5,
    random_state=42,
    max_iter=10000
)
print("Starting LassoCV fitting...")
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
features_output_path = os.path.join(output_folder, "Features_White_CV.csv")
selected_features_df.to_csv(features_output_path, index=False)
print(f"Selected features saved to {features_output_path}")

# Prepare train/test data for Cox regression
train_df = pd.concat([
    pd.DataFrame(X_train[:, selected_mask], columns=selected_features),
    y_train.reset_index(drop=True)
], axis=1)

test_df = pd.concat([
    pd.DataFrame(X_test[:, selected_mask], columns=selected_features),
    y_test.reset_index(drop=True)
], axis=1)

# Remove near-zero variance features in training set to improve Cox convergence
variances = train_df[selected_features].var()
features_to_keep = variances[variances > 1e-5].index.tolist()

print(f"Removing {len(selected_features) - len(features_to_keep)} near-zero variance features before Cox fit.")

train_df = pd.concat([train_df[features_to_keep], train_df[[event_column, time_column]]], axis=1)
test_df = pd.concat([test_df[features_to_keep], test_df[[event_column, time_column]]], axis=1)

# Fit Cox proportional hazards model with penalizer
cox_model = CoxPHFitter(penalizer=0.1)
try:
    cox_model.fit(train_df, duration_col=time_column, event_col=event_column)
    print("Cox Model Summary:")
    cox_model.print_summary()

    # Predict risk scores on test set
    risk_scores = cox_model.predict_partial_hazard(test_df[features_to_keep])

    risk_scores_df = pd.DataFrame({
        "Patient ID": data_clean.loc[y_test.index, "Patient ID"].values,
        "Risk Score": risk_scores.values
    })
    risk_scores_output_path = os.path.join(output_folder, "RiskScores_White_CV.csv")
    risk_scores_df.to_csv(risk_scores_output_path, index=False)
    print(f"Risk scores saved to {risk_scores_output_path}")
except Exception as e:
    print("Cox model failed to converge or another error occurred:", e)
