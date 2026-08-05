import pandas as pd
import numpy as np
import os

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.linear_model import LassoCV

from lifelines import CoxPHFitter, KaplanMeierFitter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =====================================================
# PATHS
# =====================================================

data_path = (
    "/Users/arnavjoshi/Desktop/LUADPredictiveAnalyticsRecurrence/"
    "luad_black_DFI_survival_data.csv"
)

output_folder = (
    "/Users/arnavjoshi/Desktop/KFoldLassoCox/LUADKFoldLassoCoxRecurrence"
)
os.makedirs(output_folder, exist_ok=True)

# =====================================================
# LOAD DATA
# =====================================================

data = pd.read_csv(data_path)
data.columns = data.columns.str.strip()

patient_id_col = "Patient ID"
event_col = "DFI"
time_col = "DFI.time"

non_gene_cols = {patient_id_col, event_col, time_col}
gene_cols = [c for c in data.columns if c not in non_gene_cols]

X_all = data[gene_cols].values
y_event = data[event_col].values
y_time = data[time_col].values
patient_ids = data[patient_id_col].values

# =====================================================
# STORAGE
# =====================================================

oof_risk_scores = np.zeros(data.shape[0])
feature_selection_counts = {}

# =====================================================
# K-FOLD SETUP
# =====================================================

kf = KFold(n_splits=10, shuffle=True, random_state=42)

# =====================================================
# K-FOLD LOOP
# =====================================================

for fold, (train_idx, test_idx) in enumerate(kf.split(X_all), 1):
    print(f"\n===== Fold {fold} =====")

    X_train, X_test = X_all[train_idx], X_all[test_idx]
    y_train_event = y_event[train_idx]
    y_train_time = y_time[train_idx]

    # Standardization
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # -----------------------------
    # LASSO FEATURE SELECTION
    # -----------------------------
    lasso = LassoCV(
        alphas=np.logspace(-4, -1, 100),
        cv=5,
        max_iter=10000,
        random_state=42
    )
    lasso.fit(X_train_scaled, y_train_event)

    selected_mask = lasso.coef_ != 0
    selected_genes = np.array(gene_cols)[selected_mask]

    print(f"Selected {len(selected_genes)} features")

    if len(selected_genes) == 0:
        print("Skipping fold due to zero selected features")
        continue

    for gene in selected_genes:
        feature_selection_counts[gene] = (
            feature_selection_counts.get(gene, 0) + 1
        )

    # -----------------------------
    # COX MODEL
    # -----------------------------
    train_df = pd.DataFrame(
        X_train_scaled[:, selected_mask],
        columns=selected_genes
    )
    train_df[event_col] = y_train_event
    train_df[time_col] = y_train_time

    test_df = pd.DataFrame(
        X_test_scaled[:, selected_mask],
        columns=selected_genes
    )

    cox = CoxPHFitter(penalizer=0.1)
    cox.fit(train_df, duration_col=time_col, event_col=event_col)

    # -----------------------------
    # OUT-OF-FOLD RISK SCORES
    # -----------------------------
    oof_risk_scores[test_idx] = (
        cox.predict_partial_hazard(test_df).values.flatten()
    )

# =====================================================
# SAVE OOF RISK SCORES
# =====================================================

risk_df = pd.DataFrame({
    "Patient ID": patient_ids,
    "Risk Score": oof_risk_scores
})

risk_df.to_csv(
    os.path.join(
        output_folder,
        "RiskScores_LUAD_Black_DFI_KFold10.csv"
    ),
    index=False
)

# =====================================================
# FEATURE SELECTION SUMMARY
# =====================================================

feature_freq_df = (
    pd.DataFrame.from_dict(
        feature_selection_counts,
        orient="index",
        columns=["Selection_Count"]
    )
    .sort_values("Selection_Count", ascending=False)
    .reset_index()
    .rename(columns={"index": "Gene"})
)

feature_freq_df.to_csv(
    os.path.join(
        output_folder,
        "SelectedFeatures_LUAD_Black_DFI_KFold10_All.csv"
    ),
    index=False
)

# =====================================================
# STABLE FEATURE SET (≥ 2 folds)
# =====================================================

stable_features_df = feature_freq_df[
    feature_freq_df["Selection_Count"] >= 2
].reset_index(drop=True)

stable_features_df.to_csv(
    os.path.join(
        output_folder,
        "StableFeatures_LUAD_Black_DFI_KFold10_Min2.csv"
    ),
    index=False
)

print(
    f"\nStable features (≥2 folds): "
    f"{stable_features_df.shape[0]}"
)

# =====================================================
# KAPLAN–MEIER USING OOF RISK SCORES
# =====================================================

median_risk = risk_df["Risk Score"].median()
high_risk = risk_df["Risk Score"] >= median_risk

kmf = KaplanMeierFitter()
plt.figure(figsize=(8, 6))

kmf.fit(
    data.loc[high_risk, time_col],
    data.loc[high_risk, event_col],
    label="High Risk"
)
kmf.plot_survival_function()

kmf.fit(
    data.loc[~high_risk, time_col],
    data.loc[~high_risk, event_col],
    label="Low Risk"
)
kmf.plot_survival_function()

plt.title("LUAD Black – Disease-Free Interval (10-Fold LASSO–Cox)")
plt.xlabel("Days")
plt.ylabel("Disease-Free Probability")
plt.tight_layout()

plt.savefig(
    os.path.join(
        output_folder,
        "KM_LUAD_Black_DFI_KFold10.png"
    ),
    dpi=300
)
plt.close()

print("10-fold recurrence analysis with stable feature selection complete (LUAD Black).")
