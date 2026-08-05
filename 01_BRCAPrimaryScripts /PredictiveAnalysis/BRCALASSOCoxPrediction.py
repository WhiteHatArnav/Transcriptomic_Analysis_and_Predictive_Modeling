import pandas as pd
import numpy as np
import os
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LassoCV
from lifelines import CoxPHFitter, KaplanMeierFitter
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt

# Input/output directories
input_folder = "/Users/arnavjoshi/Desktop/BRCAPredictiveAnalytics"
output_folder = "/Users/arnavjoshi/Desktop/BRCAPredictiveAnalytics"

# Input files and output base names
datasets = {
    "white": "brca_white_survival_data.csv",
    "black": "brca_black_survival_data.csv",
    "all": "brca_all_survival_data.csv"
}

for label, filename in datasets.items():
    print(f"\n=== Processing {label.upper()} dataset ===")
    
    input_path = os.path.join(input_folder, filename)
    data = pd.read_csv(input_path)

    # Rename first column to "Patient ID"
    data.rename(columns={data.columns[0]: "Patient ID"}, inplace=True)

    # Clean column names
    data.columns = data.columns.str.strip()

    # Identify columns
    patient_id_column = "Patient ID"
    event_column = "OS"
    time_column = "OS.time"
    gene_expression_columns = [
        col for col in data.columns if col not in [patient_id_column, event_column, time_column]
    ]

    # Drop rows with missing survival info
    data_clean = data.dropna(subset=[event_column, time_column])

    # Print OS value counts for sanity check
    print("OS value counts:")
    print(data_clean[event_column].value_counts())

    # Plot and save Kaplan-Meier curve without display
    kmf = KaplanMeierFitter()
    T = data_clean[time_column]
    E = data_clean[event_column]

    kmf.fit(T, event_observed=E)
    plt.figure(figsize=(8, 6))
    ax = kmf.plot_survival_function()
    ax.set_title(f"Kaplan-Meier Curve - {label.capitalize()}")
    ax.set_xlabel("Time")
    ax.set_ylabel("Survival Probability")
    km_plot_path = os.path.join(output_folder, f"KM_{label.capitalize()}.png")
    plt.savefig(km_plot_path)
    plt.close()
    print(f" Kaplan-Meier curve saved to {km_plot_path}")

    # Standardize gene expression
    scaler = StandardScaler()
    X = scaler.fit_transform(data_clean[gene_expression_columns])
    y = data_clean[[event_column, time_column]]

    # Train-test split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42
    )

    # LassoCV for feature selection
    lasso_cv = LassoCV(alphas=np.logspace(-4, -1, 100), cv=5, random_state=42)
    lasso_cv.fit(X_train, y_train[event_column])

    # Get selected features
    selected_mask = lasso_cv.coef_ != 0
    selected_features = np.array(gene_expression_columns)[selected_mask]
    selected_coefficients = lasso_cv.coef_[selected_mask]

    print(f"Selected {len(selected_features)} features.")

    # Save selected features
    selected_features_df = pd.DataFrame({
        "Feature": selected_features,
        "Coefficient": selected_coefficients
    })
    features_output_path = os.path.join(output_folder, f"Features_{label.capitalize()}_CV.csv")
    selected_features_df.to_csv(features_output_path, index=False)
    print(f" Selected features saved to {features_output_path}")

    # Build train/test dataframes with selected features
    train_df = pd.concat([
        pd.DataFrame(X_train[:, selected_mask], columns=selected_features),
        y_train.reset_index(drop=True)
    ], axis=1)

    test_df = pd.concat([
        pd.DataFrame(X_test[:, selected_mask], columns=selected_features),
        y_test.reset_index(drop=True)
    ], axis=1)

    # Fit Cox model
    cox_model = CoxPHFitter()
    cox_model.fit(train_df, duration_col=time_column, event_col=event_column)
    print("Cox Model Summary:")
    cox_model.print_summary()

    # Predict risk scores on test set
    risk_scores = cox_model.predict_partial_hazard(test_df[selected_features])

    # Save risk scores
    risk_scores_df = pd.DataFrame({
        "Patient ID": data_clean.loc[y_test.index, "Patient ID"].values,
        "Risk Score": risk_scores.values
    })
    risk_scores_output_path = os.path.join(output_folder, f"RiskScores_{label.capitalize()}_CV.csv")
    risk_scores_df.to_csv(risk_scores_output_path, index=False)
    print(f" Risk scores saved to {risk_scores_output_path}")
