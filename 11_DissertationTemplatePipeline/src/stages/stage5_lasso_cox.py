from pathlib import Path
import numpy as np
import pandas as pd
from tqdm import tqdm

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, KFold
from sklearn.linear_model import LassoCV
from lifelines import CoxPHFitter, KaplanMeierFitter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pipeline.session import SessionParams
from pipeline.state import PipelineState


# --------------------------------------------------
# Survival matrix builder
# --------------------------------------------------
def _build_survival_matrix(
    expr_csv: Path,
    survival_tsv: Path,
    surv_sample_col: str,
    event_col: str,
    time_col: str,
    out_csv: Path
):
    surv_df = pd.read_csv(survival_tsv, sep="\t")
    surv_df[surv_sample_col] = surv_df[surv_sample_col].astype(str)

    expr_df = pd.read_csv(expr_csv, index_col=0, low_memory=False)

    # Row 0 is Sample ID
    expr_df = expr_df.iloc[1:, :]
    expr_df = expr_df.apply(pd.to_numeric, errors="coerce")

    expr_df = expr_df.T
    expr_df = np.log2(expr_df + 1)

    expr_df.index = expr_df.index.astype(str).str[:-1]
    expr_df.index.name = surv_sample_col

    merged_df = expr_df.merge(
        surv_df[[surv_sample_col, event_col, time_col]],
        on=surv_sample_col,
        how="left"
    )

    merged_df.to_csv(out_csv, index=False)


# --------------------------------------------------
# 70:30 LASSO–Cox
# --------------------------------------------------
def _run_lasso_cox_split(df, cfg, out_dir, cohort):
    pid = cfg["sample_id_column"]
    event = cfg["event_column"]
    time = cfg["time_column"]

    df = df.dropna(subset=[event, time])
    gene_cols = [c for c in df.columns if c not in [pid, event, time]]

    if len(df) < 5:
        print(f"    Split skipped: <5 samples after survival filtering ({len(df)})")
        return

    if cfg["save_overall_km"]:
        kmf = KaplanMeierFitter()
        kmf.fit(df[time], event_observed=df[event])
        plt.figure()
        kmf.plot_survival_function()
        plt.title(f"KM Curve – {cohort}")
        plt.savefig(out_dir / f"KM_{cohort}.png")
        plt.close()

    X = StandardScaler().fit_transform(df[gene_cols])
    y = df[[event, time]]

    X_tr, X_te, y_tr, y_te = train_test_split(
        X,
        y,
        test_size=cfg["split_test_size"],
        random_state=cfg["split_random_state"],
    )

    if cfg["lasso_cv_folds"] > X_tr.shape[0]:
        print(f"    Split skipped: CV folds > samples")
        return

    lasso = LassoCV(
        alphas=np.logspace(
            cfg["lasso_alphas_logspace_start"],
            cfg["lasso_alphas_logspace_stop"],
            cfg["lasso_alphas_n"],
        ),
        cv=cfg["lasso_cv_folds"],
        max_iter=cfg["lasso_max_iter"],
        random_state=cfg["lasso_random_state"],
    )

    try:
        lasso.fit(X_tr, y_tr[event])
    except Exception:
        print(f"    Split skipped: LASSO failed")
        return

    mask = lasso.coef_ != 0
    genes = np.array(gene_cols)[mask]
    coefs = lasso.coef_[mask]

    if len(genes) == 0:
        print(f"    No features selected for {cohort}")
        return

    pd.DataFrame(
        {"Feature": genes, "Coefficient": coefs}
    ).to_csv(out_dir / f"Features_{cohort}_CV.csv", index=False)

    train_df = pd.concat(
        [pd.DataFrame(X_tr[:, mask], columns=genes), y_tr.reset_index(drop=True)],
        axis=1,
    )

    cox = CoxPHFitter(penalizer=cfg["cox_penalizer"])
    try:
        cox.fit(train_df, duration_col=time, event_col=event)
    except Exception:
        print(f"    Cox failed for {cohort}")
        return

    risks = cox.predict_partial_hazard(
        pd.DataFrame(X_te[:, mask], columns=genes)
    )

    pd.DataFrame(
        {pid: df.loc[y_te.index, pid].values, "Risk Score": risks.values}
    ).to_csv(out_dir / f"RiskScores_{cohort}_CV.csv", index=False)


# --------------------------------------------------
# K-Fold LASSO–Cox (CORE)
# --------------------------------------------------
def _run_lasso_cox_kfold(df, cfg, out_dir, cohort):
    pid = cfg["sample_id_column"]
    event = cfg["event_column"]
    time = cfg["time_column"]

    df = df.dropna(subset=[event, time])
    gene_cols = [c for c in df.columns if c not in [pid, event, time]]

    n = len(df)
    k = cfg["kfold_splits"]

    if n < k:
        print(f"    KFold skipped: n ({n}) < k ({k})")
        return

    X = StandardScaler().fit_transform(df[gene_cols])
    y = df[[event, time]]

    kf = KFold(
        n_splits=k,
        shuffle=cfg["kfold_shuffle"],
        random_state=cfg["kfold_random_state"],
    )

    fold_selected = []
    risk_rows = []

    for i, (tr, te) in enumerate(tqdm(kf.split(X), total=k, desc=f"KFold {cohort}")):
        lasso = LassoCV(
            alphas=np.logspace(
                cfg["lasso_alphas_logspace_start"],
                cfg["lasso_alphas_logspace_stop"],
                cfg["lasso_alphas_n"],
            ),
            cv=cfg["lasso_cv_folds"],
            max_iter=cfg["lasso_max_iter"],
            random_state=cfg["lasso_random_state"],
        )

        try:
            lasso.fit(X[tr], y.iloc[tr][event])
        except Exception:
            continue

        mask = lasso.coef_ != 0
        if mask.sum() == 0:
            continue

        sel_genes = np.array(gene_cols)[mask]
        sel_coefs = lasso.coef_[mask]

        for g, c in zip(sel_genes, sel_coefs):
            fold_selected.append({"Fold": i, "Feature": g, "Coefficient": c})

        train_df = pd.concat(
            [
                pd.DataFrame(X[tr][:, mask], columns=sel_genes),
                y.iloc[tr].reset_index(drop=True),
            ],
            axis=1,
        )

        cox = CoxPHFitter(penalizer=cfg["cox_penalizer"])
        try:
            cox.fit(train_df, duration_col=time, event_col=event)
        except Exception:
            continue

        risks = cox.predict_partial_hazard(
            pd.DataFrame(X[te][:, mask], columns=sel_genes)
        )

        for sid, r in zip(df.iloc[te][pid].values, risks.values):
            risk_rows.append({"Sample ID": sid, "Risk Score": r})

    if fold_selected:
        pd.DataFrame(fold_selected).to_csv(
            out_dir / f"KFold_SelectedFeatures_{cohort}.csv", index=False
        )

        counts = (
            pd.DataFrame(fold_selected)
            .groupby("Feature")["Fold"]
            .nunique()
            .reset_index(name="Folds_Selected")
        )

        stable = counts[counts["Folds_Selected"] >= cfg["stable_feature_min_folds"]]
        stable.to_csv(out_dir / f"KFold_StableFeatures_{cohort}.csv", index=False)

    if risk_rows:
        pd.DataFrame(risk_rows).to_csv(
            out_dir / f"RiskScores_{cohort}_KFold.csv", index=False
        )


# --------------------------------------------------
# Stage 5 entry point
# --------------------------------------------------
def run_stage5(session: SessionParams, state: PipelineState):
    cfg = session.stage5
    out_dir = session.output_dir / "stage5"
    out_dir.mkdir(parents=True, exist_ok=True)

    for cohort in [
        session.race_1,
        session.race_2,
        session.all_races_label,
    ]:
        print(f"Stage 5 for {cohort}")

        expr_csv = state.get(f"expr_{cohort}")
        surv_csv = out_dir / f"{session.cancer_code}_{cohort}_survival_data.csv"

        _build_survival_matrix(
            expr_csv,
            Path(cfg["survival_tsv"]),
            cfg["sample_id_column"],
            cfg["event_column"],
            cfg["time_column"],
            surv_csv,
        )

        df = pd.read_csv(surv_csv)

        _run_lasso_cox_split(df, cfg, out_dir, cohort)
        _run_lasso_cox_kfold(df, cfg, out_dir, cohort)
