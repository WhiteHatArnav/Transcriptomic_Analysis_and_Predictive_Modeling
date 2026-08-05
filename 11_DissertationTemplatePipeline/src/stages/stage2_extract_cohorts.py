
import pandas as pd
from pipeline.session import SessionParams
from pipeline.state import PipelineState

def run_stage2(session: SessionParams, state: PipelineState):
    outdir = session.output_dir / "stage2"
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(state.get("expression_matrix"), index_col=0, low_memory=False)

    sid = df.loc["Sample ID"]
    race = df.loc["Race"]
    stype = df.loc["Sample Type"]
    subtype = df.loc["Sample Subtype"]

    cohorts = {
        session.race_1: session.race_1,
        session.race_2: session.race_2,
        session.all_races_label: None
    }

    for label, rfilter in cohorts.items():
        tumor_mask = subtype.isin(session.aggressive_subtypes) & stype.isin(session.tumor_sample_types)
        control_mask = stype.isin(session.control_sample_types)

        if rfilter:
            race_mask = race.str.lower() == rfilter.lower()
            tumor_mask &= race_mask
            control_mask &= race_mask

        cols = df.columns[tumor_mask | control_mask]
        tumor_cols = set(df.columns[tumor_mask]).intersection(set(cols))

        expr = df[cols].drop(index=["Race","Sample Type","Sample Subtype"])
        expr.columns = sid[cols].values

        meta = pd.DataFrame({
            "Sample ID": sid[cols].values,
            "Condition": ["Tumor" if c in tumor_cols else "Control" for c in cols]
        })

        expr_out = outdir / f"{session.cancer_code}_{label}_expression_aggressive.csv"
        meta_out = outdir / f"{session.cancer_code}_{label}_metadata.csv"

        expr.to_csv(expr_out)
        meta.to_csv(meta_out, index=False)

        state.register(f"expr_{label}", expr_out)
        state.register(f"meta_{label}", meta_out)
