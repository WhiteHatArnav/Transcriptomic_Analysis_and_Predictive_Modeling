
from pathlib import Path
from pipeline.session import SessionParams
from pipeline.state import PipelineState
from pipeline.runner_r import run_rscript

def run_stage3(session: SessionParams, state: PipelineState):
    outdir = session.output_dir / "stage3"
    outdir.mkdir(parents=True, exist_ok=True)

    base = Path(__file__).resolve().parents[2]

    for cohort in [session.race_1, session.race_2, session.all_races_label]:
        deseq_out = outdir / f"{session.cancer_code}_{cohort}_DESeq2.csv"
        volcano_out = outdir / f"{session.cancer_code}_{cohort}_volcano.pdf"

        run_rscript(
            base / "R" / "deseq2_analysis.R",
            [state.get(f"expr_{cohort}"), state.get(f"meta_{cohort}"), deseq_out]
        )

        state.register(f"deseq_{cohort}", deseq_out)

        run_rscript(
            base / "R" / "volcano_plot.R",
            [deseq_out, volcano_out, f"{cohort} {session.cancer_code} Volcano Plot"]
        )
