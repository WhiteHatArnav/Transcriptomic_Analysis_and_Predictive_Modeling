
from pipeline.session import collect_session_params
from pipeline.state import PipelineState

from stages.stage1_expression_matrix import run_stage1
from stages.stage2_extract_cohorts import run_stage2
from stages.stage3_deseq2 import run_stage3
from stages.stage4_go_enrichment import run_stage4
from stages.stage5_lasso_cox import run_stage5

def main():
    session = collect_session_params()
    state = PipelineState()

    if session.run_stages.get("stage1", True):
        print("Running Stage 1")
        state.register("expression_matrix", run_stage1(session))

    if session.run_stages.get("stage2", True):
        print("Running Stage 2")
        run_stage2(session, state)

    if session.run_stages.get("stage3", True):
        print("Running Stage 3")
        run_stage3(session, state)

    if session.run_stages.get("stage4", False):
        print("Running Stage 4")
        run_stage4(session, state)

    if session.run_stages.get("stage5", False):
        print("Running Stage 5")
        run_stage5(session, state)

    print("Pipeline complete")

if __name__ == "__main__":
    main()
