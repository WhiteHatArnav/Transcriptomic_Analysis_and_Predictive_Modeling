
from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Any
import yaml, argparse

@dataclass
class SessionParams:
    cancer_code: str
    input_dir: Path
    output_dir: Path
    sample_sheet_tsv: Path
    clinical_tsv: Path
    subtype_tsv: Path
    subtype_sample_id_col: str
    subtype_value_col: str
    aggressive_subtypes: List[str]
    tumor_sample_types: List[str]
    control_sample_types: List[str]
    race_1: str
    race_2: str
    all_races_label: str
    run_stages: Dict[str, bool]
    stage4: Dict[str, Any]
    stage5: Dict[str, Any]

def collect_session_params():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=str, required=True)
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    def g(path, default=None):
        cur = cfg
        for p in path:
            if cur is None or p not in cur:
                return default
            cur = cur[p]
        return cur

    run_stages = g(["run_stages"], {}) or {}
    stage4 = g(["stage4"], {}) or {}
    stage5 = g(["stage5"], {}) or {}

    return SessionParams(
        cancer_code=g(["cancer_code"]),
        input_dir=Path(g(["paths","counts_dir"])),
        output_dir=Path(g(["paths","output_dir"])),
        sample_sheet_tsv=Path(g(["paths","sample_sheet_tsv"])),
        clinical_tsv=Path(g(["paths","clinical_tsv"])),
        subtype_tsv=Path(g(["paths","subtype_tsv"])),
        subtype_sample_id_col=g(["subtype","sample_id_column"]),
        subtype_value_col=g(["subtype","subtype_value_column"]),
        aggressive_subtypes=g(["subtype","aggressive_subtypes"]) or [],
        tumor_sample_types=g(["sample_types","tumor"]) or [],
        control_sample_types=g(["sample_types","control"]) or [],
        race_1=g(["races","race_1"]),
        race_2=g(["races","race_2"]),
        all_races_label=g(["races","all_label"], "All"),
        run_stages={
            "stage1": bool(run_stages.get("stage1", True)),
            "stage2": bool(run_stages.get("stage2", True)),
            "stage3": bool(run_stages.get("stage3", True)),
            "stage4": bool(run_stages.get("stage4", False)),
            "stage5": bool(run_stages.get("stage5", False)),
        },
        stage4=stage4,
        stage5=stage5
    )
