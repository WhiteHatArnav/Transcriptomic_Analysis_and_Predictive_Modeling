
import os
import pandas as pd
from tqdm import tqdm
from pipeline.session import SessionParams

UNWANTED = {"N_unmapped","N_multimapping","N_noFeature","N_ambiguous"}

def trim_tcga_sample_id(s):
    return s[:-1] if s else s

def run_stage1(session: SessionParams):
    outdir = session.output_dir / "stage1"
    outdir.mkdir(parents=True, exist_ok=True)

    sample_df = pd.read_csv(session.sample_sheet_tsv, sep="\t", dtype=str)
    file_to_case = dict(zip(sample_df["File Name"], sample_df["Case ID"]))
    file_to_sample = dict(zip(sample_df["File Name"], sample_df["Sample ID"]))
    sample_to_type = dict(zip(sample_df["Sample ID"], sample_df["Tissue Type"]))

    clinical_df = pd.read_csv(session.clinical_tsv, sep="\t", dtype=str)
    case_to_race = dict(zip(clinical_df["cases.submitter_id"], clinical_df["demographic.race"]))

    subtype_df = pd.read_csv(session.subtype_tsv, sep="\t", dtype=str)
    sample_to_subtype = dict(zip(subtype_df[session.subtype_sample_id_col], subtype_df[session.subtype_value_col]))

    expr, sample_ids, genes = {}, {}, None

    folders = [f for f in os.listdir(session.input_dir) if (session.input_dir / f).is_dir()]
    for folder in tqdm(folders, desc="Stage 1: building expression matrix"):
        fpath = session.input_dir / folder
        tsvs = [f for f in os.listdir(fpath) if f.endswith(".tsv")]
        if not tsvs:
            continue
        tsv = tsvs[0]
        if tsv not in file_to_case:
            continue

        df = pd.read_csv(fpath / tsv, sep="\t", comment="#", usecols=["gene_id","unstranded"])
        df["gene_id"] = df["gene_id"].str.split("|").str[0]

        if genes is None:
            genes = df["gene_id"].tolist()

        cid = file_to_case[tsv]
        expr[cid] = df["unstranded"].values
        sample_ids[cid] = file_to_sample[tsv]

    mat = pd.DataFrame(expr, index=genes)
    mat.index.name = "Gene ID"

    mat = pd.concat([pd.Series(sample_ids, name="Sample ID").to_frame().T, mat])
    mat = mat.loc[~mat.index.isin(UNWANTED)]

    sids = mat.loc["Sample ID"]

    mat = pd.concat([
        mat.loc[["Sample ID"]],
        pd.Series([case_to_race.get(cid,"NA") for cid in mat.columns], index=mat.columns, name="Race").to_frame().T,
        pd.Series([sample_to_type.get(s,"NA") for s in sids], index=mat.columns, name="Sample Type").to_frame().T,
        pd.Series([sample_to_subtype.get(trim_tcga_sample_id(s),"NA") for s in sids], index=mat.columns, name="Sample Subtype").to_frame().T,
        mat.drop(index=["Sample ID"])
    ])

    out = outdir / f"{session.cancer_code}_expression_matrix_with_metadata.csv"
    mat.to_csv(out)
    return out
