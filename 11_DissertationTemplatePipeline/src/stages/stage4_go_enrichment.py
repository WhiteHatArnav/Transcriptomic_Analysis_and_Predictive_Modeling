import pandas as pd
from pathlib import Path
import importlib.util

from pipeline.session import SessionParams
from pipeline.state import PipelineState
from pipeline.runner_r import run_rscript


# ------------------------------------------------------------------
# Robust g:Profiler loader (handles upstream packaging inconsistencies)
# ------------------------------------------------------------------
def _load_gprofiler():
    for name in ("gprofiler", "gprofiler_official"):
        spec = importlib.util.find_spec(name)
        if spec is not None:
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            return module.GProfiler
    raise ImportError(
        "gprofiler-official is installed but no importable module was found."
    )


GProfiler = _load_gprofiler()
# ------------------------------------------------------------------


GP_CATEGORY_MAP = {
    "GO:BP": "BPDirect",
    "GO:CC": "CCDirect",
    "GO:MF": "MFDirect",
    "KEGG": "KEGGPathway"
}


# ------------------------------------------------------------------
# Utility helpers
# ------------------------------------------------------------------
def _strip_ensembl_version(gene_ids):
    """Remove Ensembl version suffixes (ENSGxxxx.xx → ENSGxxxx)."""
    return [g.split(".")[0] for g in gene_ids]


def _extract_gene_lists(deseq_df, log2fc_thresh, pval_thresh):
    """
    Extract significantly up- and down-regulated genes using:
    - padj (only; no fallback to pvalue)
    - absolute log2FC threshold
    """

    sig = deseq_df.dropna(subset=["log2FoldChange", "padj"])
    sig = sig[sig["padj"] <= pval_thresh]

    up = sig[sig["log2FoldChange"] >= log2fc_thresh].index.astype(str)
    down = sig[sig["log2FoldChange"] <= -log2fc_thresh].index.astype(str)

    # IMPORTANT: strip Ensembl version suffixes
    up = _strip_ensembl_version(up)
    down = _strip_ensembl_version(down)

    print(f"    Selected {len(up)} up / {len(down)} down genes")

    return up, down


def _format_gprofiler(res, gene_list):
    """
    Normalize g:Profiler output across schema versions and
    produce a REVIGO-compatible enrichment table.
    """

    if res is None or res.empty:
        return pd.DataFrame()

    # --- Term ID handling ---
    if "term_id" in res.columns:
        term_id_col = res["term_id"]
    elif "native" in res.columns:
        term_id_col = res["native"]
    else:
        raise RuntimeError("No GO / pathway identifier column found.")

    # --- Gene membership handling (schema-robust) ---
    if "intersection_names" in res.columns:
        genes_col = res["intersection_names"].apply(
            lambda x: ",".join(x) if isinstance(x, list) else str(x)
        )
    elif "intersections" in res.columns:
        genes_col = res["intersections"].apply(str)
    elif "intersection" in res.columns:
        genes_col = res["intersection"].apply(
            lambda x: ",".join(x) if isinstance(x, list) else str(x)
        )
    elif "hits" in res.columns:
        genes_col = res["hits"].apply(
            lambda x: ",".join(x) if isinstance(x, list) else str(x)
        )
    else:
        # Valid case: enrichment without gene list
        genes_col = pd.Series([""] * len(res), index=res.index)

    return pd.DataFrame({
        "Category": res["source"].map(GP_CATEGORY_MAP),
        "Term": term_id_col.astype(str) + "~" + res["name"].astype(str),
        "Genes": genes_col,
        "Count": res["intersection_size"],
        "P-Value": res["p_value"],
        "Fold Enrichment": (
            res["intersection_size"] / max(1, len(gene_list))
        ) / (
            res["term_size"] / res["effective_domain_size"]
        )
    })


# ------------------------------------------------------------------
# Stage 4 main entry point
# ------------------------------------------------------------------
def run_stage4(session: SessionParams, state: PipelineState):
    outdir = session.output_dir / "stage4"
    outdir.mkdir(parents=True, exist_ok=True)

    gp = GProfiler(return_dataframe=True)

    sources = session.stage4.get(
        "enrichment_sources", ["GO:BP", "GO:CC", "GO:MF", "KEGG"]
    )
    pcut = float(session.stage4.get("pvalue_cutoff", 0.05))
    l2fc = float(session.stage4.get("log2fc_threshold", 1.0))

    for cohort in [session.race_1, session.race_2, session.all_races_label]:
        print(f"\nStage 4 enrichment for {cohort}")

        deseq_key = f"deseq_{cohort}"
        deseq_path = state.get(deseq_key)

        deseq_df = pd.read_csv(deseq_path, index_col=0)

        up_genes, down_genes = _extract_gene_lists(deseq_df, l2fc, pcut)

        for label, genes in {
            "UpRegulated": up_genes,
            "DownRegulated": down_genes
        }.items():

            if not genes:
                print(f"    Skipping {label} (0 genes)")
                continue

            res = gp.profile(
                organism=session.stage4.get("organism", "hsapiens"),
                query=genes,
                sources=sources
            )

            formatted = _format_gprofiler(res, genes)
            if formatted.empty:
                print(f"    No enrichment for {label}")
                continue

            for cat in formatted["Category"].dropna().unique():
                out_csv = outdir / f"{session.cancer_code}_{cohort}_{label}_{cat}.csv"
                formatted[formatted["Category"] == cat].to_csv(out_csv, index=False)
                print(f"    Wrote {out_csv.name}")

        # ---- Dotplots ----
        csvs = sorted(outdir.glob(f"{session.cancer_code}_{cohort}_*.csv"))
        if not csvs:
            print("    No CSVs → skipping dotplots")
            continue

        base = Path(__file__).resolve().parents[2]
        pdf_out = outdir / f"{session.cancer_code}_{cohort}_GO_PE_Dotplots.pdf"

        run_rscript(
            base / "R" / "stage4_dotplots.R",
            [
                ",".join(str(p) for p in csvs),
                str(pdf_out),
                session.cancer_code
            ]
        )

        print(f"    Dotplots written to {pdf_out}")
