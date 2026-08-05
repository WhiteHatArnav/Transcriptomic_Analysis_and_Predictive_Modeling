# Dissertation Pipeline Template (Stages 1–5)

## Overview
This repository provides a reproducible, configuration driven pipeline for TCGA style transcriptomic analysis across cancers and demographic cohorts. The pipeline is designed to support large scale RNA sequencing studies with consistent preprocessing, statistical modeling, and predictive analytics while remaining flexible to cancer type and cohort definition.

All pipeline behavior is controlled through a single `config.yaml` file, allowing the same codebase to be reused across datasets and studies without modification.

## Pipeline Stages

### Stage 1: Expression Matrix Construction
Raw count files are aggregated into a unified expression matrix. Sample level metadata are integrated at this stage, including sample identifiers and alignment quality metrics. The resulting matrix forms the foundation for all downstream analyses.

### Stage 2: Cohort Extraction
Samples are stratified into cohorts based on aggressive tumor subtypes and control definitions. Race specific cohorts and an all samples cohort are generated. This stage produces cohort specific expression matrices and corresponding metadata files used by later stages.

### Stage 3: Differential Expression Analysis
Differential expression analysis is performed using DESeq2 through an R backend. For each cohort, results include log fold changes, statistical significance values, and volcano plot visualizations.

### Stage 4: Functional Enrichment Analysis
Gene Ontology and pathway enrichment analysis is conducted using g:Profiler for Biological Process, Cellular Component, Molecular Function, and KEGG pathways. Results are summarized as dot plots and exported as structured CSV files for downstream interpretation.

### Stage 5: Predictive Analytics with LASSO–Cox
Predictive modeling is performed using LASSO regularized Cox proportional hazards models. Two modeling strategies are supported:
- A 70:30 train test split
- K fold cross validated LASSO–Cox modeling

The pipeline includes safeguards for small cohort sizes, LASSO non convergence, and Cox model fitting failures. When modeling assumptions are not satisfied, affected cohorts are skipped without terminating the pipeline.

## Running the Pipeline

```
python src/main.py --config config.yaml
```

All outputs are written under the directory specified by `output_dir` in the configuration file, organized by stage.


## Design Principles
- Configuration driven execution with no hard coded paths
- Explicit cohort level outputs for traceability
- Graceful handling of small sample sizes and modeling failures
- Modular stage separation to support partial pipeline execution
- Compatibility with TCGA and Xena style survival annotations

This template is intended as a foundation for cancer transcriptomics research and predictive modeling, and can be extended or adapted for additional analyses as needed.
