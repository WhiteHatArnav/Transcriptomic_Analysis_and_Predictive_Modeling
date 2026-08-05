# Race-Stratified Transcriptomic Analysis and Prognostic Modeling Across Multiple Cancers

## Overview

This repository contains the analysis scripts developed as part of the doctoral dissertation together with supplementary datasets generated from the completed analytical workflows. The base directory includes *Supplementary_Sample_ID_Reference_List.csv*, containing the complete list of TCGA sample IDs used for each cancer cohort to facilitate reproduction and independent validation of all analyses. It also includes *Supplementary_DEG_Gene_Symbols.csv*, a consolidated list of differentially expressed genes reported as HGNC gene symbols across all cancer, race, and regulation groups, and *Supplementary_Functional_Enrichment.xlsx*, a compiled workbook containing Gene Ontology (Biological Process, Cellular Component, and Molecular Function) and KEGG pathway enrichment results organized by cancer, race, regulation status, and functional category.

**Race-Stratified Transcriptomic Analysis and Prognostic Modeling Across Multiple Cancers**

The project investigates transcriptomic determinants of prognosis across four major cancers through differential expression analysis, functional enrichment analysis, principal component analysis, and predictive modeling using LASSO-regularized Cox proportional hazards models. Particular emphasis is placed on race-stratified transcriptomic analysis, prognostic modeling, reproducible computational workflows, and interpretation of high-dimensional gene expression datasets.

## Cancer Cohorts

- Breast Invasive Carcinoma (BRCA)
- Lung Adenocarcinoma (LUAD)
- Liver Hepatocellular Carcinoma (LIHC)
- Prostate Adenocarcinoma (PRAD)

## Analytical Workflow

1. Acquisition and integration of TCGA transcriptomic and clinical datasets
2. Cancer-specific subtype selection and demographic stratification
3. Preparation and filtering of transcriptomic expression matrices
4. Differential expression analysis using DESeq2
5. Gene Ontology (GO) and KEGG pathway enrichment analysis
6. Comparative enrichment analysis across demographic cohorts
7. Parent-term consolidation and enrichment comparison workflows
8. Principal component analysis (PCA) of cohort structure
9. LASSO-CV feature selection
10. Cox proportional hazards survival modeling
11. Risk score generation and Kaplan–Meier evaluation
12. K-fold cross-validation and 70/30 holdout validation strategies
13. Recurrence-based predictive modeling for selected endpoints
14. Development of a generalized configurable analysis pipeline

## Repository Structure

### 01_BRCAPrimaryScripts
Primary BRCA transcriptomic, enrichment, and survival modeling analyses.

### 02_LUADPrimaryScripts
Primary LUAD transcriptomic, enrichment, and survival modeling analyses.

### 03_LIHCPrimaryScripts
Primary LIHC transcriptomic, enrichment, and survival modeling analyses.

### 04_PRADPrimaryScripts
Primary PRAD transcriptomic, enrichment, and recurrence-focused predictive analyses.

### 05_KFoldPredictiveAnalysisScriptsAndOutputs
K-fold LASSO-Cox validation workflows and associated outputs.

### 06_RecurrencePrediction70_30Scripts
Recurrence endpoint predictive modeling workflows using 70/30 train-test validation.

### 07_AnalyticalKMCurvesScripts
Scripts used to generate analytical Kaplan–Meier survival evaluation plots.

### 08_PCAScriptsandOutputs
Principal component analysis workflows and associated diagnostic outputs.

### 09_GO_PE_TermRaceSetAnalysisScripts
Race-stratified Gene Ontology and pathway enrichment comparison analyses.

### 10_GO_PE_ParentTermAnalysisScripts
Parent-term consolidation and semantic reduction workflows for Gene Ontology enrichment interpretation.

### 11_DissertationTemplatePipeline
Generalized configuration-driven transcriptomic analysis pipeline developed from the dissertation methodology.

## Data Sources

- The Cancer Genome Atlas (TCGA)
- National Cancer Institute Genomic Data Commons (GDC)
- UCSC Xena Browser subtype annotations where applicable

Raw TCGA datasets are not distributed within this repository and should be obtained directly from the original data providers using the sample ID list available in the base folder of this repository.

## Software Environment

- Python
- R

Major analytical components include DESeq2, DAVID, Kaplan–Meier survival analysis, Cox proportional hazards modeling, LASSO regularization, and principal component analysis.


**Author:** Arnav Joshi
**Program:** Doctoral Program in Computational Science
**Supported by:** Cancer Prevention and Research Institute of Texas (CPRIT)
**Institution:** The University of Texas at El Paso
**Year:** 2026
