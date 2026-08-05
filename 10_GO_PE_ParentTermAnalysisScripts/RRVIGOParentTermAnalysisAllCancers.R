if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("rrvgo")


suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

# -------------------------
# Package setup (install if missing)
# -------------------------
pkg_needed <- c("rrvgo", "GOSemSim", "org.Hs.eg.db")
pkg_missing <- pkg_needed[!sapply(pkg_needed, requireNamespace, quietly = TRUE)]
if (length(pkg_missing) > 0) {
  install.packages(pkg_missing, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(rrvgo)
  library(GOSemSim)
  library(org.Hs.eg.db)
})

BiocManager::install("org.Hs.eg.db")
library(readr)
library(dplyr)
library(stringr)
library(rrvgo)
library(GOSemSim)
library(org.Hs.eg.db)
# -------------------------
# CONFIG
# -------------------------
CANCERS <- c("BRCA", "LUAD", "PRAD", "LIHC")
REGS <- c("UpReg", "DownReg")

# Only GO ontologies (rrvgo/REVIGO does not do KEGG)
ONTOLOGIES <- c("BP", "CC", "MF")

INPUT_BASE_DIR <- "/Users/arnavjoshi/Desktop/GOPETermsComparison"
OUTPUT_BASE_DIR <- "/Users/arnavjoshi/Desktop/GOParentTermAnalysis"

# Columns produced by your term-comparison CSVs
CATEGORY_COLUMNS <- c(
  "In_White",
  "In_Black",
  "In_All",
  "Common_AllThree",
  "Shared_White_Black",
  "Shared_White_All",
  "Shared_Black_All",
  "Unique_White",
  "Unique_Black",
  "Unique_All"
)

# rrvgo settings
SIM_METHOD <- "Rel"     # semantic similarity method used by rrvgo
THRESHOLD  <- 0.70      # reduction threshold (typical values: 0.5 to 0.9)

# -------------------------
# HELPERS
# -------------------------
extract_go_ids <- function(x) {
  # Extract GO:####### anywhere in the string (works for "GO:...~name", "GO:..." alone, etc.)
  ids <- str_extract(x, "GO:\\d{7}")
  ids <- ids[!is.na(ids) & ids != ""]
  unique(ids)
}

safe_write_csv <- function(df, out_path) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  write_csv(df, out_path)
}

run_rrvgo_for_terms <- function(go_ids, ontology, threshold = THRESHOLD) {
  # rrvgo needs at least 2 terms to compute a similarity matrix meaningfully
  go_ids <- unique(go_ids)
  if (length(go_ids) < 2) {
    return(NULL)
  }
  
  # Scores are required; your comparison CSVs do not include p-values.
  # Use equal scores so rrvgo can pick representatives deterministically.
  scores <- setNames(rep(1, length(go_ids)), go_ids)
  
  sim_matrix <- calculateSimMatrix(
    go_ids,
    orgdb  = "org.Hs.eg.db",
    ont    = ontology,
    method = SIM_METHOD
  )
  
  reduced <- reduceSimMatrix(
    sim_matrix,
    scores,
    threshold = threshold,
    orgdb = "org.Hs.eg.db"
  )
  
  # reduced typically includes: go, term, parent, parentTerm, score, size, members, etc.
  reduced
}

# -------------------------
# MAIN
# -------------------------
dir.create(OUTPUT_BASE_DIR, recursive = TRUE, showWarnings = FALSE)

for (cancer in CANCERS) {
  for (reg in REGS) {
    for (ont in ONTOLOGIES) {
      
      in_file <- file.path(
        INPUT_BASE_DIR,
        cancer,
        reg,
        paste0(cancer, "_", reg, "_", ont, "_TermComparison.csv")
      )
      
      if (!file.exists(in_file)) {
        message("Missing input file, skipping: ", in_file)
        next
      }
      
      df <- suppressMessages(read_csv(in_file, show_col_types = FALSE))
      
      missing_cols <- setdiff(CATEGORY_COLUMNS, colnames(df))
      if (length(missing_cols) > 0) {
        message("File is missing expected columns, skipping: ", in_file)
        message("Missing: ", paste(missing_cols, collapse = ", "))
        next
      }
      
      for (cat_col in CATEGORY_COLUMNS) {
        raw_terms <- df[[cat_col]]
        raw_terms <- raw_terms[!is.na(raw_terms) & raw_terms != ""]
        
        go_ids <- extract_go_ids(raw_terms)
        
        reduced <- run_rrvgo_for_terms(go_ids, ontology = ont, threshold = THRESHOLD)
        
        out_dir <- file.path(OUTPUT_BASE_DIR, cancer, reg, ont)
        out_file <- file.path(out_dir, paste0(cat_col, "_ParentTermSummary.csv"))
        
        if (is.null(reduced)) {
          # Write an empty-but-informative CSV for traceability
          empty_df <- tibble(
            Category = cat_col,
            Ontology = ont,
            N_GO_Terms = length(go_ids),
            Note = "Fewer than 2 GO terms; rrvgo reduction not performed."
          )
          safe_write_csv(empty_df, out_file)
          next
        }
        
        # Add some helpful metadata columns
        reduced_out <- reduced %>%
          mutate(
            Category = cat_col,
            Ontology = ont,
            Cancer = cancer,
            Regulation = reg,
            N_Input_GO_Terms = length(go_ids)
          ) %>%
          relocate(Cancer, Regulation, Ontology, Category, N_Input_GO_Terms)
        
        safe_write_csv(reduced_out, out_file)
      }
      
      message("Completed: ", cancer, " ", reg, " ", ont)
    }
  }
}

message("All done. Outputs written under: ", OUTPUT_BASE_DIR)
