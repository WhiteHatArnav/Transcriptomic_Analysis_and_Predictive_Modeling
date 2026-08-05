# -----------------------------
# 1. Install & load packages
# -----------------------------
if (!requireNamespace("UpSetR", quietly = TRUE)) install.packages("UpSetR")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("grid", quietly = TRUE)) install.packages("grid")

library(UpSetR)
library(tidyverse)
library(grid)

# -----------------------------
# 2. Define paths
# -----------------------------
summary_dir <- "/Users/arnavjoshi/Desktop/GOParentTermAnalysis/SummaryCSVs"
output_dir  <- "/Users/arnavjoshi/Desktop/GOParentTermAnalysis/UpSetPlots"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -----------------------------
# 3. Helpers
# -----------------------------
clean_col_terms <- function(x) {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & x != ""]
  unique(x)
}

# Tokenize a column name robustly (underscores are treated as separators)
name_tokens <- function(colname) {
  n <- tolower(colname)
  toks <- unlist(strsplit(n, "[^a-z0-9]+"))  # splits on underscore and any non-alnum
  toks <- toks[toks != ""]
  toks
}

# Determine which of (All, White, Black) a "region column" belongs to
# Works for:
#   Unique_All, Unique_White, Unique_Black
#   Common_White_All, Common_All_Black, Common_White_Black
#   Common_All_Three, Common_All_White_Black, etc.
membership_from_colname <- function(colname) {
  toks <- name_tokens(colname)
  
  # If "three" is present, treat as intersection of all three
  if ("three" %in% toks || ("all" %in% toks && "three" %in% toks)) {
    return(c("All", "White", "Black"))
  }
  
  members <- character(0)
  if ("all"   %in% toks) members <- c(members, "All")
  if ("white" %in% toks) members <- c(members, "White")
  if ("black" %in% toks) members <- c(members, "Black")
  members
}

df_to_sets_all_white_black <- function(df) {
  # If CSV already contains direct columns All/White/Black, use them
  orig_names <- colnames(df)
  lower_names <- tolower(orig_names)
  
  idx_all   <- which(lower_names == "all")
  idx_white <- which(lower_names == "white")
  idx_black <- which(lower_names == "black")
  
  if (length(idx_all) == 1 && length(idx_white) == 1 && length(idx_black) == 1) {
    return(list(
      All   = clean_col_terms(df[[idx_all]]),
      White = clean_col_terms(df[[idx_white]]),
      Black = clean_col_terms(df[[idx_black]])
    ))
  }
  
  # Otherwise reconstruct from region columns
  sets <- list(All = character(0), White = character(0), Black = character(0))
  
  for (i in seq_along(orig_names)) {
    colname <- orig_names[i]
    members <- membership_from_colname(colname)
    if (length(members) == 0) next
    
    terms <- clean_col_terms(df[[i]])
    if (length(terms) == 0) next
    
    for (m in members) {
      sets[[m]] <- unique(c(sets[[m]], terms))
    }
  }
  
  sets
}

# -----------------------------
# 4. Generate plots
# -----------------------------
csv_files <- list.files(
  summary_dir,
  pattern = "_GO_ParentTermSummary\\.csv$",
  full.names = TRUE
)

for (csv_path in csv_files) {
  
  file_name <- basename(csv_path)
  parts <- strsplit(file_name, "_")[[1]]
  
  cancer     <- parts[1]
  regulation <- parts[2]
  ontology   <- parts[3]
  
  plot_title <- paste(cancer, "-", regulation, ontology, "GO Parent Term Commonality")
  
  output_file <- paste0(
    output_dir, "/",
    cancer, "_", regulation, "_", ontology,
    "_UpSet_GO_ParentTerms.png"
  )
  
  df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  sets_list <- df_to_sets_all_white_black(df)
  
  # Skip if truly empty
  if (length(sets_list$All) == 0 && length(sets_list$White) == 0 && length(sets_list$Black) == 0) {
    message("Skipping (no terms found): ", csv_path)
    next
  }
  
  png(filename = output_file, width = 3000, height = 2000, res = 300)
  
  par(mar = c(6, 10, 4, 4))
  
  suppressWarnings({
    print(
      upset(
        fromList(list(All = sets_list$All, White = sets_list$White, Black = sets_list$Black)),
        sets = c("All", "White", "Black"),
        keep.order = TRUE,
        order.by = "freq",
        nintersects = NA,
        sets.bar.color = "blue",
        main.bar.color = "red",
        text.scale = c(2.2, 2.2, 1.8, 1.8, 2.2, 1.8),
        mainbar.y.label = "Number of Parent GO Terms",
        sets.x.label = "Parent GO Terms per Set"
      )
    )
  })
  
  grid.text(
    plot_title,
    x = 0.5,
    y = 0.95,
    gp = gpar(fontsize = 24, fontface = "bold")
  )
  
  dev.off()
  
  message("Saved UpSet plot: ", output_file)
}
