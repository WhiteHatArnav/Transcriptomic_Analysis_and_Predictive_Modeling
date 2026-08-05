# Load required libraries
library(biomaRt)
library(dplyr)
library(readr)

# === 1️⃣ FILE PATHS ===
files <- c(
  "/Users/arnavjoshi/Desktop/LUAD_GOandPE/WhiteRace_Downregulated_Genes.csv",
  "/Users/arnavjoshi/Desktop/LUAD_GOandPE/WhiteRace_Upregulated_Genes.csv",
  "/Users/arnavjoshi/Desktop/LUAD_GOandPE/AllRace_Downregulated_Genes.csv",
  "/Users/arnavjoshi/Desktop/LUAD_GOandPE/AllRace_Upregulated_Genes.csv"
)

# === 2️⃣ CONNECT TO ENSEMBL ===
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# === 3️⃣ FUNCTION TO PROCESS ONE FILE ===
convert_and_clean <- function(file_path) {
  
  # Read CSV
  data <- read_csv(file_path, show_col_types = FALSE)
  
  cat("\nProcessing file:", basename(file_path), "\n")
  
  # Column containing Ensembl IDs
  ensembl_col <- "gene_id"
  
  # Get unique Ensembl IDs
  ensembl_ids <- unique(na.omit(data[[ensembl_col]]))
  
  # Convert Ensembl → HGNC symbols
  mapping <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
  )
  
  # Merge and filter out blanks
  clean_data <- data %>%
    left_join(mapping, by = setNames("ensembl_gene_id", ensembl_col)) %>%
    filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
    select(hgnc_symbol)  # Only keep gene symbol column
  
  # Save output file with "_clean_genesymbols" suffix
  output_file <- sub("\\.csv$", "_clean_genesymbols.csv", file_path)
  write_csv(clean_data, output_file)
  
  # Print mapping stats
  total_ids <- nrow(data)
  mapped_ids <- nrow(clean_data)
  percent_mapped <- round(100 * mapped_ids / total_ids, 2)
  
  cat("Total Ensembl IDs:", total_ids, "\n")
  cat("Mapped to Gene Symbols:", mapped_ids, "\n")
  cat("Mapping Success Rate:", percent_mapped, "%\n")
  cat("Clean gene symbol file saved as:", output_file, "\n")
}

# === 4️⃣ PROCESS ALL FILES ===
lapply(files, convert_and_clean)
