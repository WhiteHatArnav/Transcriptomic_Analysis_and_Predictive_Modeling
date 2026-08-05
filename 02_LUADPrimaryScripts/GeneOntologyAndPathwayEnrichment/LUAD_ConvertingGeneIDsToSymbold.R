# Load required library
library(biomaRt)
library(dplyr)
library(readr)

# Define input file paths
files <- c(
  "/Users/arnavjoshi/Desktop/LUAD_GOandPE/WhiteRace_Downregulated_Genes.csv",
  "/Users/arnavjoshi/Desktop/LUAD_GOandPE/WhiteRace_Upregulated_Genes.csv",
  "/Users/arnavjoshi/Desktop/LUAD_GOandPE/AllRace_Downregulated_Genes.csv",
  "/Users/arnavjoshi/Desktop/LUAD_GOandPE/AllRace_Upregulated_Genes.csv"
)

# Connect to Ensembl (human)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Function to convert one file
convert_ensembl_to_symbol <- function(file_path) {
  # Read CSV (assuming it has a column named 'gene_id')
  data <- read_csv(file_path, show_col_types = FALSE)
  
  if (!"gene_id" %in% names(data)) {
    stop(paste("File", basename(file_path), "does not have a column named 'gene_id'"))
  }
  
  # Get unique Ensembl IDs
  ensembl_ids <- unique(na.omit(data$gene_id))
  
  # Convert Ensembl → HGNC symbol
  conversion <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = mart
  )
  
  # Merge back to original data
  merged <- data %>%
    left_join(conversion, by = c("gene_id" = "ensembl_gene_id"))
  
  # Define output file path
  output_file <- sub(".csv$", "_WithSymbols.csv", file_path)
  
  # Save output
  write_csv(merged, output_file)
  message("Saved: ", output_file)
}

# Run conversion for all files
lapply(files, convert_ensembl_to_symbol)
