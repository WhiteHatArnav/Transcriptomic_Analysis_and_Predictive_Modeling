# Load required libraries
library(DESeq2)
library(biomaRt)
library(EnhancedVolcano)
library(dplyr)

# Convert DESeqResults to a regular data frame
prad_resWhite_df <- as.data.frame(prad_resWhite)

# Strip version numbers from Ensembl IDs
prad_resWhite_df$ensembl_id <- sub("\\..*", "", rownames(prad_resWhite_df))

# Get gene names using biomaRt
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = unique(prad_resWhite_df$ensembl_id),
  mart = ensembl
)

# Join gene names to the DESeq results
prad_resWhite_df <- prad_resWhite_df %>%
  left_join(gene_map, by = c("ensembl_id" = "ensembl_gene_id"))

# Fill missing gene names with Ensembl IDs
prad_resWhite_df$hgnc_symbol[is.na(prad_resWhite_df$hgnc_symbol)] <-
  prad_resWhite_df$ensembl_id[is.na(prad_resWhite_df$hgnc_symbol)]

# Plot volcano plot
EnhancedVolcano(prad_resWhite_df,
                lab = prad_resWhite_df$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'White Race PRAD Volcano Plot')
