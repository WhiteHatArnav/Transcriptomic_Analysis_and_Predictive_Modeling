# Load required libraries
library(DESeq2)
library(biomaRt)
library(EnhancedVolcano)
library(dplyr)

# Convert DESeqResults to a regular data frame
lihc_resWhite_df <- as.data.frame(lihc_resWhite)

# Strip version numbers
lihc_resWhite_df$ensembl_id <- sub("\\..*", "", rownames(lihc_resWhite_df))

# Get gene names using biomaRt
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = unique(lihc_resWhite_df$ensembl_id),
                  mart = ensembl)

# Join gene names to the results
lihc_resWhite_df <- lihc_resWhite_df %>%
  left_join(gene_map, by = c("ensembl_id" = "ensembl_gene_id"))

# Fill missing gene names with Ensembl IDs
lihc_resWhite_df$hgnc_symbol[is.na(lihc_resWhite_df$hgnc_symbol)] <- 
  lihc_resWhite_df$ensembl_id[is.na(lihc_resWhite_df$hgnc_symbol)]

# Plot
EnhancedVolcano(lihc_resWhite_df,
                lab = lihc_resWhite_df$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'White Race LIHC Volcano Plot')
