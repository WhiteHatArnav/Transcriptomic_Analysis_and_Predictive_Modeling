# Load required libraries
library(DESeq2)
library(biomaRt)
library(EnhancedVolcano)
library(dplyr)

# Convert DESeqResults to a regular data frame
luad_resWhite_df <- as.data.frame(luad_resWhite)

# Strip version numbers
luad_resWhite_df$ensembl_id <- sub("\\..*", "", rownames(luad_resWhite_df))

# Get gene names using biomaRt
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = unique(luad_resWhite_df$ensembl_id),
                  mart = ensembl)

# Join gene names to the results
luad_resWhite_df <- luad_resWhite_df %>%
  left_join(gene_map, by = c("ensembl_id" = "ensembl_gene_id"))

# Fill missing gene names with Ensembl IDs
luad_resWhite_df$hgnc_symbol[is.na(luad_resWhite_df$hgnc_symbol)] <- 
  luad_resWhite_df$ensembl_id[is.na(luad_resWhite_df$hgnc_symbol)]

# Plot
EnhancedVolcano(luad_resWhite_df,
                lab = luad_resWhite_df$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'White Race LUAD Volcano Plot')
