# Load required libraries
library(DESeq2)
library(biomaRt)
library(EnhancedVolcano)
library(dplyr)

# Convert DESeqResults to a regular data frame
luad_resBlack_df <- as.data.frame(luad_resBlack)

# Strip version numbers
luad_resBlack_df$ensembl_id <- sub("\\..*", "", rownames(luad_resBlack_df))

# Get gene names using biomaRt
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = unique(luad_resBlack_df$ensembl_id),
                  mart = ensembl)

# Join gene names to the results
luad_resBlack_df <- luad_resBlack_df %>%
  left_join(gene_map, by = c("ensembl_id" = "ensembl_gene_id"))

# Fill missing gene names with Ensembl IDs
luad_resBlack_df$hgnc_symbol[is.na(luad_resBlack_df$hgnc_symbol)] <- 
  luad_resBlack_df$ensembl_id[is.na(luad_resBlack_df$hgnc_symbol)]

# Plot
EnhancedVolcano(luad_resBlack_df,
                lab = luad_resBlack_df$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Black Race LUAD Volcano Plot')
