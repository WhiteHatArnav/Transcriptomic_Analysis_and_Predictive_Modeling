# Load required libraries
library(DESeq2)
library(biomaRt)
library(EnhancedVolcano)
library(dplyr)

# Convert DESeqResults to a regular data frame
luad_resAll_df <- as.data.frame(luad_resAll)

# Strip version numbers
luad_resAll_df$ensembl_id <- sub("\\..*", "", rownames(luad_resAll_df))

# Get gene names using biomaRt
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = unique(luad_resAll_df$ensembl_id),
                  mart = ensembl)

# Join gene names to the results
luad_resAll_df <- luad_resAll_df %>%
  left_join(gene_map, by = c("ensembl_id" = "ensembl_gene_id"))

# Fill missing gene names with Ensembl IDs
luad_resAll_df$hgnc_symbol[is.na(luad_resAll_df$hgnc_symbol)] <- 
  luad_resAll_df$ensembl_id[is.na(luad_resAll_df$hgnc_symbol)]

# Plot
EnhancedVolcano(luad_resAll_df,
                lab = luad_resAll_df$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'All Races LUAD Volcano Plot')
