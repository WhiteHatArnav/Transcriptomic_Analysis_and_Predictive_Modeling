# library(EnhancedVolcano)
# 
# EnhancedVolcano(brca_resWhite,
#                 lab = rep('', nrow(brca_resWhite)),
#                 x = 'log2FoldChange',
#                 y = 'pvalue',
#                 title = 'White Race BRCA Volcano Plot')


# Step 1: Install BiocManager (if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Step 2: Use BiocManager to install biomaRt
BiocManager::install("biomaRt")

# Load required libraries
library(DESeq2)
library(biomaRt)
library(EnhancedVolcano)
library(dplyr)

# Convert DESeqResults to a regular data frame
brca_resWhite_df <- as.data.frame(brca_resWhite)

# Strip version numbers
brca_resWhite_df$ensembl_id <- sub("\\..*", "", rownames(brca_resWhite_df))

# Get gene names using biomaRt
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = unique(brca_resWhite_df$ensembl_id),
                  mart = ensembl)

# Load dplyr
library(dplyr)

# Join gene names to the results
brca_resWhite_df <- brca_resWhite_df %>%
  left_join(gene_map, by = c("ensembl_id" = "ensembl_gene_id"))

# Fill missing gene names with Ensembl IDs
brca_resWhite_df$hgnc_symbol[is.na(brca_resWhite_df$hgnc_symbol)] <- brca_resWhite_df$ensembl_id[is.na(brca_resWhite_df$hgnc_symbol)]

# Plot
library(EnhancedVolcano)

EnhancedVolcano(brca_resWhite_df,
                lab = brca_resWhite_df$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'White Race BRCA Volcano Plot')