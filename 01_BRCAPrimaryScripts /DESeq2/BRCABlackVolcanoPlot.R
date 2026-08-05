# library(EnhancedVolcano)
# 
# EnhancedVolcano(brca_resBlack,
#                 lab = rep('', nrow(brca_resBlack)),
#                 x = 'log2FoldChange',
#                 y = 'pvalue',
#                 title = 'Black Race BRCA Volcano Plot')

# Convert DESeqResults to a regular data frame
brca_resBlack_df <- as.data.frame(brca_resBlack)

# Strip version numbers from Ensembl IDs (from rownames)
brca_resBlack_df$ensembl_id <- sub("\\..*", "", rownames(brca_resBlack_df))

# Get gene names using biomaRt
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = unique(brca_resBlack_df$ensembl_id),
                  mart = ensembl)

# Load dplyr
library(dplyr)

# Join gene names to the results
brca_resBlack_df <- brca_resBlack_df %>%
  left_join(gene_map, by = c("ensembl_id" = "ensembl_gene_id"))

# Fill missing gene names with Ensembl IDs
brca_resBlack_df$hgnc_symbol[is.na(brca_resBlack_df$hgnc_symbol)] <- brca_resBlack_df$ensembl_id[is.na(brca_resBlack_df$hgnc_symbol)]

# Plot
library(EnhancedVolcano)

EnhancedVolcano(brca_resBlack_df,
                lab = brca_resBlack_df$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Black Race BRCA Volcano Plot')
