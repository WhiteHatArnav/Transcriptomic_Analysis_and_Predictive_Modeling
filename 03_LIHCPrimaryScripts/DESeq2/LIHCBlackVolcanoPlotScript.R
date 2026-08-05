# Convert DESeqResults to a regular data frame
lihc_resBlack_df <- as.data.frame(lihc_resBlack)

# Strip version numbers from Ensembl IDs (from rownames)
lihc_resBlack_df$ensembl_id <- sub("\\..*", "", rownames(lihc_resBlack_df))

# Get gene names using biomaRt
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = unique(lihc_resBlack_df$ensembl_id),
                  mart = ensembl)

# Load dplyr
library(dplyr)

# Join gene names to the results
lihc_resBlack_df <- lihc_resBlack_df %>%
  left_join(gene_map, by = c("ensembl_id" = "ensembl_gene_id"))

# Fill missing gene names with Ensembl IDs
lihc_resBlack_df$hgnc_symbol[is.na(lihc_resBlack_df$hgnc_symbol)] <- 
  lihc_resBlack_df$ensembl_id[is.na(lihc_resBlack_df$hgnc_symbol)]

# Plot
library(EnhancedVolcano)

EnhancedVolcano(lihc_resBlack_df,
                lab = lihc_resBlack_df$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Black Race LIHC Volcano Plot')
