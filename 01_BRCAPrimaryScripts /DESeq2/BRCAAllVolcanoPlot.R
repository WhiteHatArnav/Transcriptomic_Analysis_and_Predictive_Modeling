# library(EnhancedVolcano)
# 
# EnhancedVolcano(brca_resAll,
#                 lab = rep('', nrow(brca_resAll)),
#                 x = 'log2FoldChange',
#                 y = 'pvalue',
#                 title = 'All Races BRCA Volcano Plot')

# Convert DESeqResults to a regular data frame
brca_resAll_df <- as.data.frame(brca_resAll)

# Strip version numbers from Ensembl IDs (from rownames)
brca_resAll_df$ensembl_id <- sub("\\..*", "", rownames(brca_resAll_df))

# Load biomaRt and get gene symbols
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = unique(brca_resAll_df$ensembl_id),
                  mart = ensembl)

# Load dplyr
library(dplyr)

# Join gene names to the results
brca_resAll_df <- brca_resAll_df %>%
  left_join(gene_map, by = c("ensembl_id" = "ensembl_gene_id"))

# Replace missing gene names with Ensembl IDs
brca_resAll_df$hgnc_symbol[is.na(brca_resAll_df$hgnc_symbol)] <- brca_resAll_df$ensembl_id[is.na(brca_resAll_df$hgnc_symbol)]

# Load EnhancedVolcano and plot
library(EnhancedVolcano)

EnhancedVolcano(brca_resAll_df,
                lab = brca_resAll_df$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'All Races BRCA Volcano Plot')
