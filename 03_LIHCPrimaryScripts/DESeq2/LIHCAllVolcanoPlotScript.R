# Convert DESeqResults to a regular data frame
lihc_resAll_df <- as.data.frame(lihc_resAll)

# Strip version numbers from Ensembl IDs (from rownames)
lihc_resAll_df$ensembl_id <- sub("\\..*", "", rownames(lihc_resAll_df))

# Load biomaRt and get gene symbols
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = unique(lihc_resAll_df$ensembl_id),
                  mart = ensembl)

# Load dplyr
library(dplyr)

# Join gene names to the results
lihc_resAll_df <- lihc_resAll_df %>%
  left_join(gene_map, by = c("ensembl_id" = "ensembl_gene_id"))

# Replace missing gene names with Ensembl IDs
lihc_resAll_df$hgnc_symbol[is.na(lihc_resAll_df$hgnc_symbol)] <- lihc_resAll_df$ensembl_id[is.na(lihc_resAll_df$hgnc_symbol)]

# Load EnhancedVolcano and plot
library(EnhancedVolcano)

EnhancedVolcano(lihc_resAll_df,
                lab = lihc_resAll_df$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'All Races LIHC Volcano Plot')

