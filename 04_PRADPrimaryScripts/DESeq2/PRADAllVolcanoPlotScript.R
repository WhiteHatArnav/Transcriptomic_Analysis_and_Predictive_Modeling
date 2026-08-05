# Load required libraries
library(DESeq2)
library(biomaRt)
library(EnhancedVolcano)
library(dplyr)

# Convert DESeqResults to a regular data frame
prad_resAll_df <- as.data.frame(prad_resAll)

# Strip version numbers from Ensembl IDs
prad_resAll_df$ensembl_id <- sub("\\..*", "", rownames(prad_resAll_df))

# Get gene names using biomaRt
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = unique(prad_resAll_df$ensembl_id),
  mart = ensembl
)

# Join gene names to the DESeq results
prad_resAll_df <- prad_resAll_df %>%
  left_join(gene_map, by = c("ensembl_id" = "ensembl_gene_id"))

# Fill missing gene names with Ensembl IDs
prad_resAll_df$hgnc_symbol[is.na(prad_resAll_df$hgnc_symbol)] <-
  prad_resAll_df$ensembl_id[is.na(prad_resAll_df$hgnc_symbol)]

# Plot volcano plot
EnhancedVolcano(prad_resAll_df,
                lab = prad_resAll_df$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'All Races PRAD Volcano Plot')
