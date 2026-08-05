suppressMessages({
  library(EnhancedVolcano)
  library(biomaRt)
  library(dplyr)
})

args <- commandArgs(trailingOnly=TRUE)
res <- read.csv(args[1], row.names=1, check.names=FALSE)
res$ensembl_id <- sub("\\..*", "", rownames(res))

ens <- useEnsembl("genes","hsapiens_gene_ensembl")
map <- getBM(c("ensembl_gene_id","hgnc_symbol"),
             "ensembl_gene_id", unique(res$ensembl_id), ens)

res <- left_join(res, map, by=c("ensembl_id"="ensembl_gene_id"))
res$hgnc_symbol[is.na(res$hgnc_symbol)] <- res$ensembl_id[is.na(res$hgnc_symbol)]

pdf(args[2],8,7)
EnhancedVolcano(res, lab=res$hgnc_symbol,
                x="log2FoldChange", y="pvalue", title=args[3])
dev.off()
