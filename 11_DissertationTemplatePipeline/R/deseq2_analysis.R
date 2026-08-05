suppressMessages(library(DESeq2))

args <- commandArgs(trailingOnly=TRUE)
expr <- read.csv(args[1], row.names=1, check.names=FALSE)
if ("Sample ID" %in% rownames(expr)) {
  expr <- expr[rownames(expr)!="Sample ID", , drop=FALSE]
}

meta <- read.csv(args[2], check.names=FALSE)
rownames(meta) <- meta[["Sample ID"]]

expr <- expr[, rownames(meta), drop=FALSE]
cts <- as.matrix(expr)
storage.mode(cts) <- "integer"

dds <- DESeqDataSetFromMatrix(cts, meta, ~Condition)
dds <- DESeq(dds)
write.csv(as.data.frame(results(dds)), args[3])
