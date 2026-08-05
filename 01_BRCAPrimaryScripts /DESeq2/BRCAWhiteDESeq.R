library(DESeq2)

cts <- brca_white_expression_filtered[,-1]
rownames(cts) <- brca_white_expression_filtered$Geneid
ctsmatrixWhite <- as.matrix(cts)

coldataWhite <- White_BRCA_Metadata


ddsWhite <- DESeqDataSetFromMatrix(countData = ctsmatrixWhite,
                                   colData = coldataWhite,
                                   design=~0 + Condition)
ddsWhite <- DESeq(ddsWhite)
resultsNames(ddsWhite) # lists the coefficients
brca_resWhite <- results(ddsWhite)
BRCA_WhiteRace_DESeqTestResults <- as.data.frame(brca_resWhite)
write.csv(BRCA_WhiteRace_DESeqTestResults, "BRCAFinalWhiteRace_DESeqTestResults.csv")

