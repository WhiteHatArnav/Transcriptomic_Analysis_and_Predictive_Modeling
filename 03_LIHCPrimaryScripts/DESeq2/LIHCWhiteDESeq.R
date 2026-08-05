library(DESeq2)

cts <- lihc_white_expression_filtered[,-1]
rownames(cts) <- lihc_white_expression_filtered$Geneid
ctsmatrixWhite <- as.matrix(cts)

coldataWhite <- White_LIHC_Metadata

ddsWhite <- DESeqDataSetFromMatrix(countData = ctsmatrixWhite,
                                   colData = coldataWhite,
                                   design=~0 + Condition)
ddsWhite <- DESeq(ddsWhite)
resultsNames(ddsWhite) # lists the coefficients
lihc_resWhite <- results(ddsWhite)
LIHC_WhiteRace_DESeqTestResults <- as.data.frame(lihc_resWhite)
write.csv(LIHC_WhiteRace_DESeqTestResults, "LIHCFinalWhiteRace_DESeqTestResults.csv")
