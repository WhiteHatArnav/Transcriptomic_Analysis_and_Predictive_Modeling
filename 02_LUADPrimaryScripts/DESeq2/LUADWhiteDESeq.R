library(DESeq2)

cts <- luad_white_expression_filtered[,-1]
rownames(cts) <- luad_white_expression_filtered$Geneid
ctsmatrixWhite <- as.matrix(cts)

coldataWhite <- White_LUAD_Metadata

ddsWhite <- DESeqDataSetFromMatrix(countData = ctsmatrixWhite,
                                   colData = coldataWhite,
                                   design = ~0 + Condition)

ddsWhite <- DESeq(ddsWhite)
resultsNames(ddsWhite)  # lists the coefficients

luad_resWhite <- results(ddsWhite)
LUAD_WhiteRace_DESeqTestResults <- as.data.frame(luad_resWhite)

write.csv(LUAD_WhiteRace_DESeqTestResults, "LUADFinalWhiteRace_DESeqTestResults.csv")
