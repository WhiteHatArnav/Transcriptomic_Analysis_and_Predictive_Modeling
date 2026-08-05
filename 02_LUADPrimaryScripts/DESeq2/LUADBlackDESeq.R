library(DESeq2)

cts <- luad_black_expression_filtered[,-1]
rownames(cts) <- luad_black_expression_filtered$Geneid
ctsmatrixBlack <- as.matrix(cts)

coldataBlack <- Black_LUAD_Metadata

ddsBlack <- DESeqDataSetFromMatrix(countData = ctsmatrixBlack,
                                   colData = coldataBlack,
                                   design = ~0 + Condition)

ddsBlack <- DESeq(ddsBlack)
resultsNames(ddsBlack)  # lists the coefficients

luad_resBlack <- results(ddsBlack)
LUAD_BlackRace_DESeqTestResults <- as.data.frame(luad_resBlack)

write.csv(LUAD_BlackRace_DESeqTestResults, "LUADFinalBlackRace_DESeqTestResults.csv")
