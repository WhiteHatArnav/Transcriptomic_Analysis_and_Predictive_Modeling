library(DESeq2)

cts <- lihc_black_expression_filtered[,-1]
rownames(cts) <- lihc_black_expression_filtered$Geneid

ctsmatrixBlack <- as.matrix(cts)

coldataBlack <- Black_LIHC_Metadata

ddsBlack <- DESeqDataSetFromMatrix(countData = ctsmatrixBlack,
                                   colData = coldataBlack,
                                   design=~0 + Condition)
ddsBlack <- DESeq(ddsBlack)
resultsNames(ddsBlack) # lists the coefficients
lihc_resBlack <- results(ddsBlack)
LIHC_BlackRace_DESeqTestResults <- as.data.frame(lihc_resBlack)
write.csv(LIHC_BlackRace_DESeqTestResults, "LIHC_BlackRace_DESeqTestResults.csv")
