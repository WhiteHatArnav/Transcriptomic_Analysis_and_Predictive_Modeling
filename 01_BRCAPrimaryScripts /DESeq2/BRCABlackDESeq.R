library(DESeq2)

cts <- brca_black_expression_filtered[,-1]
rownames(cts) <- brca_black_expression_filtered$Geneid
ctsmatrixBlack <- as.matrix(cts)

coldataBlack <- Black_BRCA_Metadata


ddsBlack <- DESeqDataSetFromMatrix(countData = ctsmatrixBlack,
                                   colData = coldataBlack,
                                   design=~0 + Condition)
ddsBlack <- DESeq(ddsBlack)
resultsNames(ddsBlack) # lists the coefficients
brca_resBlack <- results(ddsBlack)
BRCA_BlackRace_DESeqTestResults <- as.data.frame(brca_resBlack)
write.csv(BRCA_BlackRace_DESeqTestResults, "BRCA_BlackRace_DESeqTestResults.csv")