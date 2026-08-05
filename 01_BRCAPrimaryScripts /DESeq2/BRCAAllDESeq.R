library(DESeq2)

cts <- brca_all_expression_filtered[,-1] 
rownames(cts) <- brca_all_expression_filtered$Geneid
ctsmatrixAll <- as.matrix(cts)

coldataAll <- All_BRCA_Metadata


ddsAll <- DESeqDataSetFromMatrix(countData = ctsmatrixAll,
                                   colData = coldataAll,
                                   design=~0 + Condition)
ddsAll <- DESeq(ddsAll)
resultsNames(ddsAll) # lists the coefficients
brca_resAll <- results(ddsAll)
BRCA_AllRace_DESeqTestResults <- as.data.frame(brca_resAll)
write.csv(BRCA_AllRace_DESeqTestResults, "BRCA_AllRace_DESeqTestResults.csv")