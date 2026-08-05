library(DESeq2)

cts <- lihc_all_expression_filtered[,-1] 
rownames(cts) <- lihc_all_expression_filtered$Geneid
ctsmatrixAll <- as.matrix(cts)

coldataAll <- All_LIHC_Metadata

ddsAll <- DESeqDataSetFromMatrix(countData = ctsmatrixAll,
                                 colData = coldataAll,
                                 design = ~0 + Condition)
ddsAll <- DESeq(ddsAll)
resultsNames(ddsAll) # lists the coefficients
lihc_resAll <- results(ddsAll)
LIHC_AllRace_DESeqTestResults <- as.data.frame(lihc_resAll)
write.csv(LIHC_AllRace_DESeqTestResults, "LIHC_AllRace_DESeqTestResults.csv")
