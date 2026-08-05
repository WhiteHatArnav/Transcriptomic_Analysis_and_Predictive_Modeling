library(DESeq2)

cts <- luad_all_expression_filtered[,-1]
rownames(cts) <- luad_all_expression_filtered$Geneid
ctsmatrixAll <- as.matrix(cts)

coldataAll <- All_LUAD_Metadata

ddsAll <- DESeqDataSetFromMatrix(countData = ctsmatrixAll,
                                 colData = coldataAll,
                                 design = ~0 + Condition)

ddsAll <- DESeq(ddsAll)
resultsNames(ddsAll)  # lists the coefficients

luad_resAll <- results(ddsAll)
LUAD_AllRace_DESeqTestResults <- as.data.frame(luad_resAll)

write.csv(LUAD_AllRace_DESeqTestResults, "LUADFinalAllRace_DESeqTestResults.csv")
