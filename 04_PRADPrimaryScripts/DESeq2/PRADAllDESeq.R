library(DESeq2)

# ---- Counts matrix from PRAD All cohort ----
cts <- prad_all_expression_filtered[, -1]
rownames(cts) <- prad_all_expression_filtered$Geneid
ctsmatrixAll <- as.matrix(cts)

# ---- Metadata (must contain columns: Sample, Condition) ----
coldataAll <- All_PRAD_Metadata
rownames(coldataAll) <- coldataAll$Sample
# Reorder metadata to match columns in counts
coldataAll <- coldataAll[colnames(ctsmatrixAll), , drop = FALSE]

# ---- DESeq2 dataset, design with no intercept (Tumor vs Normal encoded via Condition) ----
ddsAll <- DESeqDataSetFromMatrix(
  countData = ctsmatrixAll,
  colData   = coldataAll,
  design    = ~ 0 + Condition
)

ddsAll <- DESeq(ddsAll)
resultsNames(ddsAll)  # lists available coefficients

# Default results (if you want a specific contrast, set contrast = c("Condition","Tumor","Normal"))
prad_resAll <- results(ddsAll)
PRAD_AllRace_DESeqTestResults <- as.data.frame(prad_resAll)

write.csv(PRAD_AllRace_DESeqTestResults, "PRADFinalAllRace_DESeqTestResults.csv", row.names = TRUE)
