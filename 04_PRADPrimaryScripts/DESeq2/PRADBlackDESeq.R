library(DESeq2)

# ---- Counts matrix from PRAD Black cohort ----
cts <- prad_black_expression_filtered[, -1]
rownames(cts) <- prad_black_expression_filtered$Geneid
ctsmatrixBlack <- as.matrix(cts)

# ---- Metadata (must contain columns: Sample, Condition) ----
coldataBlack <- Black_PRAD_Metadata
rownames(coldataBlack) <- coldataBlack$Sample
# Reorder metadata to match columns in counts
coldataBlack <- coldataBlack[colnames(ctsmatrixBlack), , drop = FALSE]

# ---- DESeq2 dataset, design with no intercept (Tumor vs Normal encoded via Condition) ----
ddsBlack <- DESeqDataSetFromMatrix(
  countData = ctsmatrixBlack,
  colData   = coldataBlack,
  design    = ~ 0 + Condition
)

ddsBlack <- DESeq(ddsBlack)
resultsNames(ddsBlack)  # lists available coefficients

# Default results (if you want a specific contrast, set contrast= c("Condition","Tumor","Normal"))
prad_resBlack <- results(ddsBlack)
PRAD_BlackRace_DESeqTestResults <- as.data.frame(prad_resBlack)

write.csv(PRAD_BlackRace_DESeqTestResults, "PRADFinalBlackRace_DESeqTestResults.csv", row.names = TRUE)
