library(DESeq2)

# ---- Counts matrix from PRAD White cohort ----
cts <- prad_white_expression_filtered[, -1]
rownames(cts) <- prad_white_expression_filtered$Geneid
ctsmatrixWhite <- as.matrix(cts)

# ---- Metadata (must contain columns: Sample, Condition) ----
coldataWhite <- White_PRAD_Metadata
rownames(coldataWhite) <- coldataWhite$Sample
# Reorder metadata to match columns in counts
coldataWhite <- coldataWhite[colnames(ctsmatrixWhite), , drop = FALSE]

# ---- DESeq2 dataset, design with no intercept (Tumor vs Normal encoded via Condition) ----
ddsWhite <- DESeqDataSetFromMatrix(
  countData = ctsmatrixWhite,
  colData   = coldataWhite,
  design    = ~ 0 + Condition
)

ddsWhite <- DESeq(ddsWhite)
resultsNames(ddsWhite)  # lists available coefficients

# Default results (if you want a specific contrast, set contrast= c("Condition","Tumor","Normal"))
prad_resWhite <- results(ddsWhite)
PRAD_WhiteRace_DESeqTestResults <- as.data.frame(prad_resWhite)

write.csv(PRAD_WhiteRace_DESeqTestResults, "PRADFinalWhiteRace_DESeqTestResults.csv", row.names = TRUE)
