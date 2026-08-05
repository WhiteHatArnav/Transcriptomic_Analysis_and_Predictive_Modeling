# Load required libraries
library(readr)
library(ggplot2)
library(dplyr)
library(stringr)

# Define folders for Black race
BlackUpDir <- "/Users/arnavjoshi/Desktop/BRCA_GOandPE/BlackUpReg"
BlackDownDir <- "/Users/arnavjoshi/Desktop/BRCA_GOandPE/BlackDownReg"

# Get file paths
BlackUpFiles <- list.files(BlackUpDir, pattern = "\\.txt$", full.names = TRUE)
BlackDownFiles <- list.files(BlackDownDir, pattern = "\\.txt$", full.names = TRUE)

# Combine all files
BlackAllFiles <- c(BlackUpFiles, BlackDownFiles)

# Output PDF path
BlackOutputPDF <- "/Users/arnavjoshi/Desktop/BRCA_GOandPE/BRCA_BlackRace_GO_PE_Dotplots.pdf"
pdf(BlackOutputPDF, width = 10, height = 8)

# Loop through each file and generate a dotplot
for (BlackFile in BlackAllFiles) {
  # Read the enrichment data
  BlackData <- read_delim(BlackFile, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
  
  # Take top 30 rows and sort by Fold Enrichment
  BlackTop30 <- BlackData %>%
    slice(1:30) %>%
    arrange(`Fold Enrichment`)
  
  # Reorder terms for plotting
  BlackTop30$Term <- factor(BlackTop30$Term, levels = BlackTop30$Term[order(BlackTop30$`Fold Enrichment`)])
  
  # Extract a clean plot title from the filename
  BlackPlotTitle <- paste0("BRCA_", tools::file_path_sans_ext(basename(BlackFile)))
  
  # Generate dotplot
  BlackPlot <- ggplot(BlackTop30, aes(x = `Fold Enrichment`, y = Term)) +
    geom_point(aes(size = Count, color = PValue), alpha = 0.7) +
    scale_size_continuous(name = "Count", range = c(2, 10)) +
    scale_color_gradient(name = "P-Value", low = "red", high = "blue", trans = "log") +
    theme_minimal() +
    labs(title = BlackPlotTitle,
         x = "Fold Enrichment",
         y = "Term") +
    theme(
      axis.text.y = element_text(size = 8),
      legend.position = "right"
    )
  
  # Add the plot to the PDF
  print(BlackPlot)
}

# Close PDF device
dev.off()
