# Load required libraries
library(readr)
library(ggplot2)
library(dplyr)
library(stringr)

# Define folders for All races
AllUpDir <- "/Users/arnavjoshi/Desktop/BRCA_GOandPE/AllUpReg"
AllDownDir <- "/Users/arnavjoshi/Desktop/BRCA_GOandPE/AllDownReg"

# Get file paths
AllUpFiles <- list.files(AllUpDir, pattern = "\\.txt$", full.names = TRUE)
AllDownFiles <- list.files(AllDownDir, pattern = "\\.txt$", full.names = TRUE)

# Combine all files
AllAllFiles <- c(AllUpFiles, AllDownFiles)

# Output PDF path
AllOutputPDF <- "/Users/arnavjoshi/Desktop/BRCA_GOandPE/BRCA_AllRace_GO_BP_Dotplots.pdf"
pdf(AllOutputPDF, width = 10, height = 8)

# Loop through each file and generate a dotplot
for (AllFile in AllAllFiles) {
  # Read the enrichment data
  AllData <- read_delim(AllFile, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
  
  # Take top 30 rows and sort by Fold Enrichment
  AllTop30 <- AllData %>%
    slice(1:30) %>%
    arrange(`Fold Enrichment`)
  
  # Reorder terms for plotting
  AllTop30$Term <- factor(AllTop30$Term, levels = AllTop30$Term[order(AllTop30$`Fold Enrichment`)])
  
  # Extract a clean plot title from the filename
  AllPlotTitle <- paste0("BRCA_", tools::file_path_sans_ext(basename(AllFile)))
  
  # Generate dotplot
  AllPlot <- ggplot(AllTop30, aes(x = `Fold Enrichment`, y = Term)) +
    geom_point(aes(size = Count, color = PValue), alpha = 0.7) +
    scale_size_continuous(name = "Count", range = c(2, 10)) +
    scale_color_gradient(name = "P-Value", low = "red", high = "blue", trans = "log") +
    theme_minimal() +
    labs(title = AllPlotTitle,
         x = "Fold Enrichment",
         y = "Term") +
    theme(
      axis.text.y = element_text(size = 8),
      legend.position = "right"
    )
  
  # Add the plot to the PDF
  print(AllPlot)
}

# Close PDF device
dev.off()
