# Load required libraries
library(readr)
library(ggplot2)
library(dplyr)
library(stringr)

# Define folders for All race
AllUpDir <- "/Users/arnavjoshi/Desktop/LIHC_GOandPE/AllUpReg"
AllDownDir <- "/Users/arnavjoshi/Desktop/LIHC_GOandPE/AllDownReg"

# Get file paths
AllUpFiles <- list.files(AllUpDir, pattern = "\\.txt$", full.names = TRUE)
AllDownFiles <- list.files(AllDownDir, pattern = "\\.txt$", full.names = TRUE)

# Combine all files
AllAllFiles <- c(AllUpFiles, AllDownFiles)

# Output PDF (prefixed with "LIHC_")
AllOutputPDF <- "LIHC_AllRace_GO_BP_Dotplots.pdf"
pdf(AllOutputPDF, width = 10, height = 8)

# Maximum allowed label length (based on your example)
max_label_length <- nchar("GO:0000978~RNA polymerase II cis−regulatory region sequence−specific DNA binding")

# Loop through each file and generate a dotplot
for (AllFile in AllAllFiles) {
  # Read file
  AllData <- read_delim(AllFile, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
  
  # Take top 30 rows and sort by Fold Enrichment
  AllTop30 <- AllData %>%
    slice(1:30) %>%
    arrange(`Fold Enrichment`)
  
  # Truncate Term labels longer than max_label_length and add "..."
  AllTop30$Term <- sapply(AllTop30$Term, function(x) {
    if (nchar(x) > max_label_length) {
      paste0(substr(x, 1, max_label_length), "...")
    } else {
      x
    }
  })
  
  # Reorder Term factor for plotting
  AllTop30$Term <- factor(AllTop30$Term, levels = AllTop30$Term[order(AllTop30$`Fold Enrichment`)])
  
  # Extract filename for plot title and prefix with "LIHC_"
  AllPlotTitle <- paste0("LIHC_", tools::file_path_sans_ext(basename(AllFile)))
  
  # Create dotplot
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
  
  # Add plot to PDF
  print(AllPlot)
}

# Close PDF device
dev.off()
