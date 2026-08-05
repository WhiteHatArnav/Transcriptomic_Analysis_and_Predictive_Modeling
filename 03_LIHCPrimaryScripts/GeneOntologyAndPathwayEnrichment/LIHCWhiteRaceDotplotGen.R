# Load required libraries
library(readr)
library(ggplot2)
library(dplyr)
library(stringr)

# Define folders for White race
WhiteUpDir <- "/Users/arnavjoshi/Desktop/LIHC_GOandPE/WhiteUpReg"
WhiteDownDir <- "/Users/arnavjoshi/Desktop/LIHC_GOandPE/WhiteDownReg"

# Get file paths
WhiteUpFiles <- list.files(WhiteUpDir, pattern = "\\.txt$", full.names = TRUE)
WhiteDownFiles <- list.files(WhiteDownDir, pattern = "\\.txt$", full.names = TRUE)

# Combine all files
WhiteAllFiles <- c(WhiteUpFiles, WhiteDownFiles)

# Output PDF (prefixed with "LIHC_")
WhiteOutputPDF <- "LIHC_WhiteRace_GO_PE_Dotplots.pdf"
pdf(WhiteOutputPDF, width = 10, height = 8)

# Maximum allowed label length (based on your example)
max_label_length <- nchar("GO:0000978~RNA polymerase II cis−regulatory region sequence−specific DNA binding")

# Loop through each file and generate a dotplot
for (WhiteFile in WhiteAllFiles) {
  # Read file
  WhiteData <- read_delim(WhiteFile, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
  
  # Take top 30 rows and sort by Fold Enrichment
  WhiteTop30 <- WhiteData %>%
    slice(1:30) %>%
    arrange(`Fold Enrichment`)
  
  # Truncate Term labels longer than max_label_length and add "..."
  WhiteTop30$Term <- sapply(WhiteTop30$Term, function(x) {
    if (nchar(x) > max_label_length) {
      paste0(substr(x, 1, max_label_length), "...")
    } else {
      x
    }
  })
  
  # Reorder Term factor for plotting
  WhiteTop30$Term <- factor(WhiteTop30$Term, levels = WhiteTop30$Term[order(WhiteTop30$`Fold Enrichment`)])
  
  # Extract filename for plot title and prefix with "LIHC_"
  WhitePlotTitle <- paste0("LIHC_", tools::file_path_sans_ext(basename(WhiteFile)))
  
  # Create dotplot
  WhitePlot <- ggplot(WhiteTop30, aes(x = `Fold Enrichment`, y = Term)) +
    geom_point(aes(size = Count, color = PValue), alpha = 0.7) +
    scale_size_continuous(name = "Count", range = c(2, 10)) +
    scale_color_gradient(name = "P-Value", low = "red", high = "blue", trans = "log") +
    theme_minimal() +
    labs(title = WhitePlotTitle,
         x = "Fold Enrichment",
         y = "Term") +
    theme(
      axis.text.y = element_text(size = 8),
      legend.position = "right"
    )
  
  # Add plot to PDF
  print(WhitePlot)
}

# Close PDF device
dev.off()
