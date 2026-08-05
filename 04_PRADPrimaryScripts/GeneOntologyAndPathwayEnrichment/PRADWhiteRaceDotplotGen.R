# Load required libraries
library(readr)
library(ggplot2)
library(dplyr)
library(stringr)
library(tools)

# Define folders for White race (PRAD)
WhiteUpDir   <- "/Users/arnavjoshi/Desktop/PRAD_GOandPE/WhiteUpReg"
WhiteDownDir <- "/Users/arnavjoshi/Desktop/PRAD_GOandPE/WhiteDownReg"

# Get file paths (for CSV files)
WhiteUpFiles   <- list.files(WhiteUpDir, pattern = "\\.csv$", full.names = TRUE)
WhiteDownFiles <- list.files(WhiteDownDir, pattern = "\\.csv$", full.names = TRUE)

# Combine all files
WhiteAllFiles <- c(WhiteUpFiles, WhiteDownFiles)

# Output PDF (prefixed with "PRAD_")
WhiteOutputPDF <- "PRAD_WhiteRace_GO_PE_Dotplots.pdf"
pdf(WhiteOutputPDF, width = 10, height = 8)

# Maximum allowed label length (same as your example)
max_label_length <- nchar("GO:0000978~RNA polymerase II cis−regulatory region sequence−specific DNA binding")

# Loop through each file and generate a dotplot
for (WhiteFile in WhiteAllFiles) {
  
  # Read CSV file
  WhiteData <- read_csv(WhiteFile, show_col_types = FALSE)
  
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
  
  # Extract filename for plot title and prefix with "PRAD_"
  WhitePlotTitle <- paste0("PRAD_", file_path_sans_ext(basename(WhiteFile)))
  
  # Create dotplot (with proper backticks for `P-Value`)
  WhitePlot <- ggplot(WhiteTop30, aes(x = `Fold Enrichment`, y = Term)) +
    geom_point(aes(size = Count, color = `P-Value`), alpha = 0.7) +
    scale_size_continuous(name = "Count", range = c(2, 10)) +
    scale_color_gradient(name = "P-Value", low = "red", high = "blue", trans = "log") +
    theme_minimal() +
    labs(title = WhitePlotTitle, x = "Fold Enrichment", y = "Term") +
    theme(
      axis.text.y = element_text(size = 8),
      legend.position = "right"
    )
  
  # Add plot to PDF
  print(WhitePlot)
}

# Close PDF device
dev.off()
