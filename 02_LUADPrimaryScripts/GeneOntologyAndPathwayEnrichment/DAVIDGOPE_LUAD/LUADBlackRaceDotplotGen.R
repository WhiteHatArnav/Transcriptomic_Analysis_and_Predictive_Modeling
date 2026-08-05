# Load required libraries
library(readr)
library(ggplot2)
library(dplyr)
library(stringr)

# Define folders for Black race
BlackUpDir <- "/Users/arnavjoshi/Desktop/DAVIDGOPE_LUAD/BlackUpReg"
BlackDownDir <- "/Users/arnavjoshi/Desktop/DAVIDGOPE_LUAD/BlackDownReg"

# Get CSV file paths
BlackUpFiles <- list.files(
  BlackUpDir,
  pattern = "\\.csv$",
  full.names = TRUE
)

BlackDownFiles <- list.files(
  BlackDownDir,
  pattern = "\\.csv$",
  full.names = TRUE
)

# Combine all files
BlackAllFiles <- c(BlackUpFiles, BlackDownFiles)

# Output PDF
BlackOutputPDF <- paste0(
  "/Users/arnavjoshi/Desktop/DAVIDGOPE_LUAD/",
  "LUAD_BlackRace_GO_PE_Dotplots.pdf"
)

pdf(BlackOutputPDF, width = 10, height = 8)

# Loop through each file and generate a dotplot
for (BlackFile in BlackAllFiles) {
  
  # Read DAVID CSV file
  BlackData <- read_csv(
    BlackFile,
    show_col_types = FALSE,
    trim_ws = TRUE
  )
  
  # Convert relevant columns to numeric
  BlackData <- BlackData %>%
    mutate(
      Count = as.numeric(Count),
      `P-Value` = as.numeric(`P-Value`),
      `Fold Enrichment` = as.numeric(`Fold Enrichment`)
    )
  
  # Select the first 30 DAVID results and sort by fold enrichment
  BlackTop30 <- BlackData %>%
    slice_head(n = 30) %>%
    arrange(`Fold Enrichment`)
  
  # Reorder terms according to fold enrichment
  BlackTop30$Term <- factor(
    BlackTop30$Term,
    levels = BlackTop30$Term[
      order(BlackTop30$`Fold Enrichment`)
    ]
  )
  
  # Determine regulation direction from filename
  BlackRegulation <- case_when(
    str_detect(basename(BlackFile), regex("UpReg", ignore_case = TRUE)) ~
      "UpRegulated",
    str_detect(basename(BlackFile), regex("DownReg", ignore_case = TRUE)) ~
      "DownRegulated",
    TRUE ~ "UnknownRegulation"
  )
  
  # Determine DAVID category from filename
  BlackCategory <- case_when(
    str_detect(basename(BlackFile), regex("_BP\\.csv$", ignore_case = TRUE)) ~
      "BPDirect",
    str_detect(basename(BlackFile), regex("_CC\\.csv$", ignore_case = TRUE)) ~
      "CCDirect",
    str_detect(basename(BlackFile), regex("_MF\\.csv$", ignore_case = TRUE)) ~
      "MFDirect",
    str_detect(basename(BlackFile), regex("_KEGG\\.csv$", ignore_case = TRUE)) ~
      "KEGG",
    TRUE ~ "UnknownCategory"
  )
  
  # Create title
  # Example: LUAD_BlackRace_UpRegulated_CCDirect_Dotplot
  BlackPlotTitle <- paste(
    "LUAD",
    "BlackRace",
    BlackRegulation,
    BlackCategory,
    "Dotplot",
    sep = "_"
  )
  
  # Create dotplot
  BlackPlot <- ggplot(
    BlackTop30,
    aes(x = `Fold Enrichment`, y = Term)
  ) +
    geom_point(
      aes(size = Count, color = `P-Value`),
      alpha = 0.7
    ) +
    scale_size_continuous(
      name = "Count",
      range = c(2, 10)
    ) +
    scale_color_gradient(
      name = "P-Value",
      low = "red",
      high = "blue",
      trans = "log10"
    ) +
    theme_minimal() +
    labs(
      title = BlackPlotTitle,
      x = "Fold Enrichment",
      y = "Term"
    ) +
    theme(
      axis.text.y = element_text(size = 8),
      legend.position = "right"
    )
  
  # Add plot to PDF
  print(BlackPlot)
}

# Close PDF device
dev.off()

cat("PDF saved to:", BlackOutputPDF, "\n")