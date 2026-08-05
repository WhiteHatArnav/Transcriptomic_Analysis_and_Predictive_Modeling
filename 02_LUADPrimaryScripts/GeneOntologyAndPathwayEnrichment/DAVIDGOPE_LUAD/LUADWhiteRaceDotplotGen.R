# Load required libraries
library(readr)
library(ggplot2)
library(dplyr)
library(stringr)

# Define folders for White race
WhiteUpDir <- "/Users/arnavjoshi/Desktop/DAVIDGOPE_LUAD/WhiteUpReg"
WhiteDownDir <- "/Users/arnavjoshi/Desktop/DAVIDGOPE_LUAD/WhiteDownReg"

# Get CSV file paths
WhiteUpFiles <- list.files(
  WhiteUpDir,
  pattern = "\\.csv$",
  full.names = TRUE
)

WhiteDownFiles <- list.files(
  WhiteDownDir,
  pattern = "\\.csv$",
  full.names = TRUE
)

# Combine all files
WhiteAllFiles <- c(WhiteUpFiles, WhiteDownFiles)

# Output PDF
WhiteOutputPDF <- paste0(
  "/Users/arnavjoshi/Desktop/DAVIDGOPE_LUAD/",
  "LUAD_WhiteRace_GO_PE_Dotplots.pdf"
)

pdf(WhiteOutputPDF, width = 10, height = 8)

# Loop through each file and generate a dotplot
for (WhiteFile in WhiteAllFiles) {
  
  # Read DAVID CSV file
  WhiteData <- read_csv(
    WhiteFile,
    show_col_types = FALSE,
    trim_ws = TRUE
  )
  
  # Convert relevant columns to numeric
  WhiteData <- WhiteData %>%
    mutate(
      Count = as.numeric(Count),
      `P-Value` = as.numeric(`P-Value`),
      `Fold Enrichment` = as.numeric(`Fold Enrichment`)
    )
  
  # Select the first 30 DAVID results and sort by fold enrichment
  WhiteTop30 <- WhiteData %>%
    slice_head(n = 30) %>%
    arrange(`Fold Enrichment`)
  
  # Reorder terms according to fold enrichment
  WhiteTop30$Term <- factor(
    WhiteTop30$Term,
    levels = WhiteTop30$Term[
      order(WhiteTop30$`Fold Enrichment`)
    ]
  )
  
  # Determine regulation direction from filename
  WhiteRegulation <- case_when(
    str_detect(basename(WhiteFile), regex("UpReg", ignore_case = TRUE)) ~
      "UpRegulated",
    str_detect(basename(WhiteFile), regex("DownReg", ignore_case = TRUE)) ~
      "DownRegulated",
    TRUE ~ "UnknownRegulation"
  )
  
  # Determine DAVID category from filename
  WhiteCategory <- case_when(
    str_detect(basename(WhiteFile), regex("_BP\\.csv$", ignore_case = TRUE)) ~
      "BPDirect",
    str_detect(basename(WhiteFile), regex("_CC\\.csv$", ignore_case = TRUE)) ~
      "CCDirect",
    str_detect(basename(WhiteFile), regex("_MF\\.csv$", ignore_case = TRUE)) ~
      "MFDirect",
    str_detect(basename(WhiteFile), regex("_KEGG\\.csv$", ignore_case = TRUE)) ~
      "KEGG",
    TRUE ~ "UnknownCategory"
  )
  
  # Create title in the requested format
  # Example: LUAD_WhiteRace_UpRegulated_CCDirect_Dotplot
  WhitePlotTitle <- paste(
    "LUAD",
    "WhiteRace",
    WhiteRegulation,
    WhiteCategory,
    "Dotplot",
    sep = "_"
  )
  
  # Create dotplot
  WhitePlot <- ggplot(
    WhiteTop30,
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
      title = WhitePlotTitle,
      x = "Fold Enrichment",
      y = "Term"
    ) +
    theme(
      axis.text.y = element_text(size = 8),
      legend.position = "right"
    )
  
  # Add plot to PDF
  print(WhitePlot)
}

# Close PDF device
dev.off()

cat("PDF saved to:", WhiteOutputPDF, "\n")