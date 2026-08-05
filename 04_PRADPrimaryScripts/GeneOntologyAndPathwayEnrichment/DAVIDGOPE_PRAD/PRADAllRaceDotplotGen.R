# Load required libraries
library(readr)
library(ggplot2)
library(dplyr)
library(stringr)

# Define folders for All races
AllUpDir <- "/Users/arnavjoshi/Desktop/DAVIDGOPE_PRAD/AllUpReg"
AllDownDir <- "/Users/arnavjoshi/Desktop/DAVIDGOPE_PRAD/AllDownReg"

# Get CSV file paths
AllUpFiles <- list.files(
  AllUpDir,
  pattern = "\\.csv$",
  full.names = TRUE
)

AllDownFiles <- list.files(
  AllDownDir,
  pattern = "\\.csv$",
  full.names = TRUE
)

# Combine all files
AllFiles <- c(AllUpFiles, AllDownFiles)

# Define output PDF
AllOutputPDF <- paste0(
  "/Users/arnavjoshi/Desktop/DAVIDGOPE_PRAD/",
  "PRAD_AllRace_GO_PE_Dotplots.pdf"
)

pdf(AllOutputPDF, width = 10, height = 8)

# Loop through each file and generate a dotplot
for (AllFile in AllFiles) {
  
  # Read DAVID CSV file
  AllData <- read_csv(
    AllFile,
    show_col_types = FALSE,
    trim_ws = TRUE
  )
  
  # Convert relevant columns to numeric
  AllData <- AllData %>%
    mutate(
      Count = as.numeric(Count),
      `P-Value` = as.numeric(`P-Value`),
      `Fold Enrichment` = as.numeric(`Fold Enrichment`)
    )
  
  # Select the first 30 DAVID results
  AllTop30 <- AllData %>%
    slice_head(n = 30) %>%
    arrange(`Fold Enrichment`)
  
  # Reorder terms by fold enrichment
  AllTop30$Term <- factor(
    AllTop30$Term,
    levels = AllTop30$Term[
      order(AllTop30$`Fold Enrichment`)
    ]
  )
  
  # Determine regulation direction from filename
  AllRegulation <- case_when(
    str_detect(
      basename(AllFile),
      regex("UpReg", ignore_case = TRUE)
    ) ~ "UpRegulated",
    
    str_detect(
      basename(AllFile),
      regex("DownReg", ignore_case = TRUE)
    ) ~ "DownRegulated",
    
    TRUE ~ "UnknownRegulation"
  )
  
  # Determine DAVID category from filename
  AllCategory <- case_when(
    str_detect(
      basename(AllFile),
      regex("_BP\\.csv$", ignore_case = TRUE)
    ) ~ "BPDirect",
    
    str_detect(
      basename(AllFile),
      regex("_CC\\.csv$", ignore_case = TRUE)
    ) ~ "CCDirect",
    
    str_detect(
      basename(AllFile),
      regex("_MF\\.csv$", ignore_case = TRUE)
    ) ~ "MFDirect",
    
    str_detect(
      basename(AllFile),
      regex("_KEGG\\.csv$", ignore_case = TRUE)
    ) ~ "KEGG",
    
    TRUE ~ "UnknownCategory"
  )
  
  # Create plot title
  # Example: PRAD_AllRace_UpRegulated_CCDirect_Dotplot
  AllPlotTitle <- paste(
    "PRAD",
    "AllRace",
    AllRegulation,
    AllCategory,
    "Dotplot",
    sep = "_"
  )
  
  # Create dotplot
  AllPlot <- ggplot(
    AllTop30,
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
      title = AllPlotTitle,
      x = "Fold Enrichment",
      y = "Term"
    ) +
    theme(
      axis.text.y = element_text(size = 8),
      legend.position = "right"
    )
  
  # Add plot to PDF
  print(AllPlot)
}

# Close PDF device
dev.off()

cat("PDF saved to:", AllOutputPDF, "\n")