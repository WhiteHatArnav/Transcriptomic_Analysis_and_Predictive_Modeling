# Load required libraries
library(readr)
library(ggplot2)
library(dplyr)
library(stringr)

# Define folders for Black race
BlackUpDir <- "/Users/arnavjoshi/Desktop/DAVIDGOPE_PRAD/BlackUpReg"
BlackDownDir <- "/Users/arnavjoshi/Desktop/DAVIDGOPE_PRAD/BlackDownReg"

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

# Define output PDF
BlackOutputPDF <- paste0(
  "/Users/arnavjoshi/Desktop/DAVIDGOPE_PRAD/",
  "PRAD_BlackRace_GO_PE_Dotplots.pdf"
)

# Use a wider PDF to accommodate long GO term labels
pdf(
  BlackOutputPDF,
  width = 14,
  height = 10
)

# Loop through each file and generate a dotplot
for (BlackFile in BlackAllFiles) {
  
  # Read DAVID CSV file
  BlackData <- read_csv(
    BlackFile,
    show_col_types = FALSE,
    trim_ws = TRUE
  )
  
  # Convert relevant columns to numeric and remove unusable rows
  BlackData <- BlackData %>%
    mutate(
      Count = as.numeric(Count),
      `P-Value` = as.numeric(`P-Value`),
      `Fold Enrichment` = as.numeric(`Fold Enrichment`)
    ) %>%
    filter(
      !is.na(Term),
      !is.na(Count),
      !is.na(`P-Value`),
      !is.na(`Fold Enrichment`),
      `P-Value` > 0
    )
  
  # Select the first 30 DAVID results and sort by fold enrichment
  BlackTop30 <- BlackData %>%
    slice_head(n = 30) %>%
    arrange(`Fold Enrichment`)
  
  # Skip empty files safely
  if (nrow(BlackTop30) == 0) {
    warning(
      paste(
        "No usable rows found in:",
        basename(BlackFile)
      )
    )
    next
  }
  
  # Wrap long GO terms to prevent compression of the plotting area
  BlackTop30 <- BlackTop30 %>%
    mutate(
      WrappedTerm = str_wrap(
        as.character(Term),
        width = 55
      )
    )
  
  # Preserve fold-enrichment ordering after wrapping
  BlackTop30$WrappedTerm <- factor(
    BlackTop30$WrappedTerm,
    levels = BlackTop30$WrappedTerm[
      order(BlackTop30$`Fold Enrichment`)
    ]
  )
  
  # Determine regulation direction from filename
  BlackRegulation <- case_when(
    str_detect(
      basename(BlackFile),
      regex("UpReg", ignore_case = TRUE)
    ) ~ "UpRegulated",
    
    str_detect(
      basename(BlackFile),
      regex("DownReg", ignore_case = TRUE)
    ) ~ "DownRegulated",
    
    TRUE ~ "UnknownRegulation"
  )
  
  # Determine DAVID category from filename
  BlackCategory <- case_when(
    str_detect(
      basename(BlackFile),
      regex("_BP\\.csv$", ignore_case = TRUE)
    ) ~ "BPDirect",
    
    str_detect(
      basename(BlackFile),
      regex("_CC\\.csv$", ignore_case = TRUE)
    ) ~ "CCDirect",
    
    str_detect(
      basename(BlackFile),
      regex("_MF\\.csv$", ignore_case = TRUE)
    ) ~ "MFDirect",
    
    str_detect(
      basename(BlackFile),
      regex("_KEGG\\.csv$", ignore_case = TRUE)
    ) ~ "KEGG",
    
    TRUE ~ "UnknownCategory"
  )
  
  # Create plot title
  # Example: PRAD_BlackRace_DownRegulated_MFDirect_Dotplot
  BlackPlotTitle <- paste(
    "PRAD",
    "BlackRace",
    BlackRegulation,
    BlackCategory,
    "Dotplot",
    sep = "_"
  )
  
  # Create dotplot
  BlackPlot <- ggplot(
    BlackTop30,
    aes(
      x = `Fold Enrichment`,
      y = WrappedTerm
    )
  ) +
    geom_point(
      aes(
        size = Count,
        color = `P-Value`
      ),
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
      axis.text.y = element_text(
        size = 7,
        lineheight = 0.9
      ),
      axis.text.x = element_text(
        size = 9
      ),
      axis.title.x = element_text(
        size = 12
      ),
      axis.title.y = element_text(
        size = 12
      ),
      plot.title = element_text(
        size = 14,
        hjust = 0.5
      ),
      legend.position = "right",
      plot.margin = margin(
        t = 10,
        r = 30,
        b = 10,
        l = 10
      )
    )
  
  # Add plot to PDF
  print(BlackPlot)
}

# Close PDF device
dev.off()

cat(
  "PDF saved to:\n",
  BlackOutputPDF,
  "\n"
)