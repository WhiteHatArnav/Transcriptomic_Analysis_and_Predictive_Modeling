library(readr)
library(ggplot2)
library(dplyr)
library(stringr)
library(tools)

args <- commandArgs(trailingOnly = TRUE)

csv_files <- strsplit(args[1], ",")[[1]]
output_pdf <- args[2]
cancer_code <- args[3]

pdf(output_pdf, width = 10, height = 8)

max_label_length <- nchar("GO:0000978~RNA polymerase II cis−regulatory region sequence−specific DNA binding")

for (csv_file in csv_files) {
  data <- read_csv(csv_file, show_col_types = FALSE)
  if (nrow(data) == 0) next

  top30 <- data %>%
    slice(1:30) %>%
    arrange(`Fold Enrichment`)

  top30$Term <- sapply(top30$Term, function(x) {
    if (nchar(x) > max_label_length) {
      paste0(substr(x, 1, max_label_length), "...")
    } else {
      x
    }
  })

  top30$Term <- factor(top30$Term, levels = top30$Term[order(top30$`Fold Enrichment`)])

  plot_title <- paste0(cancer_code, "_", file_path_sans_ext(basename(csv_file)))

  p <- ggplot(top30, aes(x = `Fold Enrichment`, y = Term)) +
    geom_point(aes(size = Count, color = `P-Value`), alpha = 0.7) +
    scale_size_continuous(name = "Count", range = c(2, 10)) +
    scale_color_gradient(name = "P-Value", low = "red", high = "blue", trans = "log") +
    theme_minimal() +
    labs(title = plot_title, x = "Fold Enrichment", y = "Term") +
    theme(axis.text.y = element_text(size = 8), legend.position = "right")

  print(p)
}

dev.off()
