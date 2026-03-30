library(argparse)
library(ggplot2)
library(readr)
library(dplyr)

setup_parser <- function() {
  parser <- ArgumentParser(description = "Generate a volcano plot from differential expression data")
  
  parser$add_argument("--input", help = "Input file (csv/tsv)", required = TRUE)
  parser$add_argument("--output", help = "Output plot file", default = "volcano_plot.png")
  parser$add_argument("--pcol", help = "Column name for p-values", default = "pvalue")
  parser$add_argument("--fccol", help = "Column name for log2 fold change", default = "log2FoldChange")
  parser$add_argument("--sep", help = "Separator: ',' or '\\t'", default = ",")
  parser$add_argument("--pthresh", help = "p-value threshold", type = "double", default = 0.05)
  parser$add_argument("--fcthresh", help = "log2FC threshold", type = "double", default = 1)
  
  return(parser)
}


read_data <- function(path, sep) {
  if (sep == "\\t") {
    data <- read_tsv(path)
  } else {
    data <- read_csv(path)
  }
  return(data)
}


prepare_data <- function(data, pcol, fccol, pthresh, fcthresh) {
  
  data <- data %>%
    mutate(
      neg_log10_p = -log10(.data[[pcol]]),
      significant = case_when(
        .data[[pcol]] < pthresh & .data[[fccol]] > fcthresh  ~ "Up",
        .data[[pcol]] < pthresh & .data[[fccol]] < -fcthresh ~ "Down",
        TRUE ~ "Not significant"
      )
    )
  
  return(data)
}


create_volcano_plot <- function(data, fccol, pthresh, fcthresh) {
  
  p <- ggplot(data, aes(
    x = .data[[fccol]],
    y = neg_log10_p,
    color = significant
  )) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c(
      "Up" = "red",
      "Down" = "blue",
      "Not significant" = "grey"
    )) +
    theme_minimal() +
    labs(
      title = "Volcano Plot",
      x = "Log2 Fold Change",
      y = "-log10(p-value)"
    ) +
    geom_vline(xintercept = c(-fcthresh, fcthresh), linetype = "dashed") +
    geom_hline(yintercept = -log10(pthresh), linetype = "dashed")
  
  return(p)
}


save_plot <- function(plot, output_path) {
  ggsave(output_path, plot = plot, width = 8, height = 6)
}

main <- function() {
  parser <- setup_parser()
  args <- parser$parse_args()
  
  data <- read_data(args$input, args$sep)
  data_prepared <- prepare_data(
    data,
    args$pcol,
    args$fccol,
    args$pthresh,
    args$fcthresh
  )
  
  plot <- create_volcano_plot(
    data_prepared,
    args$fccol,
    args$pthresh,
    args$fcthresh
  )
  
  save_plot(plot, args$output)
  
  cat("Volcano plot saved to:", args$output, "\n")
}

main()