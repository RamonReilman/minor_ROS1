library(argparse)
library(biomaRt)
library(dplyr)
setup_parser <- function() {
  parser <- ArgumentParser(description = "Add count data together and annotate gene names")
  parser$add_argument("--input", help = "Folder containing the star read counts.", required = T)
  parser$add_argument("--output", help = "Folder where the count data will be outputted", default = "./output.csv")

  return(parser)
}


read_sample <- function(path, file.name) {
  sample.name <- strsplit(file.name, ".", fixed = TRUE)[[1]][1]
  sample <- read.table(paste(path, file.name, sep = "/"), header = FALSE, sep="\t", 
                       row.names = NULL, skip = 4)
  names(sample)[2] <- sample.name

  return(sample[c(1, 2)])
}

merge_counts <- function(path) {
  file.names <- list.files(path, pattern = "*ReadsPerGene.out.tab")
  counts <- read_sample(path, file.names[1])
  
  for (file.name in file.names[2:length(file.names)]) {
  sample <- read_sample(path, file.name)
  counts <- merge(counts, sample, by = 1)
  }
  rownames(counts) <- counts$V1
  counts <- counts[-1]
  return(data.frame(counts))
  
}

annotate_genes <- function(counts){
  ensembl_ids <- gsub("\\..*$", "", rownames(counts))

  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

  gene_info <- getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol", "external_gene_name"),
      filters = "ensembl_gene_id",
      values = ensembl_ids,
      mart = ensembl
  )
  return(gene_info)
}

flatten_gene_info <- function(gene_info) {
  gene_info_unique <- gene_info %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    hgnc_symbol = paste(unique(hgnc_symbol[!is.na(hgnc_symbol) & hgnc_symbol != ""]), 
                       collapse = "/"),
    external_gene_name = first(na.omit(external_gene_name))
  ) %>%
  ungroup()
  
  return(gene_info_unique)
}

merge_counts_with_symbols <- function(counts, gene_info_unique) {
  ensembl_ids <- gsub("\\..*$", "", rownames(counts))
  
  counts_df <- data.frame(
    ensembl_gene_id = ensembl_ids,
    counts,
    row.names = NULL
  )
  
  merged_df <- merge(
    counts_df,
    gene_info_unique[, c("ensembl_gene_id", "hgnc_symbol")],
    by = "ensembl_gene_id",
    all.x = TRUE
  )
  
  merged_df$hgnc_symbol[is.na(merged_df$hgnc_symbol)] <- 
    merged_df$ensembl_gene_id[is.na(merged_df$hgnc_symbol)]
  
  merged_df$hgnc_symbol[merged_df$hgnc_symbol == ""] <- 
    merged_df$ensembl_gene_id[merged_df$hgnc_symbol == ""]
  
  rownames(merged_df) <- make.unique(merged_df$hgnc_symbol, sep = "_")
  
  final_counts <- merged_df[, !names(merged_df) %in% c("ensembl_gene_id", "hgnc_symbol")]
  
  return(final_counts)
}



main <- function() {
  parser <- setup_parser()
  args <- parser$parse_args()
  count_path <- args$input
  counts <- merge_counts(count_path)
  gene_info <- annotate_genes(counts)
  gene_info_unique <- flatten_gene_info(gene_info)

  counts_annotated <- merge_counts_with_symbols(counts, gene_info_unique)


  output_file <- file.path(args$output)
  write.csv(counts_annotated, output_file)
}

main()