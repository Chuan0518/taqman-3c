#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  cat(paste0(
    "Usage:\n",
    "  Rscript inst/scripts/design_taqman_cli.R <SEQUENCE> <ANCHOR_START> <ANCHOR_END> <ANCHOR_CUT> <TEST_REGIONS_TSV> [GENOME_FASTA] [MAX_HITS]\n\n",
    "TSV columns: region_id, start, end, cut_site\n"
  ))
  quit(status = 1)
}

source("R/design_taqman.R")

sequence <- args[[1]]
anchor_region <- c(as.integer(args[[2]]), as.integer(args[[3]]))
anchor_cut_site <- as.integer(args[[4]])
test_regions <- read.table(args[[5]], sep = "\t", header = TRUE, stringsAsFactors = FALSE)

uniqueness_genome_fasta <- ifelse(length(args) >= 6, args[[6]], NA)
if (is.na(uniqueness_genome_fasta) || uniqueness_genome_fasta == "") uniqueness_genome_fasta <- NULL
uniqueness_max_hits <- ifelse(length(args) >= 7, as.integer(args[[7]]), 1)

assays <- design_3c_taqman_assays(
  sequence = sequence,
  anchor_region = anchor_region,
  anchor_cut_site = anchor_cut_site,
  test_regions = test_regions,
  uniqueness_genome_fasta = uniqueness_genome_fasta,
  uniqueness_max_hits = uniqueness_max_hits
)
print_assays(assays)
