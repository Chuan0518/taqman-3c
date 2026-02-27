#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat(paste0(
    "Usage:\n",
    "  Rscript inst/scripts/design_from_bedpe_cli.R <SEQUENCE> <BEDPE_PATH> [SEQUENCE_GENOMIC_START] [REGION_PADDING] [GENOME_FASTA] [MAX_HITS]\n\n",
    "BEDPE col7 anchor flag: A=first interval is anchor, B=second interval is anchor\n"
  ))
  quit(status = 1)
}

source("R/design_taqman.R")

sequence <- args[[1]]
bedpe_path <- args[[2]]
sequence_genomic_start <- ifelse(length(args) >= 3, as.integer(args[[3]]), 1)
region_padding <- ifelse(length(args) >= 4, as.integer(args[[4]]), 0)
uniqueness_genome_fasta <- ifelse(length(args) >= 5, args[[5]], NA)
if (is.na(uniqueness_genome_fasta) || uniqueness_genome_fasta == "") uniqueness_genome_fasta <- NULL
uniqueness_max_hits <- ifelse(length(args) >= 6, as.integer(args[[6]]), 1)

assays <- design_3c_from_bedpe(
  sequence = sequence,
  bedpe_path = bedpe_path,
  sequence_genomic_start = sequence_genomic_start,
  region_padding = region_padding,
  uniqueness_genome_fasta = uniqueness_genome_fasta,
  uniqueness_max_hits = uniqueness_max_hits
)

print_assays(assays)
