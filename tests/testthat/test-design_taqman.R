test_that("sequence validation works", {
  expect_equal(taqman3c:::validate_sequence("acgt"), "ACGT")
  expect_error(taqman3c:::validate_sequence("ACGN"), "非法碱基")
})

test_that("reverse complement works", {
  expect_equal(taqman3c:::revcomp("ATGC"), "GCAT")
})

test_that("local FASTA uniqueness counting works", {
  f <- tempfile(fileext = ".fa")
  writeLines(c(
    ">chr1",
    "ACGTACGTACGT",
    ">chr2",
    "TTTTACGTAAAA"
  ), con = f)

  seqs <- taqman3c:::read_fasta_sequences(f)
  expect_equal(length(seqs), 2)
  expect_equal(taqman3c:::count_exact_hits("ACGT", seqs), 4)
  expect_equal(taqman3c:::count_exact_hits("CCCC", seqs), 0)
})

test_that("BEDPE anchor flag parsing works", {
  bed <- tempfile(fileext = ".bedpe")
  writeLines(c(
    "chr1	5000	6000	chr1	10000	11000	A",
    "chr1	7000	8000	chr1	12000	13000	A"
  ), con = bed)

  x <- taqman3c:::read_bedpe_with_anchor(bed)
  expect_equal(nrow(x), 2)
  expect_equal(x$anchor_flag[1], "A")

  parsed <- taqman3c:::parse_bedpe_to_3c_regions(
    bedpe = x,
    sequence_n = 20000,
    sequence_genomic_start = 1,
    region_padding = 0
  )

  expect_equal(parsed$anchor_region[1], 5001)
  expect_equal(parsed$anchor_region[2], 6000)
  expect_equal(parsed$test_regions$start[1], 10001)
  expect_equal(parsed$test_regions$end[1], 11000)
})


test_that("direct sequence input can include trailing newline", {
  expect_equal(taqman3c:::validate_sequence("GATC\n"), "GATC")
})


test_that("bedpe wrapper defaults genome.fa uniqueness path", {
  f <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t0\t10\tchr1\t20\t30\tA", con = f)
  # only verify signature-level behavior by ensuring missing genome.fa path triggers error
  expect_error(
    design_3c_from_bedpe(sequence = "GATCGATCGATCGATCGATCGATCGATC", bedpe_path = f),
    "uniqueness_genome_fasta 文件不存在"
  )
})
