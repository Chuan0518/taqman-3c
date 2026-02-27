test_that("sequence validation works", {
  expect_equal(taqman3c:::validate_sequence("acgt"), "ACGT")
  expect_error(taqman3c:::validate_sequence("ACGN"), "invalid bases")
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

test_that("3C assay API returns data.frame", {
  seq <- paste0(
    "ATGCGTACGTTACGCGTATCGATCGATGCTAGCTAGCTAGCGTATCGATCGATCGTACGATCGTAGCTAGCTAGC",
    "GTACGTAGCTAGCATCGATCGTACGATCGATCGTAGCTAGCATCGATCGATCGTAGCTAGCTAGCATCGATCGTA",
    "ATGCGTACGTTACGCGTATCGATCGATGCTAGCTAGCTAGCGTATCGATCGATCGTACGATCGTAGCTAGCTAGC",
    "GTACGTAGCTAGCATCGATCGTACGATCGATCGTAGCTAGCATCGATCGATCGTAGCTAGCTAGCATCGATCGTA"
  )

  test_regions <- data.frame(
    region_id = c("enh1", "enh2"),
    start = c(250, 500),
    end = c(430, 700),
    cut_site = c(390, 640),
    stringsAsFactors = FALSE
  )

  out <- design_3c_taqman_assays(
    sequence = seq,
    anchor_region = c(40, 220),
    anchor_cut_site = 180,
    test_regions = test_regions,
    distance_to_cut = c(20, 120),
    amplicon_size_max = 400
  )

  expect_true(is.data.frame(out))
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
