#' 3C TaqMan qPCR 设计核心（内部）
#'
#' 实现面向3C的设计流程：
#' 在锚定片段使用1条引物+1条探针，并在多个测试片段上组合测试引物。
#' 所有引物在同一条基因组链窗口中筛选，
#' 探针序列输出为相对锚定引物链的反义序列。
#'
#' @param sequence DNA序列（仅允许 A/C/G/T）。
#' @param anchor_region 整数向量 c(start, end)，1-based 且闭区间。
#' @param anchor_cut_site anchor_region 内的酶切位点坐标。
#' @param test_regions 包含 `region_id`, `start`, `end`, `cut_site` 的数据框。
#' @param max_results 返回结果的最大条数。
#' @param primer_len_range 引物长度范围，默认 c(18, 30)。
#' @param probe_len_range 探针长度范围，默认 c(20, 30)。
#' @param primer_tm_range 引物Tm范围，默认 c(58, 62)。
#' @param probe_tm_delta 探针Tm高于引物Tm的范围，默认 c(8, 10)。
#' @param primer_gc_range 引物GC范围，默认 c(30, 80)。
#' @param amplicon_size_best 伪扩增子优选范围，默认 c(60, 150)。
#' @param amplicon_size_max 伪扩增子最大允许长度，默认 300。
#' @param distance_to_cut 引物到酶切位点距离范围，默认 c(50, 150)。
#' @param primer_tm_diff_max Anchor/Test引物最大Tm差，默认 3。
#' @param probe_gap_range 锚定引物与探针间距范围，默认 c(1, 10)。
#' @param secondary_enzyme_sites 不允许出现在拼接模板中的次级酶切位点序列向量
#'   （精确匹配）。
#' @param uniqueness_genome_fasta 用于全基因组唯一性过滤的本地FASTA路径。
#' @param uniqueness_max_hits 引物在正负链允许的最大精确命中次数。
#' @param enforce_probe_uniqueness 是否同时要求探针唯一性。
#' @return 按分数升序排序的候选结果 data.frame。
#'
design_3c_core <- function(sequence,
                                    anchor_region,
                                    anchor_cut_site,
                                    test_regions,
                                    max_results = 20,
                                    primer_len_range = c(18, 30),
                                    probe_len_range = c(20, 30),
                                    primer_tm_range = c(58, 62),
                                    probe_tm_delta = c(8, 10),
                                    primer_gc_range = c(30, 80),
                                    amplicon_size_best = c(60, 150),
                                    amplicon_size_max = 300,
                                    distance_to_cut = c(50, 150),
                                    primer_tm_diff_max = 3,
                                    probe_gap_range = c(1, 10),
                                    secondary_enzyme_sites = character(0),
                                    uniqueness_genome_fasta = NULL,
                                    uniqueness_max_hits = 1,
                                    enforce_probe_uniqueness = FALSE) {
  seq <- validate_sequence(sequence)
  check_region(anchor_region, nchar(seq), "anchor_region")
  check_cut_site(anchor_cut_site, anchor_region, "anchor_cut_site")
  check_test_regions(test_regions, nchar(seq))

  uniqueness_ctx <- init_uniqueness_context(uniqueness_genome_fasta, uniqueness_max_hits)

  anchor_primers <- generate_primers_near_cut(
    seq, anchor_region, anchor_cut_site,
    primer_len_range = primer_len_range,
    tm_range = primer_tm_range,
    gc_range = primer_gc_range,
    distance_to_cut = distance_to_cut
  )

  if (length(anchor_primers) == 0) return(data.frame())

  candidates <- list()
  idx <- 1

  for (a_primer in anchor_primers) {
    if (!passes_uniqueness(a_primer$sequence, uniqueness_ctx)) next

    probe_candidates <- generate_anchor_probes(
      target = seq,
      anchor_primer = a_primer,
      anchor_cut_site = anchor_cut_site,
      probe_len_range = probe_len_range,
      probe_gap_range = probe_gap_range,
      probe_tm_min = a_primer$tm + probe_tm_delta[1],
      probe_tm_max = a_primer$tm + probe_tm_delta[2],
      gc_range = primer_gc_range
    )

    if (length(probe_candidates) == 0) next

    for (i in seq_len(nrow(test_regions))) {
      tr <- test_regions[i, ]
      test_region <- c(as.integer(tr$start), as.integer(tr$end))
      test_cut <- as.integer(tr$cut_site)
      test_primers <- generate_primers_near_cut(
        seq, test_region, test_cut,
        primer_len_range = primer_len_range,
        tm_range = primer_tm_range,
        gc_range = primer_gc_range,
        distance_to_cut = distance_to_cut
      )
      if (length(test_primers) == 0) next

      for (t_primer in test_primers) {
        if (!passes_uniqueness(t_primer$sequence, uniqueness_ctx)) next
        if (abs(a_primer$tm - t_primer$tm) > primer_tm_diff_max) next

        pseudo_amplicon <- (anchor_cut_site - a_primer$start + 1) + (test_cut - t_primer$start + 1)
        if (pseudo_amplicon > amplicon_size_max) next

        chimera <- paste0(
          substr(seq, a_primer$start, anchor_cut_site - 1),
          substr(seq, t_primer$start, test_cut - 1)
        )
        if (contains_any_motif(chimera, secondary_enzyme_sites)) next

        for (probe in probe_candidates) {
          if (enforce_probe_uniqueness && !passes_uniqueness(probe$sequence, uniqueness_ctx)) next

          score <- score_3c_candidate(
            anchor_primer = a_primer,
            test_primer = t_primer,
            probe = probe,
            pseudo_amplicon = pseudo_amplicon,
            amplicon_size_best = amplicon_size_best
          )

          candidates[[idx]] <- list(
            region_id = as.character(tr$region_id),
            anchor_primer_seq = a_primer$sequence,
            anchor_primer_start = a_primer$start,
            anchor_primer_end = a_primer$end,
            anchor_primer_tm = a_primer$tm,
            anchor_primer_gc = a_primer$gc_percent,
            test_primer_seq = t_primer$sequence,
            test_primer_start = t_primer$start,
            test_primer_end = t_primer$end,
            test_primer_tm = t_primer$tm,
            test_primer_gc = t_primer$gc_percent,
            probe_seq_antisense = probe$sequence,
            probe_start = probe$start,
            probe_end = probe$end,
            probe_tm = probe$tm,
            probe_gc = probe$gc_percent,
            probe_gap_from_anchor_primer = probe$gap,
            pseudo_amplicon_size = pseudo_amplicon,
            score = score
          )
          idx <- idx + 1
        }
      }
    }
  }

  if (length(candidates) == 0) return(data.frame())

  df <- as.data.frame(do.call(rbind, lapply(candidates, as.data.frame)), stringsAsFactors = FALSE)
  num_cols <- setdiff(colnames(df), c("region_id", "anchor_primer_seq", "test_primer_seq", "probe_seq_antisense"))
  df[num_cols] <- lapply(df[num_cols], as.numeric)
  df <- df[order(df$score), , drop = FALSE]
  head(df, max_results)
}

#' 以可读格式打印候选结果
#' @param assays data.frame returned by [design_3c_from_bedpe()].
#' @export
print_assays <- function(assays) {
  if (nrow(assays) == 0) {
    cat("当前约束下未找到可用候选。\n")
    return(invisible(NULL))
  }

  for (i in seq_len(nrow(assays))) {
    row <- assays[i, ]
    cat(sprintf("候选 #%d | 区域=%s | 分数=%.2f | 伪扩增子=%d bp\n",
                i, row$region_id, row$score, row$pseudo_amplicon_size))
    cat(sprintf("  Anchor引物 (%d-%d): %s | Tm=%.1f | GC=%.1f%%\n",
                row$anchor_primer_start, row$anchor_primer_end,
                row$anchor_primer_seq, row$anchor_primer_tm, row$anchor_primer_gc))
    cat(sprintf("  Test引物   (%d-%d): %s | Tm=%.1f | GC=%.1f%%\n",
                row$test_primer_start, row$test_primer_end,
                row$test_primer_seq, row$test_primer_tm, row$test_primer_gc))
    cat(sprintf("  探针(反义) (%d-%d): %s | Tm=%.1f | GC=%.1f%% | 间距=%d\n\n",
                row$probe_start, row$probe_end,
                row$probe_seq_antisense, row$probe_tm, row$probe_gc,
                row$probe_gap_from_anchor_primer))
  }

  invisible(NULL)
}

validate_sequence <- function(seq) {
  out <- toupper(trimws(seq))
  if (nchar(out) == 0) stop("输入序列为空。")
  bad <- unique(strsplit(gsub("[ACGT]", "", out), "")[[1]])
  bad <- bad[nzchar(bad)]
  if (length(bad) > 0) {
    stop(sprintf("序列包含非法碱基: %s", paste(sort(bad), collapse = "")))
  }
  out
}

revcomp <- function(seq) {
  chars <- strsplit(toupper(seq), "")[[1]]
  comp <- c(A = "T", C = "G", G = "C", T = "A")
  paste0(rev(comp[chars]), collapse = "")
}

read_fasta_sequences <- function(fasta_path) {
  lines <- readLines(fasta_path, warn = FALSE)
  if (length(lines) == 0) stop("FASTA 文件为空")

  seqs <- character(0)
  current <- character(0)
  for (line in lines) {
    if (startsWith(line, ">")) {
      if (length(current) > 0) {
        seqs <- c(seqs, paste0(current, collapse = ""))
        current <- character(0)
      }
      next
    }
    current <- c(current, gsub("\\s+", "", line))
  }
  if (length(current) > 0) seqs <- c(seqs, paste0(current, collapse = ""))
  if (length(seqs) == 0) stop("FASTA 中未找到有效序列")

  lapply(seqs, validate_sequence)
}

count_exact_hits <- function(query, genome_sequences) {
  q <- validate_sequence(query)
  q_rc <- revcomp(q)
  total <- 0
  for (chrom in genome_sequences) {
    total <- total + count_fixed_substring(chrom, q)
    if (q_rc != q) total <- total + count_fixed_substring(chrom, q_rc)
  }
  total
}

count_fixed_substring <- function(text, pattern) {
  hits <- gregexpr(pattern, text, fixed = TRUE)[[1]]
  if (length(hits) == 1 && hits[1] == -1) return(0)
  length(hits)
}

init_uniqueness_context <- function(uniqueness_genome_fasta, uniqueness_max_hits) {
  if (is.null(uniqueness_genome_fasta)) {
    return(list(enabled = FALSE))
  }

  if (!file.exists(uniqueness_genome_fasta)) {
    stop("uniqueness_genome_fasta 文件不存在")
  }
  if (uniqueness_max_hits < 1) {
    stop("uniqueness_max_hits 必须 >= 1")
  }

  list(
    enabled = TRUE,
    genome_sequences = read_fasta_sequences(uniqueness_genome_fasta),
    max_hits = as.integer(uniqueness_max_hits),
    cache = new.env(parent = emptyenv())
  )
}

passes_uniqueness <- function(oligo_seq, uniqueness_ctx) {
  if (!isTRUE(uniqueness_ctx$enabled)) return(TRUE)

  key <- paste0("k_", oligo_seq)
  if (!exists(key, envir = uniqueness_ctx$cache, inherits = FALSE)) {
    hits <- count_exact_hits(oligo_seq, uniqueness_ctx$genome_sequences)
    assign(key, hits, envir = uniqueness_ctx$cache)
  }

  hits <- get(key, envir = uniqueness_ctx$cache, inherits = FALSE)
  hits <= uniqueness_ctx$max_hits
}

gc_percent <- function(seq) {
  chars <- strsplit(seq, "")[[1]]
  gc <- sum(chars %in% c("G", "C"))
  100 * gc / length(chars)
}

wallace_tm <- function(seq) {
  chars <- strsplit(seq, "")[[1]]
  at <- sum(chars %in% c("A", "T"))
  gc <- length(chars) - at
  2 * at + 4 * gc
}

has_homopolymer <- function(seq, max_run = 4) {
  chars <- strsplit(seq, "")[[1]]
  if (length(chars) < 2) return(FALSE)
  run <- 1
  for (i in 2:length(chars)) {
    if (chars[i] == chars[i - 1]) {
      run <- run + 1
      if (run > max_run) return(TRUE)
    } else {
      run <- 1
    }
  }
  FALSE
}

simple_self_complementarity_score <- function(seq) {
  rc <- revcomp(seq)
  best <- 0
  max_k <- min(9, nchar(seq))
  for (k in 4:max_k) {
    for (i in 1:(nchar(seq) - k + 1)) {
      sub <- substr(seq, i, i + k - 1)
      if (grepl(sub, rc, fixed = TRUE)) best <- max(best, k)
    }
  }
  best
}

tail_gc_count <- function(seq, n = 5) {
  tail <- substr(seq, max(1, nchar(seq) - n + 1), nchar(seq))
  chars <- strsplit(tail, "")[[1]]
  sum(chars %in% c("G", "C"))
}

count_base <- function(seq, base) {
  chars <- strsplit(seq, "")[[1]]
  sum(chars == base)
}

contains_any_motif <- function(seq, motifs) {
  if (length(motifs) == 0) return(FALSE)
  for (m in motifs) {
    if (nchar(m) > 0 && grepl(m, seq, fixed = TRUE)) return(TRUE)
  }
  FALSE
}

check_region <- function(region, seq_len, name) {
  if (length(region) != 2) stop(sprintf("%s 必须是 c(start, end)", name))
  if (region[1] < 1 || region[2] > seq_len || region[1] >= region[2]) {
    stop(sprintf("%s 超出序列范围", name))
  }
}

check_cut_site <- function(cut_site, region, name) {
  if (length(cut_site) != 1 || cut_site <= region[1] || cut_site > region[2]) {
    stop(sprintf("%s must lie within region and downstream of region start", name))
  }
}

check_test_regions <- function(test_regions, seq_len) {
  needed <- c("region_id", "start", "end", "cut_site")
  if (!all(needed %in% colnames(test_regions))) {
    stop("test_regions 必须包含列: region_id, start, end, cut_site")
  }
  for (i in seq_len(nrow(test_regions))) {
    r <- c(as.integer(test_regions$start[i]), as.integer(test_regions$end[i]))
    check_region(r, seq_len, sprintf("test_regions row %d", i))
    check_cut_site(as.integer(test_regions$cut_site[i]), r, sprintf("test_regions row %d cut_site", i))
  }
}

generate_primers_near_cut <- function(target,
                                      region,
                                      cut_site,
                                      primer_len_range,
                                      tm_range,
                                      gc_range,
                                      distance_to_cut) {
  out <- list()
  idx <- 1

  region_seq_start <- region[1]
  region_seq_end <- region[2]

  for (len in primer_len_range[1]:primer_len_range[2]) {
    start_min <- max(region_seq_start, cut_site - distance_to_cut[2] - len + 1)
    start_max <- min(region_seq_end - len + 1, cut_site - distance_to_cut[1] - len + 1)
    if (start_min > start_max) next

    for (start in start_min:start_max) {
      end <- start + len - 1
      if (end >= cut_site) next

      seq <- substr(target, start, end)
      tm <- wallace_tm(seq)
      gc <- gc_percent(seq)

      if (tm < tm_range[1] || tm > tm_range[2]) next
      if (gc < gc_range[1] || gc > gc_range[2]) next
      if (tail_gc_count(seq, 5) > 2) next
      if (has_homopolymer(seq, max_run = 3)) next
      if (grepl("GGGG", seq, fixed = TRUE)) next
      if (simple_self_complementarity_score(seq) >= 7) next

      out[[idx]] <- list(
        sequence = seq,
        start = start,
        end = end,
        tm = tm,
        gc_percent = gc
      )
      idx <- idx + 1
    }
  }

  out
}

generate_anchor_probes <- function(target,
                                   anchor_primer,
                                   anchor_cut_site,
                                   probe_len_range,
                                   probe_gap_range,
                                   probe_tm_min,
                                   probe_tm_max,
                                   gc_range) {
  out <- list()
  idx <- 1

  primer_end <- anchor_primer$end
  probe_start_min <- primer_end + probe_gap_range[1]
  probe_start_max <- primer_end + probe_gap_range[2]

  for (len in probe_len_range[1]:probe_len_range[2]) {
    for (start in probe_start_min:probe_start_max) {
      end <- start + len - 1
      if (end >= anchor_cut_site) next

      sense_seq <- substr(target, start, end)
      antisense_probe <- revcomp(sense_seq)
      tm <- wallace_tm(antisense_probe)
      gc <- gc_percent(antisense_probe)

      if (tm < probe_tm_min || tm > probe_tm_max) next
      if (gc < gc_range[1] || gc > gc_range[2]) next
      if (startsWith(antisense_probe, "G")) next
      if (count_base(antisense_probe, "C") <= count_base(antisense_probe, "G")) next
      if (has_homopolymer(antisense_probe, max_run = 3)) next
      if (grepl("GGGG", antisense_probe, fixed = TRUE)) next
      if (simple_self_complementarity_score(antisense_probe) >= 7) next

      out[[idx]] <- list(
        sequence = antisense_probe,
        start = start,
        end = end,
        tm = tm,
        gc_percent = gc,
        gap = start - primer_end
      )
      idx <- idx + 1
    }
  }

  out
}

score_3c_candidate <- function(anchor_primer,
                               test_primer,
                               probe,
                               pseudo_amplicon,
                               amplicon_size_best) {
  tm_delta <- abs(anchor_primer$tm - test_primer$tm)
  primer_tm_target <- 60
  primer_tm_penalty <- abs(anchor_primer$tm - primer_tm_target) + abs(test_primer$tm - primer_tm_target)

  avg_primer_tm <- (anchor_primer$tm + test_primer$tm) / 2
  probe_tm_penalty <- abs(probe$tm - (avg_primer_tm + 9))

  if (pseudo_amplicon < amplicon_size_best[1]) {
    amp_penalty <- amplicon_size_best[1] - pseudo_amplicon
  } else if (pseudo_amplicon > amplicon_size_best[2]) {
    amp_penalty <- pseudo_amplicon - amplicon_size_best[2]
  } else {
    amp_penalty <- 0
  }

  tm_delta * 2 + primer_tm_penalty + probe_tm_penalty + amp_penalty * 0.1
}



#' 从BEDPE对直接进行3C设计
#'
#' BEDPE 默认列含义：
#' chr1 start1 end1 chr2 start2 end2 anchor_flag
#' 其中 `anchor_flag` 为 `A`（第一段为Anchor）或 `B`（第二段为Anchor）。
#'
#' @param sequence DNA序列（仅允许 A/C/G/T）。
#' @param bedpe_path 本地 BEDPE 文件路径。
#' @param sequence_genomic_start 与 sequence 第1位对应的基因组1-based坐标。
#' @param region_padding BEDPE区间左右各扩展的碱基数。
#' @param uniqueness_genome_fasta 默认唯一性过滤使用的本地基因组FASTA路径。
#' @param uniqueness_max_hits 允许的最大精确命中次数。
#' @param enforce_probe_uniqueness 是否同时要求探针唯一性。
#' @param ... 透传给内部BEDPE设计核心参数。
#' @return 候选结果 data.frame。
#' @export
#'
design_3c_from_bedpe <- function(sequence,
                                 bedpe_path,
                                 sequence_genomic_start = 1,
                                 region_padding = 0,
                                 uniqueness_genome_fasta = "genome.fa",
                                 uniqueness_max_hits = 1,
                                 enforce_probe_uniqueness = FALSE,
                                 ...) {
  bedpe <- read_bedpe_with_anchor(bedpe_path)

  parsed <- parse_bedpe_to_3c_regions(
    bedpe = bedpe,
    sequence_n = nchar(validate_sequence(sequence)),
    sequence_genomic_start = sequence_genomic_start,
    region_padding = region_padding
  )

  design_3c_core(
    sequence = sequence,
    anchor_region = parsed$anchor_region,
    anchor_cut_site = parsed$anchor_cut_site,
    test_regions = parsed$test_regions,
    uniqueness_genome_fasta = uniqueness_genome_fasta,
    uniqueness_max_hits = uniqueness_max_hits,
    enforce_probe_uniqueness = enforce_probe_uniqueness,
    ...
  )
}

read_bedpe_with_anchor <- function(bedpe_path) {
  if (!file.exists(bedpe_path)) stop("bedpe_path 文件不存在")
  x <- read.table(bedpe_path, sep = "	", header = FALSE, stringsAsFactors = FALSE, comment.char = "")
  if (ncol(x) < 7) stop("BEDPE 至少需要7列: chr1 start1 end1 chr2 start2 end2 anchor_flag")
  names(x)[1:7] <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "anchor_flag")
  x$anchor_flag <- toupper(trimws(as.character(x$anchor_flag)))
  if (any(!x$anchor_flag %in% c("A", "B"))) stop("BEDPE 第7列 anchor_flag 必须是 A 或 B")
  x
}

parse_bedpe_to_3c_regions <- function(bedpe, sequence_n, sequence_genomic_start = 1, region_padding = 0) {
  if (nrow(bedpe) == 0) stop("BEDPE 文件没有数据行")
  if (region_padding < 0) stop("region_padding 必须 >= 0")

  to_local <- function(start0, end0) {
    # BEDPE: start is 0-based, end is 1-based exclusive.
    start1 <- as.integer(start0) + 1
    end1 <- as.integer(end0)
    c(start = start1 - sequence_genomic_start + 1,
      end = end1 - sequence_genomic_start + 1)
  }

  first <- bedpe[1, ]
  if (first$anchor_flag == "A") {
    anchor_g <- c(start = as.integer(first$start1) + 1, end = as.integer(first$end1))
  } else {
    anchor_g <- c(start = as.integer(first$start2) + 1, end = as.integer(first$end2))
  }
  anchor_g[1] <- anchor_g[1] - region_padding
  anchor_g[2] <- anchor_g[2] + region_padding

  anchor_region <- c(
    anchor_g[1] - sequence_genomic_start + 1,
    anchor_g[2] - sequence_genomic_start + 1
  )
  anchor_cut_site <- anchor_region[2]

  test_regions <- vector("list", nrow(bedpe))
  for (i in seq_len(nrow(bedpe))) {
    row <- bedpe[i, ]
    if (row$anchor_flag == "A") {
      test_g <- c(start = as.integer(row$start2) + 1, end = as.integer(row$end2))
    } else {
      test_g <- c(start = as.integer(row$start1) + 1, end = as.integer(row$end1))
    }

    test_g[1] <- test_g[1] - region_padding
    test_g[2] <- test_g[2] + region_padding

    local <- c(
      start = test_g[1] - sequence_genomic_start + 1,
      end = test_g[2] - sequence_genomic_start + 1
    )

    test_regions[[i]] <- data.frame(
      region_id = paste0("bedpe_", i),
      start = local[1],
      end = local[2],
      cut_site = local[2],
      stringsAsFactors = FALSE
    )
  }

  test_regions <- do.call(rbind, test_regions)

  if (anchor_region[1] < 1 || anchor_region[2] > sequence_n) {
    stop("BEDPE 解析出的 Anchor 区域超出 sequence 边界，请检查 sequence_genomic_start")
  }
  if (any(test_regions$start < 1 | test_regions$end > sequence_n)) {
    stop("至少一个 Test 区域超出 sequence 边界，请检查 sequence_genomic_start")
  }

  list(
    anchor_region = as.integer(anchor_region),
    anchor_cut_site = as.integer(anchor_cut_site),
    test_regions = test_regions
  )
}
