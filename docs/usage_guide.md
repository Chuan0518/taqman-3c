# taqman3c 使用说明（详细版）

版本：0.4.0  
适用场景：3C / 3C-qPCR / TaqMan Anchor-Test 引物探针设计

---

## 1. 包简介

`taqman3c` 是一个 R 包，用于 3C 场景下的 TaqMan qPCR 设计：

- 在 Anchor（恒定区域）上设计一条引物 + 一条探针；
- 在多个 Test（候选互作区域）上设计测试引物；
- 按规则过滤并评分候选；
- 可选本地 FASTA 全基因组唯一性过滤；
- 支持 BEDPE 输入（第7列 A/B 指定 Anchor 在哪一侧）。

---

## 2. 安装与环境准备

### 2.1 必需软件

- R（建议 >= 4.2）
- 推荐：`devtools`
- 测试依赖：`testthat`（在 DESCRIPTION 的 Suggests 中）

### 2.2 安装方式（开发模式）

在仓库根目录执行：

```r
install.packages("devtools")
devtools::load_all(".")
```

### 2.3 安装到本地库（可选）

```r
install.packages("devtools")
devtools::install(".")
library(taqman3c)
```

### 2.4 运行测试

```r
devtools::test()
```

或：

```bash
Rscript -e "devtools::test()"
```

---

## 3. 输入数据要求

### 3.1 序列（sequence）

- 必须是 DNA 字符串：仅允许 `A/C/G/T`。
- 所有坐标最终都映射到这条序列。

### 3.2 直接 API 模式（非 BEDPE）

需要提供：

- `anchor_region = c(start, end)`（1-based，闭区间）
- `anchor_cut_site`（位于 anchor_region 内）
- `test_regions` 数据框，包含列：
  - `region_id`
  - `start`
  - `end`
  - `cut_site`

### 3.3 BEDPE 模式

默认 BEDPE 列：

1. chr1
2. start1
3. end1
4. chr2
5. start2
6. end2
7. anchor_flag

其中 `anchor_flag` 含义：

- `A`：第一个区间是 Anchor，第2个区间是 Test；
- `B`：第二个区间是 Anchor，第1个区间是 Test。

> 坐标约定：BEDPE 输入按常见格式处理（start 0-based，end 1-based/exclusive）。

示例行：

```
chr1    5000    6000    chr1    10000   11000   A
```

---

## 4. 核心 API：`design_3c_taqman_assays()`

### 4.1 最小示例

```r
seq <- "ACGT..."

test_regions <- data.frame(
  region_id = c("enh1", "enh2"),
  start = c(800, 1200),
  end = c(1100, 1500),
  cut_site = c(980, 1360)
)

assays <- design_3c_taqman_assays(
  sequence = seq,
  anchor_region = c(100, 400),
  anchor_cut_site = 320,
  test_regions = test_regions
)

print_assays(assays)
```

### 4.2 常用参数详解

- `max_results`：返回条目上限。
- `primer_len_range`：引物长度范围，默认 `c(18, 30)`。
- `probe_len_range`：探针长度范围，默认 `c(20, 30)`。
- `primer_tm_range`：引物 Tm 范围，默认 `c(58, 62)`。
- `probe_tm_delta`：探针 Tm 比引物高的范围，默认 `c(8, 10)`。
- `primer_gc_range`：引物 GC% 范围，默认 `c(30, 80)`。
- `distance_to_cut`：引物距离酶切位点的范围，默认 `c(50, 150)`。
- `probe_gap_range`：探针与 anchor 引物距离，默认 `c(1, 10)`。
- `amplicon_size_best`：伪扩增子优选区间，默认 `c(60, 150)`。
- `amplicon_size_max`：伪扩增子最大上限，默认 `300`。
- `secondary_enzyme_sites`：次级酶切位点 motif 列表，命中即过滤。

---

## 5. BEDPE API：`design_3c_from_bedpe()`

### 5.1 快速示例

```r
assays <- design_3c_from_bedpe(
  sequence = seq,
  bedpe_path = "pairs.bedpe",
  sequence_genomic_start = 1,
  region_padding = 0
)
```

### 5.2 参数解释

- `bedpe_path`：BEDPE 文件路径。
- `sequence_genomic_start`：`sequence` 第1个碱基对应的基因组坐标（1-based）。
- `region_padding`：对 BEDPE 区间左右扩展的碱基数。
- `...`：透传给 `design_3c_taqman_assays()` 的参数（如唯一性过滤）。

### 5.3 A/B 锚点规则注意事项

- 同一个 BEDPE 文件应语义一致；
- 若第7列错误（非 A/B），函数会报错；
- 若映射后区间超出序列边界，会报错提示检查 `sequence_genomic_start`。

---

## 6. 本地全基因组唯一性过滤

### 6.1 使用示例

```r
assays <- design_3c_taqman_assays(
  sequence = seq,
  anchor_region = c(100, 400),
  anchor_cut_site = 320,
  test_regions = test_regions,
  uniqueness_genome_fasta = "genome.fa",
  uniqueness_max_hits = 1,
  enforce_probe_uniqueness = FALSE
)
```

### 6.2 规则说明

- 精确匹配计数（exact match）；
- 在正负链都计数；
- `uniqueness_max_hits=1` 表示仅允许全基因组唯一；
- 默认只过滤引物；可选过滤探针；
- 内部有缓存，避免重复计数影响性能。

---

## 7. CLI 用法

### 7.1 BEDPE CLI

```bash
Rscript inst/scripts/design_from_bedpe_cli.R \
  "ACGT..." pairs.bedpe 1 0 genome.fa 1
```

参数顺序：

1. `SEQUENCE`
2. `BEDPE_PATH`
3. `SEQUENCE_GENOMIC_START`（可选，默认1）
4. `REGION_PADDING`（可选，默认0）
5. `GENOME_FASTA`（可选）
6. `MAX_HITS`（可选，默认1）

### 7.2 区间模式 CLI

```bash
Rscript inst/scripts/design_taqman_cli.R \
  "ACGT..." 100 400 320 test_regions.tsv genome.fa 1
```

---

## 8. 输出结果字段说明

典型输出字段：

- `region_id`
- `anchor_primer_seq`, `anchor_primer_start`, `anchor_primer_end`, `anchor_primer_tm`, `anchor_primer_gc`
- `test_primer_seq`, `test_primer_start`, `test_primer_end`, `test_primer_tm`, `test_primer_gc`
- `probe_seq_antisense`, `probe_start`, `probe_end`, `probe_tm`, `probe_gc`
- `probe_gap_from_anchor_primer`
- `pseudo_amplicon_size`
- `score`（越小越优）

---

## 9. 常见报错与排查

### 9.1 `Sequence includes invalid bases`

说明序列中有非 A/C/G/T 字符。

### 9.2 `BEDPE column 7 anchor_flag must be A or B`

第七列必须严格为 A 或 B（大小写不敏感，会自动转大写）。

### 9.3 `...outside sequence bounds`

BEDPE 映射后坐标超出 `sequence` 长度。请检查：

- `sequence_genomic_start` 是否正确；
- BEDPE 坐标是否对应同一参考组装版本；
- 是否需要调整 `region_padding`。

### 9.4 结果为空（0 行）

可能原因：

- 过滤参数过严（Tm/GC/distance/uniqueness）；
- anchor/test 区间太短或离 cut-site 不满足约束；
- 探针条件（5'禁G、C>G、gap）限制导致无可行探针。

建议逐步放宽：

- `distance_to_cut`
- `primer_tm_range`
- `uniqueness_max_hits`
- `amplicon_size_max`

---

## 10. 建议工作流（生产使用）

1. 使用 BEDPE 批量生成候选；
2. 开启本地唯一性过滤作为第一层筛选；
3. 对 top 候选做离线二级结构评估；
4. 再进行外部工具/实验验证；
5. 最后确定引物探针组合并建实验板。

---

## 11. 版本记录（文档）

- v0.4.0：支持 BEDPE（A/B Anchor 语义）+ 本地唯一性过滤。

