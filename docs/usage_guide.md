# taqman3c 使用说明（详细版）

版本：0.4.2

## 1. 目标

本包只支持一种工作模式：**BEDPE -> design_3c_from_bedpe()（唯一入口）**。

## 2. 安装

```r
install.packages("devtools")
devtools::load_all(".")
```

## 3. 输入

### 3.1 sequence

- 必须为 DNA 字符串（A/C/G/T）
- 可直接输入如 `"GATC\n"`（会自动 trim）

### 3.2 BEDPE

默认列：
`chr1 start1 end1 chr2 start2 end2 anchor_flag`

第7列：
- `A`：第一段是锚定区（Anchor）
- `B`：第二段是锚定区（Anchor）

## 4. 主函数

```r
assays <- design_3c_from_bedpe(
  sequence = "GATC\n",
  bedpe_path = "pairs.bedpe",
  sequence_genomic_start = 1,
  region_padding = 0
)
```

默认参数：
- `uniqueness_genome_fasta = "genome.fa"`
- `uniqueness_max_hits = 1`
- `enforce_probe_uniqueness = FALSE`

## 5. CLI

```bash
Rscript inst/scripts/design_from_bedpe_cli.R   "ACGT..." pairs.bedpe 1 0 genome.fa 1
```

## 6. 输出

主要字段：
- `region_id`
- `anchor_primer_seq`, `test_primer_seq`
- `probe_seq_antisense`
- `pseudo_amplicon_size`
- `score`

## 7. 常见问题

- `uniqueness_genome_fasta 文件不存在`：请提供 `genome.fa` 或显式传入路径。
- `BEDPE 第7列 anchor_flag 必须是 A 或 B`：检查第7列。
- `映射后的区间超出 sequence 边界`：检查 `sequence_genomic_start` 与 BEDPE 坐标系。
