# taqman3c（R 包）：3C-qPCR BEDPE 引物/探针设计

本包现在只保留 **BEDPE 输入模式**：`design_3c_from_bedpe()`。

## 唯一入口

```r
assays <- design_3c_from_bedpe(
  sequence = "GATC\n",   # 可直接输入字符串，末尾换行会自动清理
  bedpe_path = "pairs.bedpe",
  sequence_genomic_start = 1,
  region_padding = 0
)
print_assays(assays)
```

默认行为：
- 默认开启本地全基因组唯一性过滤；
- 默认唯一性 FASTA 路径为当前目录 `genome.fa`；
- 默认 `uniqueness_max_hits = 1`。

## BEDPE 格式

默认 7 列：

`chr1 start1 end1 chr2 start2 end2 anchor_flag`

第 7 列规则：
- `A`：第一个区间是 Anchor；
- `B`：第二个区间是 Anchor。

示例：

`chr1 5000 6000 chr1 10000 11000 A`

## CLI

```bash
Rscript inst/scripts/design_from_bedpe_cli.R   "ACGT..." pairs.bedpe 1 0 genome.fa 1
```

## 详细文档

- `docs/usage_guide.md`
- `docs/usage_guide.pdf`

可重新生成 PDF：

```bash
python3 scripts/generate_usage_pdf.py
```
