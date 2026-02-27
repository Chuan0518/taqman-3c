# taqman3c（R 包）：3C-qPCR Anchor/Test 引物-探针设计

本包按你提出的 3C 特异逻辑实现：

- **Anchor/Constant + Test 配对策略**：在恒定片段上设计一条 anchor primer + 一条 probe；在多个测试片段上批量设计 test primer。
- **单向（同链）设计**：anchor/test primer 都从同一条基因组链（正义链窗口）筛选。
- **探针反义链构型**：probe 放在 anchor primer 与 anchor 酶切位点之间，输出序列为该区段的反义链（antisense）。

## 主要规则（默认）

### 引物
- 长度：18–30 nt
- Tm：58–62°C
- GC：30%–80%
- 距离酶切位点：50–150 bp
- 同次比较中 anchor/test primer Tm 差：≤3°C
- 3'端最后 5 nt 中 G/C 数量 ≤2
- 避免同聚物（>3）与连续 `GGGG`
- 过滤明显自互补

### 探针
- 长度：20–30 nt
- Tm 比对应 anchor primer 高 8–10°C
- 5' 端禁止 G
- C 数必须 > G 数
- 与 anchor primer 距离 1–10 nt
- 不允许跨越/覆盖 anchor 酶切位点
- 过滤同聚物与明显自互补

### 扩增子（伪扩增子）
- 推荐：60–150 bp
- 上限：300 bp
- 可提供次级酶切位点 motif（如 `GAATTC`），在拼接后的嵌合模板中命中则剔除

## API

```r
devtools::load_all(".")

seq <- "ACGT..."
anchor_region <- c(100, 400)
anchor_cut_site <- 320

test_regions <- data.frame(
  region_id = c("enh1", "enh2"),
  start = c(800, 1200),
  end = c(1100, 1500),
  cut_site = c(980, 1360)
)

assays <- design_3c_taqman_assays(
  sequence = seq,
  anchor_region = anchor_region,
  anchor_cut_site = anchor_cut_site,
  test_regions = test_regions,
  secondary_enzyme_sites = c("GAATTC", "CCCGGG")
)

print_assays(assays)
```





## BEDPE 输入设计（新增）

如果你输入的是 BEDPE，可直接用 `design_3c_from_bedpe()`：

- 默认 BEDPE 列：`chr1 start1 end1 chr2 start2 end2 anchor_flag`
- 第 7 列 `anchor_flag` 规则：
  - `A`：第一个坐标对（start1/end1）是 Anchor
  - `B`：第二个坐标对（start2/end2）是 Anchor

> 例如：`chr1 5000 6000 chr1 10000 11000 A` 表示 Anchor 区间为第一个区间。

```r
assays <- design_3c_from_bedpe(
  sequence = seq,
  bedpe_path = "pairs.bedpe",
  sequence_genomic_start = 1,
  region_padding = 0,
  uniqueness_genome_fasta = "genome.fa",
  uniqueness_max_hits = 1
)
```

命令行：

```bash
Rscript inst/scripts/design_from_bedpe_cli.R   "ACGT..." pairs.bedpe 1 0 genome.fa 1
```

## 本地全基因组唯一性验证（已内置）

可以直接传入本地基因组 FASTA，程序会对引物做**精确匹配计数（正负链）**并过滤：

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

- `uniqueness_genome_fasta`: 本地 FASTA 路径（不联网）。
- `uniqueness_max_hits`: 允许命中次数上限（默认 1，表示全基因组唯一）。
- `enforce_probe_uniqueness`: 是否对探针也做同样唯一性过滤（默认仅过滤引物）。

## 说明
- 包内暂不直接集成 BLAST/BLAT 联网校验；建议将输出候选导出后做全基因组唯一性验证。
- `design_taqman_assays()` 保留为兼容别名（内部调用 `design_3c_taqman_assays()`）。



## 详细 Usage PDF

已提供详细文档：`docs/usage_guide.pdf`
文档源文件：`docs/usage_guide.md`
如需重新生成：

```bash
python3 scripts/generate_usage_pdf.py
```
