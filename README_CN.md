# Smelt

WGBS 甲基化下游数据分析 CLI 工具。读取 BISMARK 比对结果，计算位点 / 滑动窗口 / 基因组元件的甲基化水平，进行差异甲基化分析，生成汇总统计。输出 Parquet 格式供外部分析。

## 安装

```bash
git clone https://github.com/zxd06407-afk/Smelt.git && cd Smelt
uv sync
uv pip install -e .
```

**环境要求：** Python >= 3.10

## 快速入门

```bash
# 位点级甲基化分析（从 BISMARK BAM）
smelt site --input sample.bam --sample-name tumor1 -o tumor1.site.parquet

# 滑动窗口（2000bp 窗口，500bp 步长）
smelt window --input tumor1.site.parquet -w 2000 -s 500 -o tumor1.window.parquet

# 差异甲基化分析（两组比较）
smelt dmr \
  --samples tumor1=tumor1.window.parquet tumor2=tumor2.window.parquet \
           normal1=normal1.window.parquet normal2=normal2.window.parquet \
  --group1 tumor1,tumor2 --group2 normal1,normal2 \
  -o dmr_results.parquet

# 全基因组统计
smelt stats --input tumor1.site.parquet
```

## 支持的输入格式

| 格式 | 来源 | 说明 |
|------|------|------|
| SAM/BAM | BISMARK 比对 | XM 标签提供 context（Z=CpG, H=CHH, X=CHG）和甲基化状态；alignment flag 提供链信息 |
| CX_report.txt | `coverage2cytosine` | 完整信息：chr, pos, strand, meth, unmeth, context, trinucleotide |

## 子命令

| 命令 | 功能 |
|------|------|
| `smelt site` | 逐位点甲基化水平，含 CpG/CHH/CHG 分类 |
| `smelt window` | 固定窗口滑动聚合 |
| `smelt element` | GTF feature 类型甲基化 |
| `smelt metaplot` | Gene body ± 侧翼等比分箱 |
| `smelt custom` | BED 自定义区间甲基化 |
| `smelt dmr` | 窗口级 Fisher 精确检验 + BH 校正，双阈值筛选 |
| `smelt stats` | 全基因组/染色体按 context 分层统计 |

详细参数请参见 `smelt <command> --help`。

## 输出

输出为 Parquet 格式。使用 `--sample-name` 时包含 `sample` 列。
完整列说明见 `docs/superpowers/specs/2026-04-29-smelt-design.md`。

## 基准验证

Smelt 已与 WGBS 领域金标准工具 **methylKit**（R/Bioconductor）在拟南芥 met1 突变体 vs 野生型数据上进行了基准对比：

| 指标 | 结果 |
|------|------|
| DMR 重叠率（Jaccard） | CpG/CHG/CHH 均达 95% |
| DMR 甲基化差异相关性 | Pearson r = 0.998–1.000 |
| 窗口甲基化率相关性 | Pearson r = 0.982 |

完整报告：`benchmark_data/BENCHMARK_REPORT.md`

## 引用

Smelt — A WGBS Methylation Downstream Analysis Tool. (2026). https://github.com/zxd06407-afk/Smelt

## 许可证

MIT License。详见 `LICENSE` 文件。
