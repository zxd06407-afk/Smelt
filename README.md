# Smelt

WGBS methylation downstream analysis CLI tool. Reads BISMARK alignment output, computes methylation at sites / sliding windows / genomic elements, performs DMR calling, and generates summary statistics. Outputs Parquet for external analysis.

## Installation

```bash
git clone https://github.com/zxd06407-afk/Smelt.git && cd Smelt
uv sync
uv pip install -e .
```

**Requirements:** Python >= 3.10

## Quickstart

```bash
# Site-level methylation from BISMARK BAM
smelt site --input sample.bam --sample-name tumor1 -o tumor1.site.parquet

# Sliding windows (2000bp window, 500bp step)
smelt window --input tumor1.site.parquet -w 2000 -s 500 -o tumor1.window.parquet

# Differential methylation between groups
smelt dmr \
  --samples tumor1=tumor1.window.parquet tumor2=tumor2.window.parquet \
           normal1=normal1.window.parquet normal2=normal2.window.parquet \
  --group1 tumor1,tumor2 --group2 normal1,normal2 \
  -o dmr_results.parquet

# Genome-wide statistics
smelt stats --input tumor1.site.parquet
```

## Supported Input Formats

| Format | Source | Description |
|--------|--------|-------------|
| SAM/BAM | BISMARK alignment | XM tag provides context (Z=CpG, H=CHH, X=CHG) and methylation; alignment flag provides strand |
| CX_report.txt | `coverage2cytosine` | Complete: chr, pos, strand, meth, unmeth, context, trinucleotide |

## Subcommands

| Command | Purpose |
|---------|---------|
| `smelt site` | Per-cytosine methylation with context (CpG/CHH/CHG) |
| `smelt window` | Fixed-window sliding-window aggregation |
| `smelt element` | GTF feature-type methylation |
| `smelt metaplot` | Gene body +/- flanking equal-ratio binning |
| `smelt custom` | BED custom interval methylation |
| `smelt dmr` | Window-level Fisher exact + BH correction, dual threshold |
| `smelt stats` | Genome/chromosome summary by context |

See `smelt <command> --help` for detailed options.

## Output

Outputs are Parquet files. Include a `sample` column when `--sample-name` is used.
Column specifications are in `docs/superpowers/specs/2026-04-29-smelt-design.md`.

## Benchmark Validation

Smelt has been validated against **methylKit** (R/Bioconductor), the gold-standard WGBS differential methylation tool, using Arabidopsis thaliana met1 mutant vs wild-type data:

| Metric | Result |
|--------|--------|
| DMR overlap (Jaccard) | 95% across CpG/CHG/CHH |
| DMR delta correlation | Pearson r = 0.998–1.000 |
| Window methylation correlation | Pearson r = 0.982 |

Full report: `benchmark_data/BENCHMARK_REPORT.md`

## Citation

Smelt — A WGBS Methylation Downstream Analysis Tool. (2026). https://github.com/zxd06407-afk/Smelt

## License

MIT License. See `LICENSE` for details.
