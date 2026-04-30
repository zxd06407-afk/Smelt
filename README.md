# Smelt

WGBS methylation downstream analysis CLI tool. Reads BISMARK alignment output, computes methylation at sites / sliding windows / genomic elements, performs DMR calling, and generates summary statistics. Outputs TSV for external plotting.

## Installation

```bash
git clone <repo-url> && cd Smelt
uv sync
uv pip install -e .
```

## Quickstart

```bash
# Site-level methylation from BISMARK BAM
smelt site --input sample.bam --sample-name tumor1 -o tumor1.site.tsv

# Sliding windows (2000bp window, 500bp step)
smelt window --input tumor1.site.tsv -w 2000 -s 500 -o tumor1.window.tsv

# Differential methylation between groups
smelt dmr \
  --samples tumor1=tumor1.window.tsv tumor2=tumor2.window.tsv \
           normal1=normal1.window.tsv normal2=normal2.window.tsv \
  --group1 tumor1,tumor2 --group2 normal1,normal2 \
  -o dmr_results.tsv

# Genome-wide statistics
smelt stats --input tumor1.site.tsv
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

All outputs are tab-separated TSV files. Per-sample outputs include a `sample` column. See `docs/superpowers/specs/2026-04-29-smelt-design.md` for full output specifications and column descriptions.
