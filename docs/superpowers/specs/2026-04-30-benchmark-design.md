# Smelt Benchmark — Design Spec

## Overview

Benchmark Smelt against methylKit (R/Bioconductor) — the most widely used WGBS differential methylation tool — using Arabidopsis thaliana met1 mutant vs WT data. Compare methylation levels, window aggregation, DMR calling, and GTF element annotation.

## Data

- **Source:** PRJEB9919 (Stroud et al. 2013, Cell)
- **Genome:** TAIR10 (135 Mb)
- **Samples:** Col-0 rep1 (ERR965674), Col-0 rep2 (ERR965675), met1-3 rep1 (ERR965676), met1-3 rep2 (ERR965677)
- **Input for both tools:** BISMARK `*bismark.cov.gz` files (from `bismark_methylation_extractor`)

## Comparison Tool

- **methylKit** (R/Bioconductor, v1.28+)
- Key functions: `read.bismark()`, `filterByCoverage()`, `tileMethylCounts()`, `calculateDiffMeth()`, `annotate.WithGenicParts()`
- Same statistical approach: Fisher exact test + BH correction per context

## Parameter Alignment

| Semantic | Smelt | methylKit | Aligned value |
|----------|-------|-----------|---------------|
| Min read depth | `--min-depth 5` | `filterByCoverage(lo.count=5)` | 5 |
| Window size | `--window 2000` | `tileMethylCounts(win.size=2000)` | 2000 bp |
| Window step | `--step 500` | `tileMethylCounts(step.size=500)` | 500 bp |
| DMR statistical test | Fisher exact | `calculateDiffMeth(method="Fisher")` | Fisher |
| DMR FDR threshold | `--q-threshold 0.05` | `getMethylDiff(qvalue=0.05)` | 0.05 |
| DMR delta threshold | `--delta-threshold 0.2` | `getMethylDiff(difference=20)` | 20% (=0.2) |
| CpG merge | dyad merge (N+/N+1-) | `read.bismark()` handles per-strand | per-tool default |

---

## Benchmark Modules

### 1. Site-Level Methylation (`smelt site` vs `methylKit::read.bismark()`)

**Input:** Same BISMARK `*bismark.cov.gz` files.

**Metrics:**
- Number of covered cytosines per context (CpG/CHH/CHG)
- Per-site methylation ratio correlation (Pearson r, Spearman ρ)
- Per-context correlation

**Output:** `benchmark/01_site.json`

### 2. Window-Level Methylation (`smelt window` vs `methylKit::tileMethylCounts()`)

**Metrics:**
- Window count per sample per context
- Window methylation ratio correlation (Pearson r, Spearman ρ)
- Per-context correlation

**Output:** `benchmark/02_window.json`

### 3. DMR Calling (`smelt dmr` vs `methylKit::calculateDiffMeth()`)

**Method:** Both tools run Fisher exact test on window-level meth/unmeth counts + BH correction.

**Metrics:**
- Number of significant DMRs per context (q<0.05, |delta|>0.2)
- Overlap: Jaccard index of DMR regions
- Direction agreement: hyper/hypo concordance in overlapping DMRs
- Delta correlation in overlapping windows

**Output:** `benchmark/03_dmr.json`

### 4. Element Annotation (`smelt element` vs `methylKit::annotate.WithGenicParts()`)

**Input:** Same GTF annotation file.

**Metrics:**
- Feature-type methylation levels correlation
- Per-feature methylation ratio agreement

**Output:** `benchmark/04_element.json`

### 5. Genome/Chromosome Statistics (`smelt stats`)

**Metrics:**
- Per-chromosome mean methylation ratios
- Genome-wide weighted means per context
- Comparison with methylKit's `getData()` aggregate statistics

**Output:** `benchmark/05_stats.json`

---

## Benchmark Script

R script at `scripts/benchmark_methylkit.R` runs the methylKit pipeline.
Python script at `scripts/benchmark.py` extended to compare Smelt vs methylKit outputs.
Results saved to `benchmark/benchmark_report.json`.

---

## Implementation Order

1. Install methylKit R package
2. Write R benchmark script (site → window → DMR → element)
3. Run methylKit pipeline
4. Write Python comparison script
5. Run comparison and generate report
