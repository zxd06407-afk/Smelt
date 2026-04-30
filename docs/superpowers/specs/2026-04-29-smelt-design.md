# Smelt — WGBS Methylation Downstream Analysis Tool: Design Spec

## Overview

Smelt is a CLI tool for downstream analysis of WGBS (Whole Genome Bisulfite Sequencing) data. It reads BISMARK output, computes methylation levels at sites / sliding windows / genomic elements, performs differential methylation analysis, and generates genome- and chromosome-level statistics. Results are exported as TSV tables for external plotting.

- **Language:** Python 3
- **CLI Framework:** Typer
- **Dependency Management:** uv
- **Logging:** structlog
- **Versioning:** Semantic versioning, starting at 0.1.0

---

## Architecture

```
Smelt/
├── smelt/
│   ├── __init__.py
│   ├── cli.py            # Typer CLI entry point, subcommand dispatch
│   ├── io.py             # Read BISMARK cov.gz, SAM/BAM, GTF, BED, FASTA
│   ├── site.py           # Site-level methylation with sequence context
│   ├── window.py         # Fixed-window sliding-window aggregation
│   ├── element.py        # GTF feature type aggregation
│   ├── metaplot.py       # Gene body ± flanking region binning
│   ├── custom.py         # BED custom interval aggregation
│   ├── dmr.py            # Window-level Fisher exact test + BH correction
│   ├── stats.py          # Genome/chromosome-level statistics by context
│   ├── filter.py         # Coverage depth and site-count filtering
│   └── utils.py          # Coordinate conversion, sequence context helpers
├── tests/
│   ├── data/             # Test fixtures (constructed + simulated)
│   └── ...
├── scripts/
│   └── benchmark.sh      # Manual benchmark run
├── pyproject.toml
└── README.md
```

---

## Subcommands

### 1. `smelt site` — Site-Level Methylation

Read BISMARK output, classify each cytosine by sequence context (CpG / CHH / CHG), and output per-site methylation ratios.

**Supported input formats:**
- `--input`: BISMARK SAM/BAM (`*.bam`, `*.sam`) or `CX_report.txt` from `coverage2cytosine`
- SAM/BAM: context from XM tag (Z=CpG, H=CHH, X=CHG), strand from alignment flag
- CX_report.txt: context and strand read directly from file columns

**Processing:**
- Context and strand come from the input file — no FASTA-based inference
- CpG dyad merging: by default (`--merge-cpg-strands`), adjacent CpG partners at (N, +) and (N+1, -) are merged into a single dyad at position N with combined counts. CHH/CHG remain strand-separated. Use `--no-merge-cpg-strands` to disable.
- Coordinate system: BISMARK uses 1-based positions. Convert to internal 0-based representation at the I/O boundary.

**Filtering:**
- `--min-depth` (default 5): minimum read depth per cytosine

**Output:** TSV — `chr, pos, context, strand, meth, unmeth, total, ratio, sample`

---

### 2. `smelt window` — Sliding Window Methylation

Aggregate methylation signals in fixed-size sliding windows across the genome.

**Input:** `smelt site` output or directly a cov.gz (auto-runs site computation internally)

**Processing:**
- Per chromosome: slide a window of `--window` bp (default 2000) with step `--step` bp (default 500)
- Within each window, aggregate per context (CpG/CHH/CHG):
  - Total methylated / unmethylated counts
  - Number of sites (cytosines)
  - Mean methylation ratio
- Coordinate system: output windows are 0-based half-open `[start, end)`

**Filtering:**
- `--min-sites` (default 10): minimum number of cytosines per window per context

**Output:** TSV — `chr, start, end, context, n_sites, meth, unmeth, total, ratio`

---

### 3. `smelt element` — GTF Feature Aggregation

Calculate methylation levels for genomic features defined in a GTF file.

**Input:** site file + GTF

**Processing:**
- Read GTF column 3 (feature type). Filter by `--features` (comma-separated, e.g. `exon,CDS,gene`)
- Each GTF row is processed independently; no feature-level merge is performed
- For each feature row, aggregate site-level methylation per context
- GTF uses 1-based closed coordinates; convert internally to 0-based half-open before intersecting with sites

**Output:** TSV — `chr, feature_type, feature_id, context, n_sites, meth, unmeth, ratio`

---

### 4. `smelt metaplot` — Gene Body ± Flanking Region Profile

Compute methylation signal distribution across gene bodies and flanking regions, for metaplot visualization.

**Input:** site file + GTF or BED defining gene body intervals

**Processing:**
- Extract gene body intervals (start to end of each gene)
- Extend upstream by `--upstream` bp (default 2000) and downstream by `--downstream` bp (default 2000)
- Equal-ratio binning:
  - Gene body: `--body-bins` (default 50) equal bins per gene, normalized by gene length
  - Upstream: `--up-bins` (default 20) bins
  - Downstream: `--down-bins` (default 20) bins
- For each bin, aggregate methylation ratio across all genes (mean across genes)
- BED input is 0-based; GTF input is 1-based. Normalize to internal 0-based representation.

**Output:** TSV matrix — `region (upstream/body/downstream), bin, sample, ratio`

---

### 5. `smelt custom` — Custom Intervals (BED)

Calculate methylation levels for user-defined genomic intervals.

**Input:** site file + BED file

**Processing:**
- Read BED intervals (0-based half-open)
- Intersect with site-level data, aggregate by context
- Handle BED columns: standard BED3/6, also read an optional `name` column (col 4) for interval labeling

**Output:** TSV — `interval_id, chr, start, end, context, n_sites, meth, unmeth, ratio`

**Note:** BED is 0-based half-open. No coordinate conversion needed.

---

### 6. `smelt dmr` — Differential Methylation Regions

Perform window-level differential methylation analysis between two groups.

**Input:** Multiple window files (one per sample) + grouping parameters. Each sample provides its own window file; the sample name is extracted from the filename or specified explicitly.

**Input format:**
- `--samples`: space-separated list of `sample_name=file_path` pairs
  - Example: `--samples tumor1=tumor1.window.tsv tumor2=tumor2.window.tsv normal1=normal1.window.tsv`

**Processing:**
- Load all window files, tag each row with its sample name
- For each window × context combination, build a 2×2 contingency table:
  - Group 1: sum(meth) vs sum(unmeth) across all samples in the group
  - Group 2: sum(meth) vs sum(unmeth) across all samples in the group
- Fisher's exact test (two-sided) → raw p-value
- Multiple testing correction: Benjamini-Hochberg (FDR) **per sequence context** (CpG/CHH/CHG corrected separately)
  - Option `--fdr-all-contexts` to pool all contexts for global BH correction
- Calculate delta: mean_ratio_group2 − mean_ratio_group1
- DMR classification:
  - `hyper` if delta > 0 and passes thresholds (group2 hypermethylated)
  - `hypo` if delta < 0 and passes thresholds (group2 hypomethylated)

**Extreme-value handling:**
- Windows where all counts are zero (both groups): skip, output to excluded file
- Windows where total count = 0 in one group but > 0 in the other: run Fisher, but flag with `low-coverage` in excluded file
- Windows where any group has fewer than `--min-samples-per-group` samples with coverage: skip

**Significance thresholds:**
- `--q-threshold` (default 0.05): FDR-adjusted q-value cutoff
- `--delta-threshold` (default 0.2): minimum absolute methylation difference (0–1 scale)
- Both thresholds must be met for a DMR to be called significant

**Group specification:**
- `--group1 sample_a,sample_b` — sample names matching those in `--samples`
- `--group2 sample_c,sample_d`
- All sample names in groups must appear in `--samples`

**Output:**
- `sample.dmr.tsv`: `chr, start, end, context, p_value, q_value, delta, direction, ratio_group1, ratio_group2`
- `sample.dmr_excluded.tsv`: windows excluded with reason

---

### 7. `smelt stats` — Genome & Chromosome Statistics

Generate summary statistics at genome-wide and per-chromosome levels.

**Input:** site file or window file

**Processing:**
- Per chromosome × context: mean methylation ratio, total cytosine count, mean coverage depth
- Genome-wide totals and means

**Output:** TSV — `chr, context, mean_ratio, n_sites, mean_depth`

---

## Coordinate System Handling

Bioinformatics file formats disagree on coordinate conventions. Smelt enforces a strict internal rule:

### Internal convention: 0-based half-open `[start, end)`

| Format | Convention | Conversion at I/O boundary |
|--------|-----------|---------------------------|
| BISMARK cov.gz | 1-based position | `pos_zero = pos_one − 1` |
| GTF | 1-based closed `[start, end]` | `start_zero = start_one − 1; end_zero = end_one` |
| GFF | 1-based closed `[start, end]` | Same as GTF |
| BED | 0-based half-open `[start, end)` | No conversion |
| SAM/BAM | 1-based position | `pos_zero = pos_one − 1` |

### Verification rules:
- When intersecting intervals from different formats, always convert to internal representation first
- When writing output, use 0-based half-open consistently
- Add conversion helpers in `utils.py`: `to_zero_based()`, `to_one_based()`

---

## Chromosome Naming

Chromosome names must match exactly across all input files. If a chromosome name in a GTF/BED file does not appear in the coverage data (or vice versa), Smelt exits with an error listing the mismatched names. No silent data loss.

Option `--strict-names` is always on by default. No fuzzy matching (e.g. `1` ↔ `chr1`) is performed.

---

## Output File Naming

Default behavior: each subcommand generates a descriptive filename from the input basename and parameters.

- `smelt site input.cov.gz` → `input.site.tsv`
- `smelt window input.site.tsv --window 2000 --step 500` → `input.window2000s500.tsv`
- `smelt dmr input.window2000s500.tsv --group1 ... --group2 ...` → `input.dmr.tsv`

Use `--output <path>` to override.

---

## Dependency & Parallel Strategy

- Default: single-threaded, per-chromosome iteration (memory-friendly)
- Optional: `--threads N` enables multi-process parallelism, dispatching chromosomes to a process pool
- Chromosome-level processing: each chromosome's data is loaded, processed, and written independently

---

## Testing Strategy

### Unit Tests (CI)
- Hand-constructed minimal data (tens of rows per format)
- Test each function in isolation
- Cover edge cases: zero counts, single CpG, missing context, boundary positions

### Integration Tests (CI, < 2 minutes)
- Simulated dataset with ground truth:
  - 2 artificial chromosomes (chrA ~500kb, chrB ~300kb)
  - 4 samples (2 per group), mean depth 10x
  - Cytosine counts: ~2000 CpG, ~1500 CHH, ~800 CHG
  - 20 genes with varied gene body lengths (500bp–50kb)
  - 5 custom BED intervals
  - 8 pre-embedded DMRs (5 CpG, 2 CHH, 1 CHG, deltas from 0.15 to 0.8)
- Validate: DMR detection rate, methylation ratio numerical precision (< 1e-6 deviation)
- Data generation script included in repo for reproducibility

### Benchmark Tests (Manual, pre-release)
- Public benchmark dataset (from published WGBS tool-comparison papers); fallback: Arabidopsis thaliana WGBS (met1 mutant vs WT)
- Cross-tool comparison with methylKit / DSS
- Metrics: DMR overlap rate, runtime, peak memory
- Run via `scripts/benchmark.sh`

---

## Dependencies (Python)

- `typer` — CLI framework
- `pandas` — tabular data handling
- `numpy` — numerical operations
- `scipy` — Fisher exact test, BH correction
- `pysam` — SAM/BAM reading (optional, for BISMARK alignment input)
- `pyranges` or `pybedtools` — genomic interval operations
- `structlog` — structured logging
- `pytest` — testing
