# Smelt Review Fixes — Design Spec

## Overview

Follow-up modifications based on comprehensive project review. Addresses 8 issues covering documentation, design-implementation consistency, input format architecture, and biological correctness.

---

## 1. README.md

Write a complete README covering:
- Project description (WGBS methylation downstream analysis CLI)
- Installation (`uv sync`, `pip install`)
- Quickstart example
- 7 subcommands with one-line descriptions and typical usage
- Supported input formats (BAM, SAM, CX_report.txt)
- Output format descriptions
- Reference to the full design spec

---

## 2. Chromosome Name Validation

**Location:** `smelt/cli.py` — `element`, `custom`, `metaplot` commands

**Change:** After loading both site file and GTF/BED, compute the intersection of chromosome names:

```python
site_chroms = set(sites["chr"].unique())
annot_chroms = set(annot_df["chr"].unique())
common = site_chroms & annot_chroms
if not common:
    raise typer.BadParameter(
        f"No chromosomes in common between input files.\n"
        f"  Site file: {sorted(site_chroms)}\n"
        f"  Annotation: {sorted(annot_chroms)}"
    )
```

`window` and `stats` take a single site/window file — no cross-file check needed.

---

## 3. Design Doc Update: element Feature Merge

**Location:** `docs/superpowers/specs/2026-04-29-smelt-design.md`

**Change:** Update the element section to state that each GTF row is processed independently (no feature-level merge). Remove the "merge consecutive rows" requirement.

---

## 4. Add `--sample-name` to Output Modules

**Location:** `smelt/cli.py` + each compute module

**Modules affected:** `site`, `window`, `element`, `metaplot`, `custom`, `stats`

**Change:**
- Each CLI command gains a `--sample-name` option (default: inferred from input filename stem)
- The name is passed to the compute function as `sample_name`
- Output DataFrame gains a `sample` column with this value
- `dmr` is excluded — its output is already between-group, not per-sample

**Default inference:** `Path(input_file).stem.replace(".cov", "").replace(".CX_report", "")`

---

## 5. stats: Add Genome-Wide Summary Row

**Location:** `smelt/stats.py`

**Change:** After computing per-chromosome stats, append a genome-wide row:

```python
genome_row = {
    "chr": "genome",
    "context": ctx,
    "mean_ratio": weighted_mean,  # weighted by n_sites
    "n_sites": total_sites,
    "mean_depth": overall_mean_depth,
}
```

Output: same columns, one extra row per context with `chr = "genome"`.

---

## 6. Input Format Architecture

### Removed Support
- `cov.gz` as standalone input (no strand, no context)

### New Supported Inputs

| Format | Source | Strand | Context |
|--------|--------|--------|---------|
| SAM/BAM | BISMARK alignment | alignment flag (0x10) | XM tag (Z=CG, H=CHH, X=CHG) |
| CX_report.txt | `coverage2cytosine` | column 3 | column 6 |

### Changes Required

**`smelt/io.py`:**
- Remove `read_cov()` or deprecate it
- Update `read_bismark_sam()` to extract strand from alignment flag and context from XM tag directly
- Add `read_cx_report()` function

**`smelt/site.py`:**
- `compute_site_methylation()` no longer determines context from FASTA
- Context and strand come from the input DataFrame directly (via `io.py` readers)
- `fasta` parameter removed; `context_df` parameter removed
- CpG strand merging uses correct strand info

**`smelt/cli.py`:**
- `smelt site --input` accepts `.bam`, `.sam`, `.CX_report.txt`
- `--fasta` option removed from `site` command

### CX_report.txt Format

BISMARK `coverage2cytosine` output:
```
<chromosome> <position> <strand> <count methylated> <count unmethylated> <context> <trinucleotide context>
```

Columns: chr, pos (1-based), strand (+/-), meth, unmeth, context (CG/CHH/CHG), trinucleotide

Reader converts to internal 0-based format.

---

## 7. DMR Exclusion Category Column

**Location:** `smelt/dmr.py`

**Change:** Add `category` column to excluded DataFrame:

| category | Criteria |
|----------|----------|
| `low_coverage` | Fewer than `min_samples_per_group` samples with coverage in either group |
| `all_zero` | Both groups have zero total coverage |
| `no_signal` | Tested but failed q-value or delta thresholds |
| `fisher_failed` | Fisher exact test returned NaN |

Existing `reason` column preserved with human-readable detail.

---

## 8. CpG Dyad Merging Fix

**Location:** `smelt/utils.py` — `merge_cpg_strands()`

**Problem:** Current merge groups by `(chr, pos)`, but CpG dyad partners are at adjacent positions (N on + strand, N+1 on - strand).

**Fix:** New algorithm:
1. Separate CpG records by strand
2. For each + strand CpG at position N, look for - strand CpG at position N+1
3. Merge found pairs: `meth = meth_N + meth_{N+1}`, `unmeth = unmeth_N + unmeth_{N+1}`, position = N
4. Unmatched CpG sites (no partner) pass through as-is
5. CHH/CHG sites pass through unchanged

Strand information is now reliable (from BAM/CX_report), so this merging is correct.

---

## Implementation Order

1. Input format refactor (#6) — foundational, everything depends on it
2. CpG dyad fix (#8) — depends on correct strand from #6
3. Sample name column (#4)
4. Stats genome row (#5)
5. Chromosome validation (#2)
6. DMR exclusion category (#7)
7. Design doc update (#3)
8. README (#1)
