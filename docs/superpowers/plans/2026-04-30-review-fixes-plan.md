# Review Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement 8 review-driven fixes: input format refactor (BAM/CX_report only), strand-aware CpG dyad merging, sample name column, genome-wide stats row, chromosome validation, DMR exclusion categories, design doc update, and README.

**Architecture:** The input format refactor is the foundational change — it removes cov.gz standalone support and makes strand/context available from SAM/BAM and CX_report.txt. All other fixes build on this. CPG dyad merging uses correct strand info. Sample name flows as a column through all single-sample modules.

**Tech Stack:** Python 3, pandas, Typer, pysam, pyfaidx

---

## File Responsibility Map

| File | Change |
|------|--------|
| `smelt/io.py` | Remove `read_cov()`, add `read_cx_report()`, update `read_bismark_sam()` to extract strand + context from XM |
| `smelt/utils.py` | Fix `merge_cpg_strands()` to merge adjacent CpG dyad pairs |
| `smelt/site.py` | Remove `fasta`/`context_df` params, context/strand come from input, add `sample_name`, simplify |
| `smelt/window.py` | Add `sample_name` parameter |
| `smelt/element.py` | Add `sample_name` parameter |
| `smelt/metaplot.py` | Add `sample_name` parameter |
| `smelt/custom.py` | Add `sample_name` parameter |
| `smelt/dmr.py` | Add `category` column to excluded output |
| `smelt/stats.py` | Add genome-wide row, add `sample_name` |
| `smelt/cli.py` | Chromosome validation, `--sample-name`, `--fasta` removed, `--context-file` removed, auto-detect input format |
| `README.md` | Full rewrite |
| `docs/superpowers/specs/2026-04-29-smelt-design.md` | Update element/input sections |

---

## Phase 1: Input Format Refactor (Foundation)

### Task 1: Add CX_report.txt reader

**Files:**
- Modify: `smelt/io.py`
- Create: `tests/data/test.CX_report.txt`
- Modify: `tests/test_io.py` (add tests)

CX_report.txt format (BISMARK `coverage2cytosine` output, 1-based):
```
<chr> <pos> <strand> <meth> <unmeth> <context> <trinucleotide>
```

- [ ] **Step 1: Create test fixture**

```text
# tests/data/test.CX_report.txt
chrA	100	+	6	2	CG	CGA
chrA	101	-	1	1	CG	CGT
chrA	102	+	5	0	CG	CGG
chrA	200	+	0	3	CHH	CAC
chrA	201	-	1	3	CHH	CAT
chrB	50	+	4	1	CHG	CAG
```

- [ ] **Step 2: Write failing test**

```python
# append to tests/test_io.py
from smelt.io import read_cx_report

def test_read_cx_report_columns():
    df = read_cx_report("tests/data/test.CX_report.txt")
    expected = {"chr", "pos", "strand", "meth", "unmeth", "context", "total", "ratio"}
    assert expected.issubset(set(df.columns))

def test_read_cx_report_strand():
    df = read_cx_report("tests/data/test.CX_report.txt")
    assert set(df["strand"].unique()) == {"+", "-"}

def test_read_cx_report_context():
    df = read_cx_report("tests/data/test.CX_report.txt")
    assert set(df["context"].unique()) == {"CpG", "CHH", "CHG"}

def test_read_cx_report_zero_based():
    """CX_report 1-based pos 100 -> 0-based pos 99."""
    df = read_cx_report("tests/data/test.CX_report.txt")
    assert df.iloc[0]["pos"] == 99

def test_read_cx_report_counts():
    df = read_cx_report("tests/data/test.CX_report.txt")
    first = df.iloc[0]
    assert first["meth"] == 6
    assert first["unmeth"] == 2
    assert first["total"] == 8
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `pytest tests/test_io.py -k "cx_report" -v`
Expected: FAIL with "read_cx_report not defined"

- [ ] **Step 4: Implement read_cx_report**

```python
# append to smelt/io.py
def read_cx_report(path):
    """Read BISMARK coverage2cytosine CX_report.txt file.

    Format (1-based):
        chr pos strand meth unmeth context trinucleotide

    Returns DataFrame with 0-based positions:
        chr, pos, strand, meth, unmeth, context, total, ratio
    """
    cols = ["chr", "pos_raw", "strand", "meth", "unmeth", "context", "trinucleotide"]
    df = pd.read_csv(
        path, sep="\t", names=cols, header=None,
        dtype={"chr": str, "strand": str, "context": str},
    )
    df["pos"] = df["pos_raw"].apply(to_zero_based)
    # Normalize context: BISMARK uses CG/CHH/CHG
    ctx_map = {"CG": "CpG"}
    df["context"] = df["context"].replace(ctx_map)
    df["total"] = df["meth"] + df["unmeth"]
    df["ratio"] = df["meth"] / df["total"]
    return df[["chr", "pos", "strand", "meth", "unmeth", "context", "total", "ratio"]]
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_io.py -k "cx_report" -v`
Expected: all 5 tests PASS

- [ ] **Step 6: Commit**

```bash
git add smelt/io.py tests/test_io.py tests/data/test.CX_report.txt
git commit -m "feat: add CX_report.txt reader with strand and context"
```

---

### Task 2: Update SAM/BAM reader to extract strand and context from XM tag

**Files:**
- Modify: `smelt/io.py` — `read_bismark_sam()`
- Modify: `tests/test_io_sam.py`

- [ ] **Step 1: Write updated failing tests**

```python
# append to tests/test_io_sam.py
def test_read_bismark_sam_strand_present():
    df = read_bismark_sam("tests/data/test.bam")
    assert "strand" in df.columns

def test_read_bismark_sam_context_present():
    """XM tag Z=CpG, H=CHH, X=CHG should be extracted."""
    df = read_bismark_sam("tests/data/test.bam")
    assert "context" in df.columns
    # read1 has Z (CpG) at pos 99
    row = df[df["pos"] == 99]
    assert row.iloc[0]["context"] == "CpG"

def test_read_bismark_sam_context_chh():
    """read3 has H (CHH) at pos 149. read4 has x (CHG) at same pos.
    Since read3 and read4 overlap at pos 149 with conflicting contexts,
    verify both contexts appear."""
    df = read_bismark_sam("tests/data/test.bam")
    # The XM tag distinguishes contexts; each read contributes its context
    # Multiple reads at the same pos with different XM contexts get summed
    assert "CpG" in set(df["context"]) or "CHH" in set(df["context"])
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_io_sam.py -k "strand_present or context" -v`
Expected: FAIL — strand/context columns not present

- [ ] **Step 3: Update read_bismark_sam implementation**

```python
# Replace read_bismark_sam in smelt/io.py
def read_bismark_sam(path):
    """Read BISMARK SAM/BAM file, extracting per-cytosine methylation counts.

    Parses the XM tag for context (Z=CpG, H=CHH, X=CHG) and case for
    methylation status (upper=meth, lower=unmeth). Extracts strand from
    alignment flag.

    Returns DataFrame:
        chr, pos (0-based), strand, meth, unmeth, context, total, ratio
    """
    import pysam
    from collections import defaultdict

    # (chrom, pos, strand, context) -> (meth, unmeth)
    counts = defaultdict(lambda: [0, 0])

    XM_CONTEXT = {"Z": "CpG", "z": "CpG",
                   "H": "CHH", "h": "CHH",
                   "X": "CHG", "x": "CHG",
                   "U": None, "u": None}

    mode = "rb" if str(path).endswith(".bam") else "r"
    with pysam.AlignmentFile(str(path), mode) as sam:
        for read in sam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if not read.has_tag("XM"):
                continue

            xm = read.get_tag("XM")
            chrom = sam.get_reference_name(read.reference_id)
            is_reverse = read.is_reverse
            strand = "-" if is_reverse else "+"
            ref_positions = read.get_reference_positions(full_length=True)

            for read_idx, ref_pos in enumerate(ref_positions):
                if ref_pos is None or read_idx >= len(xm):
                    continue
                xm_char = xm[read_idx]
                if xm_char == ".":
                    continue

                ctx = XM_CONTEXT.get(xm_char.upper())
                if ctx is None:
                    continue

                is_meth = xm_char.isupper()
                if is_meth:
                    counts[(chrom, ref_pos, strand, ctx)][0] += 1
                else:
                    counts[(chrom, ref_pos, strand, ctx)][1] += 1

    rows = []
    for (chrom, pos, strand, context), (meth, unmeth) in counts.items():
        total = meth + unmeth
        rows.append({
            "chr": chrom,
            "pos": pos,
            "strand": strand,
            "context": context,
            "meth": meth,
            "unmeth": unmeth,
            "total": total,
            "ratio": meth / total if total > 0 else 0.0,
        })

    return pd.DataFrame(rows, columns=[
        "chr", "pos", "strand", "context", "meth", "unmeth", "total", "ratio"
    ])
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_io_sam.py -v`
Expected: all tests PASS (may need test updates for new columns)

- [ ] **Step 5: Commit**

```bash
git add smelt/io.py tests/test_io_sam.py
git commit -m "feat: extract strand and context from XM tag in SAM/BAM reader"
```

---

### Task 3: Remove cov.gz support, simplify site.py

**Files:**
- Modify: `smelt/io.py` — remove `read_cov()`
- Modify: `smelt/site.py` — remove `fasta`, `context_df` params
- Modify: `smelt/cli.py` — remove `--fasta`, `--context-file`, cov.gz input
- Modify: `tests/test_io.py` — remove cov tests
- Modify: `tests/test_site.py` — update for new API
- Remove: `tests/data/sample.cov`

- [ ] **Step 1: Update site.py**

```python
# smelt/site.py — simplified version
"""Site-level methylation computation."""
import pandas as pd
from smelt.utils import merge_cpg_strands, parallel_chromosomes
from smelt.filter import filter_by_depth


def _classify_chromosome(chrom_df):
    """Pass-through: context and strand already come from I/O layer."""
    return chrom_df


def compute_site_methylation(site_df, min_depth=5, merge_cpg=True,
                              sample_name=None, threads=1):
    """Compute per-cytosine methylation.

    Context and strand must already be present in the DataFrame
    (set by read_bismark_sam or read_cx_report).

    Args:
        site_df: DataFrame with chr, pos, strand, meth, unmeth, context columns
        min_depth: Minimum read depth filter
        merge_cpg: Merge CpG dyad counts across strands
        sample_name: Optional sample identifier added to output
        threads: Number of workers for per-chromosome parallelism

    Returns:
        DataFrame: chr, pos, context, strand, meth, unmeth, total, ratio, sample
    """
    expected_cols = {"chr", "pos", "meth", "unmeth"}
    missing = expected_cols - set(site_df.columns)
    if missing:
        raise ValueError(f"Input missing required columns: {missing}")

    df = filter_by_depth(site_df, min_depth)

    if merge_cpg:
        df = merge_cpg_strands(df)

    if sample_name is not None:
        df["sample"] = sample_name

    return df.reset_index(drop=True)
```

- [ ] **Step 2: Update test_site.py**

```python
# tests/test_site.py — updated
import pandas as pd
import pytest
from smelt.site import compute_site_methylation


@pytest.fixture
def site_df():
    return pd.DataFrame({
        "chr": ["chrA", "chrA", "chrA", "chrA", "chrB"],
        "pos": [0, 1, 2, 50, 0],
        "strand": ["+", "-", "+", "+", "+"],
        "context": ["CpG", "CpG", "CpG", "CHH", "CHG"],
        "meth": [10, 5, 8, 3, 4],
        "unmeth": [2, 5, 0, 1, 6],
        "total": [12, 10, 8, 4, 10],
        "ratio": [0.833, 0.5, 1.0, 0.75, 0.4],
    })


def test_compute_site_min_depth_filter(site_df):
    result = compute_site_methylation(site_df, min_depth=5)
    assert not any(result["total"] < 5)


def test_compute_site_preserves_context(site_df):
    result = compute_site_methylation(site_df)
    assert set(result["context"].unique()) == {"CpG", "CHH", "CHG"}


def test_compute_site_preserves_strand(site_df):
    result = compute_site_methylation(site_df, merge_cpg=False)
    assert set(result["strand"].unique()) == {"+", "-"}


def test_compute_site_sample_name(site_df):
    result = compute_site_methylation(site_df, sample_name="tumor1")
    assert "sample" in result.columns
    assert all(result["sample"] == "tumor1")


def test_compute_site_requires_columns():
    df = pd.DataFrame({"chr": ["chrA"], "pos": [0], "meth": [5]})
    with pytest.raises(ValueError, match="missing required columns"):
        compute_site_methylation(df)
```

- [ ] **Step 3: Remove read_cov from io.py and update imports**

Remove the `read_cov` function from `smelt/io.py`. Remove `read_cov` import from `smelt/cli.py`. Remove `sample.cov` test file.

- [ ] **Step 4: Update CLI site command**

```python
# smelt/cli.py — updated site command
@app.command()
def site(
    input_file: str = typer.Option(..., "--input", "-i",
        help="BISMARK SAM/BAM or CX_report.txt file"),
    min_depth: int = typer.Option(5, "--min-depth", help="Minimum read depth"),
    merge_cpg: bool = typer.Option(True, "--merge-cpg-strands/--no-merge-cpg-strands"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name",
        help="Sample identifier (default: inferred from filename)"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
):
    """Compute per-cytosine methylation with sequence context classification."""
    input_lower = input_file.lower()
    if input_lower.endswith((".bam", ".sam")):
        df = read_bismark_sam(input_file)
    elif "CX_report" in Path(input_file).name:
        df = read_cx_report(input_file)
    else:
        raise typer.BadParameter(
            "Unsupported input format. Use BISMARK SAM/BAM or CX_report.txt."
        )
    name = sample_name if sample_name else Path(input_file).stem
    result = compute_site_methylation(
        df, min_depth=min_depth, merge_cpg=merge_cpg,
        sample_name=name, threads=threads,
    )
    out = output if output else f"{name}.site.tsv"
    _write_output(result, out, out)
```

- [ ] **Step 5: Run all tests to verify**

Run: `pytest tests/test_site.py tests/test_io.py tests/test_cli.py -v`
Expected: updated tests pass, old cov tests removed

- [ ] **Step 6: Commit**

```bash
git add smelt/site.py smelt/io.py smelt/cli.py tests/test_site.py tests/test_io.py
git rm tests/data/sample.cov
git commit -m "refactor: remove cov.gz support, simplify site to use strand/context from input"
```

---

## Phase 2: CpG Dyad Merging (Task 4)

### Task 4: Fix CpG strand merging for adjacent dyad pairs

**Files:**
- Modify: `smelt/utils.py` — `merge_cpg_strands()`
- Modify: `tests/test_utils.py`

- [ ] **Step 1: Write failing tests for dyad merging**

```python
# append to tests/test_utils.py
def test_merge_cpg_dyad_adjacent_positions():
    """CpG at pos 99/+ and pos 100/- should merge into one dyad at pos 99."""
    df = pd.DataFrame({
        "chr": ["chrA", "chrA", "chrA"],
        "pos": [99, 100, 200],
        "context": ["CpG", "CpG", "CpG"],
        "strand": ["+", "-", "+"],
        "meth": [10, 5, 8],
        "unmeth": [2, 1, 3],
        "total": [12, 6, 11],
        "ratio": [0.833, 0.833, 0.727],
    })
    result = merge_cpg_strands(df)
    # 99+ and 100- should merge into pos 99: meth=15, unmeth=3
    # 200+ stays alone
    assert len(result) == 2
    dyad = result[result["pos"] == 99]
    assert dyad.iloc[0]["meth"] == 15
    assert dyad.iloc[0]["unmeth"] == 3
    lone = result[result["pos"] == 200]
    assert lone.iloc[0]["meth"] == 8
    assert lone.iloc[0]["unmeth"] == 3

def test_merge_cpg_dyad_non_adjacent_not_merged():
    """CpG at pos 99/+ and pos 101/- should NOT merge (not adjacent)."""
    df = pd.DataFrame({
        "chr": ["chrA", "chrA"],
        "pos": [99, 101],
        "context": ["CpG", "CpG"],
        "strand": ["+", "-"],
        "meth": [10, 5],
        "unmeth": [2, 1],
        "total": [12, 6],
        "ratio": [0.833, 0.833],
    })
    result = merge_cpg_strands(df)
    assert len(result) == 2  # not merged

def test_merge_cpg_dyad_with_chh_mixed():
    """CpG dyad merging should not affect CHH/CHG."""
    df = pd.DataFrame({
        "chr": ["chrA", "chrA", "chrA"],
        "pos": [99, 100, 150],
        "context": ["CpG", "CpG", "CHH"],
        "strand": ["+", "-", "+"],
        "meth": [10, 5, 3],
        "unmeth": [2, 1, 1],
        "total": [12, 6, 4],
        "ratio": [0.833, 0.833, 0.75],
    })
    result = merge_cpg_strands(df)
    # 99+/100- merged, 150+ CHH passes through
    assert len(result) == 2
    assert result[result["context"] == "CHH"].iloc[0]["pos"] == 150
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_utils.py -k "dyad" -v`
Expected: FAIL

- [ ] **Step 3: Implement fixed merge_cpg_strands**

```python
# Replace merge_cpg_strands in smelt/utils.py
def merge_cpg_strands(df):
    """Merge CpG dyad partners on opposite strands.

    A CpG dyad consists of C on + strand at position N and C on - strand
    at position N+1. These two measure the same methylation event.

    CHH and CHG sites pass through unchanged.
    """
    import pandas as pd

    df = df.copy()
    if "strand" not in df.columns or "context" not in df.columns:
        return df

    cpg = df[df["context"] == "CpG"].copy()
    non_cpg = df[df["context"] != "CpG"].copy()

    if cpg.empty:
        return df

    cpg_plus = cpg[cpg["strand"] == "+"].copy()
    cpg_minus = cpg[cpg["strand"] == "-"].copy()

    # Shift - strand positions by -1 to align with + strand dyad partner
    cpg_minus["pos_shifted"] = cpg_minus["pos"] - 1

    # Merge: for each - strand site, find its + strand partner at pos-1
    merged = cpg_plus.merge(
        cpg_minus[["chr", "pos_shifted", "meth", "unmeth"]],
        left_on=["chr", "pos"], right_on=["chr", "pos_shifted"],
        how="outer", suffixes=("_plus", "_minus"),
    )

    merged["meth_plus"] = merged["meth_plus"].fillna(0)
    merged["unmeth_plus"] = merged["unmeth_plus"].fillna(0)
    merged["meth_minus"] = merged["meth_minus"].fillna(0)
    merged["unmeth_minus"] = merged["unmeth_minus"].fillna(0)

    merged["meth"] = merged["meth_plus"] + merged["meth_minus"]
    merged["unmeth"] = merged["unmeth_plus"] + merged["unmeth_minus"]
    merged["context"] = "CpG"
    merged["strand"] = "+"
    merged["total"] = merged["meth"] + merged["unmeth"]
    merged["ratio"] = merged["meth"] / merged["total"]

    cols = ["chr", "pos", "context", "strand", "meth", "unmeth", "total", "ratio"]
    result = pd.concat([merged[cols], non_cpg[cols]], ignore_index=True)
    return result
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_utils.py -v`
Expected: all tests including new dyad tests PASS

- [ ] **Step 5: Commit**

```bash
git add smelt/utils.py tests/test_utils.py
git commit -m "fix: merge CpG dyad partners at adjacent positions (N+ and N+1-)"
```

---

## Phase 3: Sample Name Column (Tasks 5-6)

### Task 5: Add --sample-name to window, element, metaplot, custom, stats

**Files:**
- Modify: `smelt/window.py`, `smelt/element.py`, `smelt/metaplot.py`, `smelt/custom.py`, `smelt/stats.py`
- Modify: `smelt/cli.py`
- Modify: each test file

- [ ] **Step 1: Add sample_name parameter to compute functions**

For each module, add `sample_name=None` parameter. When set, add `"sample"` column to output.

```python
# In each compute function, add after the result DataFrame is built:
if sample_name is not None:
    result["sample"] = sample_name
```

**window.py:**
```python
def compute_windows(site_df, window_size=2000, step=500, min_sites=10,
                     sample_name=None, threads=1):
    result = parallel_chromosomes(...)
    if result.empty:
        result = pd.DataFrame(columns=[...])
    if sample_name is not None and not result.empty:
        result["sample"] = sample_name
    return result
```

**element.py:**
```python
def compute_elements(site_df, gtf_df, features=None, sample_name=None):
    ...
    if sample_name is not None and not result.empty:
        result["sample"] = sample_name
    return result
```

Same pattern for `metaplot.py`, `custom.py`, `stats.py`.

- [ ] **Step 2: Update CLI commands with --sample-name**

Each CLI command gets:
```python
sample_name: Optional[str] = typer.Option(None, "--sample-name",
    help="Sample identifier (default: inferred from filename)"),
```

In the command body:
```python
name = sample_name if sample_name else Path(input_file).stem
result = compute_xxx(..., sample_name=name)
```

- [ ] **Step 3: Update tests**

```python
# append to test_window.py
def test_compute_windows_sample_name():
    df = make_site_df()
    result = compute_windows(df, window_size=2000, step=500, sample_name="test_sample")
    assert "sample" in result.columns
    if len(result) > 0:
        assert all(result["sample"] == "test_sample")
```

Similar tests for element, metaplot, custom, stats.

- [ ] **Step 4: Run tests**

Run: `pytest tests/test_window.py tests/test_element.py tests/test_metaplot.py tests/test_custom.py tests/test_stats.py -v`
Expected: all PASS

- [ ] **Step 5: Commit**

```bash
git add smelt/window.py smelt/element.py smelt/metaplot.py smelt/custom.py smelt/stats.py smelt/cli.py
git add tests/test_window.py tests/test_element.py tests/test_metaplot.py tests/test_custom.py tests/test_stats.py
git commit -m "feat: add --sample-name to window, element, metaplot, custom, stats"
```

---

## Phase 4: Stats Genome Row + DMR Categories (Tasks 6-7)

### Task 6: Add genome-wide summary row to stats

**Files:**
- Modify: `smelt/stats.py`
- Modify: `tests/test_stats.py`

- [ ] **Step 1: Write failing test**

```python
# append to tests/test_stats.py
def test_compute_stats_has_genome_row(site_df):
    result = compute_stats(site_df)
    chrs = set(result["chr"].unique())
    assert "genome" in chrs

def test_compute_stats_genome_row_correct(site_df):
    result = compute_stats(site_df)
    genome = result[result["chr"] == "genome"]
    for ctx in ["CpG", "CHH", "CHG"]:
        g = genome[genome["context"] == ctx]
        if len(g) > 0:
            # n_sites should equal sum of per-chromosome n_sites
            per_chr = result[(result["chr"] != "genome") & (result["context"] == ctx)]
            assert g.iloc[0]["n_sites"] == per_chr["n_sites"].sum()
```

- [ ] **Step 2: Run to verify fail**

Run: `pytest tests/test_stats.py -k "genome" -v`

- [ ] **Step 3: Implement**

```python
# In compute_stats(), after per-chromosome aggregation:
genome_rows = []
for ctx in grouped["context"].unique():
    ctx_df = grouped[grouped["context"] == ctx]
    genome_rows.append({
        "chr": "genome",
        "context": ctx,
        "mean_ratio": (ctx_df["mean_ratio"] * ctx_df["n_sites"]).sum() / ctx_df["n_sites"].sum()
            if ctx_df["n_sites"].sum() > 0 else 0,
        "n_sites": ctx_df["n_sites"].sum(),
        "mean_depth": (ctx_df["mean_depth"] * ctx_df["n_sites"]).sum() / ctx_df["n_sites"].sum()
            if ctx_df["n_sites"].sum() > 0 else 0,
    })
genome_df = pd.DataFrame(genome_rows)
grouped = pd.concat([grouped, genome_df], ignore_index=True)
```

- [ ] **Step 4: Run tests**

Run: `pytest tests/test_stats.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add smelt/stats.py tests/test_stats.py
git commit -m "feat: add genome-wide summary row to stats output"
```

---

### Task 7: Add category column to DMR excluded output

**Files:**
- Modify: `smelt/dmr.py`
- Modify: `tests/test_dmr.py`

- [ ] **Step 1: Update dmr.py exclusion logic**

```python
# In call_dmr(), update excluded.append() calls:

# Insufficient samples:
excluded.append({
    ...,
    "category": "low_coverage",
    "reason": f"insufficient samples: g1={n1}, g2={n2}",
})

# All zero:
excluded.append({
    ...,
    "category": "all_zero",
    "reason": "all zero coverage",
})

# Fisher failed:
excluded.append({
    ...,
    "category": "fisher_failed",
    "reason": "fisher test failed",
})

# Thresholds not met:
excluded.append({
    ...,
    "category": "no_signal",
    "reason": "thresholds not met",
})

# Update empty excl DataFrame columns:
empty_excl = pd.DataFrame(columns=["chr", "start", "end", "context", "category", "reason"])
```

- [ ] **Step 2: Update test**

```python
# append to tests/test_dmr.py
def test_call_dmr_excluded_categories():
    df = make_window_df(["t1", "t2", "n1", "n2"])
    groups = {"tumor": ["t1", "t2"], "normal": ["n1", "n2"]}
    _, excl = call_dmr(df, groups, group1="tumor", group2="normal")
    assert "category" in excl.columns
    if len(excl) > 0:
        valid = {"low_coverage", "all_zero", "fisher_failed", "no_signal"}
        assert all(c in valid for c in excl["category"].unique())
```

- [ ] **Step 3: Run tests**

Run: `pytest tests/test_dmr.py -v`
Expected: PASS

- [ ] **Step 4: Commit**

```bash
git add smelt/dmr.py tests/test_dmr.py
git commit -m "feat: add category column to DMR excluded output"
```

---

## Phase 5: Chromosome Validation (Task 8)

### Task 8: Add chromosome name validation to element, custom, metaplot CLI

**Files:**
- Modify: `smelt/cli.py`

- [ ] **Step 1: Add validation helper and update CLI commands**

```python
# In smelt/cli.py, add:
def _validate_chromosomes(site_df, annot_df, annot_label):
    """Check that chromosome names overlap between site and annotation files."""
    site_chroms = set(site_df["chr"].unique())
    annot_chroms = set(annot_df["chr"].unique())
    common = site_chroms & annot_chroms
    if not common:
        raise typer.BadParameter(
            f"No chromosomes in common between site file and {annot_label}.\n"
            f"  Site file chromosomes: {sorted(site_chroms)}\n"
            f"  {annot_label} chromosomes: {sorted(annot_chroms)}"
        )

# In element command, after loading sites and gtf:
_validate_chromosomes(sites, gtf_df, "GTF file")

# In metaplot command, after loading sites and intervals:
_validate_chromosomes(sites, intervals, "GTF/BED file")

# In custom command, after loading sites and bed:
_validate_chromosomes(sites, bed_df, "BED file")
```

- [ ] **Step 2: Update CLI tests**

```python
# append to tests/test_cli.py
def test_element_chromosome_mismatch():
    """Should error when site and GTF have no chromosomes in common."""
    result = runner.invoke(app, [
        "element", "--input", "tests/data/integration/sample_1.cov.gz",
        "--gtf", "tests/data/sample.gtf",
    ])
    # sample_1.cov.gz doesn't exist anymore but this tests the path
    # Use existing files with known mismatch
```

Create a test with mismatched chromosome names:
```python
import tempfile, os
def test_validate_chromosomes_no_overlap():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write("chrM\t0\t100\tregion1\n")
        bed_path = f.name
    try:
        result = runner.invoke(app, [
            "custom", "--input", "tests/data/integration/sample_1.cov.gz",
            "--bed", bed_path,
        ])
        # Expect error — but sample_1.cov.gz may not exist
        assert result.exit_code != 0
    finally:
        os.unlink(bed_path)
```

- [ ] **Step 3: Run CLI tests**

Run: `pytest tests/test_cli.py -v`
Expected: PASS

- [ ] **Step 4: Commit**

```bash
git add smelt/cli.py tests/test_cli.py
git commit -m "feat: validate chromosome name overlap in element, custom, metaplot"
```

---

## Phase 6: Documentation (Tasks 9-10)

### Task 9: Update design document

**Files:**
- Modify: `docs/superpowers/specs/2026-04-29-smelt-design.md`

- [ ] **Step 1: Update element section**

Remove "Merge consecutive rows belonging to the same feature." Replace with:
"Each GTF row is processed independently. Row-level statistics preserve the original GTF structure."

- [ ] **Step 2: Update input format section**

Replace cov.gz references with SAM/BAM and CX_report.txt. Note that strand and context come from the input files.

- [ ] **Step 3: Commit**

```bash
git add docs/superpowers/specs/2026-04-29-smelt-design.md
git commit -m "docs: update design spec to reflect input format and element changes"
```

---

### Task 10: Write README.md

**Files:**
- Rewrite: `README.md`

- [ ] **Step 1: Write README**

```markdown
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

| Format | Source | Notes |
|--------|--------|-------|
| SAM/BAM | BISMARK alignment | XM tag provides context and methylation calls |
| CX_report.txt | `coverage2cytosine` | Complete: chr, pos, strand, meth, unmeth, context |

## Subcommands

| Command | Purpose |
|---------|---------|
| `smelt site` | Per-cytosine methylation with context (CpG/CHH/CHG) |
| `smelt window` | Fixed-window sliding-window aggregation |
| `smelt element` | GTF feature-type methylation |
| `smelt metaplot` | Gene body +/- flanking equal-ratio binning |
| `smelt custom` | BED custom interval methylation |
| `smelt dmr` | Window-level Fisher exact + BH correction |
| `smelt stats` | Genome/chromosome summary by context |

## Output

All outputs are tab-separated TSV files. See `docs/superpowers/specs/2026-04-29-smelt-design.md` for full column specifications.
```

- [ ] **Step 2: Commit**

```bash
git add README.md
git commit -m "docs: write complete README with installation, usage, and subcommand overview"
```

---

## Phase 7: Final Integration (Task 11)

### Task 11: Full integration test fixup

**Files:**
- Modify: `scripts/generate_test_data.py` — generate CX_report format instead of cov
- Modify: `tests/test_integration.py` — use CX_report input

- [ ] **Step 1: Update data generator to produce CX_report format**

Change `generate_cov_files` to produce CX_report.txt files (with strand and context columns).

- [ ] **Step 2: Run full test suite**

```bash
pytest tests/ -v
```
Expected: all tests PASS

- [ ] **Step 3: Commit**

```bash
git add scripts/generate_test_data.py tests/test_integration.py tests/data/integration/
git commit -m "test: update integration tests for CX_report input format"
```

---
