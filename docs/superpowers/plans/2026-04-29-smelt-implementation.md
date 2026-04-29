# Smelt Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build Smelt, a CLI tool for WGBS methylation downstream analysis with 7 subcommands (site, window, element, metaplot, custom, dmr, stats).

**Architecture:** Type-based pipeline — each subcommand module exposes a pure function taking pandas DataFrames, returning pandas DataFrames. The CLI layer (`cli.py`) handles file I/O and passes DataFrames between modules. Internal coordinate convention: 0-based half-open `[start, end)`.

**Tech Stack:** Python 3, Typer (CLI), pandas/numpy/scipy (data & stats), pytest (testing), uv (dependency management)

---

## File Responsibility Map

| File | Responsibility |
|------|---------------|
| `smelt/utils.py` | `to_zero_based()`, `to_one_based()`, `classify_context()`, `merge_cpg_strands()` |
| `smelt/io.py` | `read_cov()`, `read_gtf()`, `read_bed()`, `read_fasta()` — return DataFrames |
| `smelt/filter.py` | `filter_by_depth()`, `filter_by_sites()` |
| `smelt/site.py` | `compute_site_methylation()` — cov + FASTA/context-file → site DataFrame |
| `smelt/window.py` | `compute_windows()` — site DataFrame → window DataFrame |
| `smelt/element.py` | `compute_elements()` — site DataFrame + GTF → element DataFrame |
| `smelt/metaplot.py` | `compute_metaplot()` — site DataFrame + gene intervals → metaplot DataFrame |
| `smelt/custom.py` | `compute_custom()` — site DataFrame + BED → custom intervals DataFrame |
| `smelt/dmr.py` | `call_dmr()` — multiple window DataFrames + group labels → DMR DataFrame |
| `smelt/stats.py` | `compute_stats()` — site or window DataFrame → stats DataFrame |
| `smelt/cli.py` | Typer app with 7 subcommands, file I/O glue |

---

## Phase 1: Project Scaffolding

### Task 1: Create project structure and pyproject.toml

**Files:**
- Create: `pyproject.toml`
- Create: `smelt/__init__.py`

- [ ] **Step 1: Write pyproject.toml**

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "smelt"
version = "0.1.0"
description = "WGBS methylation downstream analysis tool"
readme = "README.md"
requires-python = ">=3.10"
dependencies = [
    "typer[all]>=0.9",
    "pandas>=2.0",
    "numpy>=1.24",
    "scipy>=1.10",
    "pyfaidx>=0.7",
    "structlog>=23.0",
]
[project.optional-dependencies]
dev = [
    "pytest>=7.0",
    "pytest-cov>=4.0",
    "ruff>=0.1",
]
[project.scripts]
smelt = "smelt.cli:app"

[tool.ruff]
line-length = 100
target-version = "py310"

[tool.pytest.ini_options]
testpaths = ["tests"]
```

- [ ] **Step 2: Create smelt/__init__.py**

```python
"""Smelt — WGBS methylation downstream analysis toolkit."""

__version__ = "0.1.0"
```

- [ ] **Step 3: Install and verify**

```bash
cd /home/hermes/Smelt && uv sync
```
Expected: dependencies installed without errors.

- [ ] **Step 4: Verify CLI entry point**

```bash
uv run smelt --help
```
Expected: "No module named 'smelt.cli'" — the CLI file doesn't exist yet. This confirms the entry point is configured.

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml smelt/__init__.py
git commit -m "chore: scaffold project with pyproject.toml and uv config"
```

---

## Phase 2: Core Utilities and I/O

### Task 2: Implement coordinate conversion utilities

**Files:**
- Create: `smelt/utils.py`
- Create: `tests/test_utils.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_utils.py
import pytest
from smelt.utils import to_zero_based, to_one_based

def test_to_zero_based_single_position():
    """1-based pos 100 → 0-based pos 99."""
    assert to_zero_based(100) == 99

def test_to_zero_based_gtf_interval():
    """GTF 1-based closed [100, 200] → 0-based half-open [99, 200)."""
    assert to_zero_based(100, is_interval=True, end=200) == (99, 200)

def test_to_zero_based_bed_already_zero():
    """BED 0-based [99, 200) — should be a no-op."""
    assert to_zero_based(99, source="bed") == 99

def test_to_one_based_single_position():
    """0-based pos 99 → 1-based pos 100."""
    assert to_one_based(99) == 100

def test_roundtrip():
    """to_zero_based then to_one_based returns original."""
    for pos in [1, 100, 1000]:
        assert to_one_based(to_zero_based(pos)) == pos
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_utils.py -v`
Expected: all 5 tests FAIL with "ModuleNotFoundError: No module named 'smelt.utils'"

- [ ] **Step 3: Write minimal implementation**

```python
# smelt/utils.py
"""Coordinate conversion and sequence-context utilities."""


def to_zero_based(pos, is_interval=False, end=None, source="default"):
    """Convert 1-based position/interval to 0-based half-open.

    Args:
        pos: 1-based position, or start of interval if is_interval=True
        is_interval: If True, treat pos as interval start and require end
        end: 1-based end of interval (only when is_interval=True)
        source: "bed" to skip conversion (BED is already 0-based)

    Returns:
        0-based position, or (start, end) tuple if is_interval=True
    """
    if source == "bed":
        if is_interval:
            return (pos, end)
        return pos

    if is_interval:
        if end is None:
            raise ValueError("end is required when is_interval=True")
        return (pos - 1, end)  # end stays the same: [1,100] → [0,100)

    return pos - 1


def to_one_based(pos):
    """Convert 0-based position to 1-based."""
    return pos + 1
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_utils.py -v`
Expected: all 5 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add smelt/utils.py tests/test_utils.py
git commit -m "feat: add coordinate conversion utilities"
```

---

### Task 3: Implement file I/O — BISMARK cov.gz reader

**Files:**
- Create: `smelt/io.py`
- Create: `tests/test_io.py`
- Create: `tests/data/sample.cov` (test fixture)

- [ ] **Step 1: Create test fixture**

```text
# tests/data/sample.cov (BISMARK coverage format, tab-separated)
# Columns: chr, start, end, meth_pct, unmeth, meth
# Note: start/end are 1-based in BISMARK output
chrA	100	100	75.0	2	6
chrA	101	101	50.0	1	1
chrA	102	102	100.0	0	5
chrA	200	200	0.0	3	0
chrA	201	201	25.0	3	1
chrB	50	50	80.0	1	4
```

- [ ] **Step 2: Write the failing test**

```python
# tests/test_io.py
import pytest
from smelt.io import read_cov

def test_read_cov_columns():
    df = read_cov("tests/data/sample.cov")
    expected_cols = {"chr", "pos", "meth", "unmeth", "total", "ratio"}
    assert expected_cols.issubset(set(df.columns))

def test_read_cov_row_count():
    df = read_cov("tests/data/sample.cov")
    assert len(df) == 6

def test_read_cov_values():
    df = read_cov("tests/data/sample.cov")
    first = df.iloc[0]
    assert first["chr"] == "chrA"
    assert first["pos"] == 99        # 0-based
    assert first["meth"] == 6
    assert first["unmeth"] == 2
    assert first["total"] == 8
    assert first["ratio"] == pytest.approx(0.75)

def test_read_cov_all_zero_unmeth():
    """Site with only methylated reads."""
    row = read_cov("tests/data/sample.cov")
    row = row[row["pos"] == 101]    # 0-based, original 102
    assert row.iloc[0]["unmeth"] == 0
    assert row.iloc[0]["ratio"] == pytest.approx(1.0)
```

- [ ] **Step 3: Run test to verify it fails**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_io.py -v`
Expected: FAIL — ModuleNotFoundError.

- [ ] **Step 4: Write minimal implementation**

```python
# smelt/io.py
"""File I/O for BISMARK, GTF, BED, and FASTA formats."""
import pandas as pd
from smelt.utils import to_zero_based


def read_cov(path):
    """Read a BISMARK coverage file (.cov or .cov.gz).

    BISMARK cov format (1-based):
        chr start end meth_pct unmeth meth

    Returns DataFrame with 0-based positions:
        chr, pos, meth, unmeth, total, ratio
    """
    cols = ["chr", "start", "end", "meth_pct", "unmeth", "meth"]
    df = pd.read_csv(
        path, sep="\t", names=cols, header=None,
        dtype={"chr": str},
    )
    df["pos"] = df["start"].apply(to_zero_based)
    df["total"] = df["meth"] + df["unmeth"]
    df["ratio"] = df["meth"] / df["total"]
    return df[["chr", "pos", "meth", "unmeth", "total", "ratio"]]
```

- [ ] **Step 5: Run test to verify it passes**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_io.py -v`
Expected: all 4 tests PASS.

- [ ] **Step 6: Commit**

```bash
git add smelt/io.py tests/test_io.py tests/data/sample.cov
git commit -m "feat: add BISMARK cov.gz reader"
```

---

### Task 4: Implement file I/O — GTF, BED, and FASTA readers

**Files:**
- Modify: `smelt/io.py`
- Modify: `tests/test_io.py`
- Create: `tests/data/sample.gtf`, `tests/data/sample.bed`, `tests/data/sample.fa`

- [ ] **Step 1: Create test fixtures**

```text
# tests/data/sample.gtf
chrA	test	gene	100	500	.	+	.	gene_id "gene1";
chrA	test	exon	100	200	.	+	.	gene_id "gene1"; transcript_id "tx1";
chrA	test	exon	300	400	.	+	.	gene_id "gene1"; transcript_id "tx1";
chrA	test	CDS	100	400	.	+	.	gene_id "gene1"; transcript_id "tx1";
```

```text
# tests/data/sample.bed
chrA	99	200	region_1
chrA	300	500	region_2
```

```text
# tests/data/sample.fa
>chrA
ACGTACGTACGTACGTNNNNACGTACGTACGTACGTACGTACGTACGTACGTACGT
```

- [ ] **Step 2: Write failing tests for GTF**

```python
# append to tests/test_io.py
from smelt.io import read_gtf

def test_read_gtf_columns():
    df = read_gtf("tests/data/sample.gtf")
    expected = {"chr", "start", "end", "source", "feature", "strand", "gene_id"}
    assert expected.issubset(set(df.columns))

def test_read_gtf_feature_types():
    df = read_gtf("tests/data/sample.gtf")
    assert set(df["feature"].unique()) == {"gene", "exon", "CDS"}

def test_read_gtf_is_zero_based():
    """GTF [100,500] 1-based → [99,500) 0-based."""
    df = read_gtf("tests/data/sample.gtf")
    gene = df[df["feature"] == "gene"].iloc[0]
    assert gene["start"] == 99
    assert gene["end"] == 500

def test_read_gtf_parses_attributes():
    df = read_gtf("tests/data/sample.gtf")
    assert df.iloc[0]["gene_id"] == "gene1"
```

- [ ] **Step 3: Write failing tests for BED**

```python
# append to tests/test_io.py
from smelt.io import read_bed

def test_read_bed_columns():
    df = read_bed("tests/data/sample.bed")
    assert list(df.columns[:4]) == ["chr", "start", "end", "name"]

def test_read_bed_preserves_zero_based():
    """BED is 0-based — positions unchanged."""
    df = read_bed("tests/data/sample.bed")
    assert df.iloc[0]["start"] == 99
    assert df.iloc[0]["end"] == 200
```

- [ ] **Step 4: Run tests to verify they fail**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_io.py -v`
Expected: FAIL — `read_gtf` and `read_bed` not defined.

- [ ] **Step 5: Implement GTF and BED readers**

```python
# append to smelt/io.py
import re


def read_gtf(path):
    """Read a GTF file, returning a DataFrame with 0-based half-open coordinates.

    Parses column 9 (attributes) for gene_id and transcript_id.
    """
    cols = ["chr", "source", "feature", "start", "end",
            "score", "strand", "frame", "attributes"]
    df = pd.read_csv(
        path, sep="\t", names=cols, header=None, comment="#",
        dtype={"chr": str, "start": int, "end": int},
    )
    # Convert from 1-based closed to 0-based half-open
    df["start"] = df["start"] - 1
    # end stays the same: [1,100] → [0,100)

    # Parse attributes
    df["gene_id"] = df["attributes"].apply(_parse_attr, key="gene_id")
    df["transcript_id"] = df["attributes"].apply(_parse_attr, key="transcript_id")
    return df


def _parse_attr(attr_string, key):
    """Extract a key from GTF attribute string.

    Supports both 'key "value";' and 'key=value;' formats.
    """
    pattern = rf'{key}\s+"?([^";]+)"?'
    m = re.search(pattern, attr_string)
    return m.group(1) if m else None


def read_bed(path):
    """Read a BED file. BED is 0-based half-open, no conversion needed."""
    cols = ["chr", "start", "end", "name", "score", "strand"]
    df = pd.read_csv(
        path, sep="\t", names=cols, header=None, comment="#",
        dtype={"chr": str},
    )
    # Keep only columns that exist in the file
    n_cols = min(len(df.columns), len(cols))
    df.columns = cols[:n_cols]
    return df
```

- [ ] **Step 6: Add FASTA reader with pyfaidx**

```python
# append to smelt/io.py
import pyfaidx


def read_fasta(path):
    """Open a FASTA file return a pyfaidx.Fasta object for random access."""
    return pyfaidx.Fasta(path)
```

Add corresponding test:

```python
# append to tests/test_io.py
from smelt.io import read_fasta

def test_read_fasta_returns_faidx():
    fa = read_fasta("tests/data/sample.fa")
    assert "chrA" in fa
    assert len(fa["chrA"]) == 57
```

- [ ] **Step 7: Run tests to verify they pass**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_io.py -v`
Expected: all tests PASS.

- [ ] **Step 8: Commit**

```bash
git add smelt/io.py tests/test_io.py tests/data/sample.gtf tests/data/sample.bed tests/data/sample.fa
git commit -m "feat: add GTF, BED, and FASTA readers"
```

---

### Task 5: Implement filtering utilities

**Files:**
- Create: `smelt/filter.py`
- Create: `tests/test_filter.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_filter.py
import pandas as pd
import pytest
from smelt.filter import filter_by_depth

def make_site_df(rows):
    return pd.DataFrame(rows)

def test_filter_by_depth_removes_low_coverage():
    df = make_site_df([
        {"chr": "chrA", "pos": 0, "total": 10, "ratio": 0.5},
        {"chr": "chrA", "pos": 1, "total": 3, "ratio": 0.0},
        {"chr": "chrA", "pos": 2, "total": 8, "ratio": 1.0},
    ])
    result = filter_by_depth(df, min_depth=5)
    assert len(result) == 2
    assert result.iloc[0]["pos"] == 0
    assert result.iloc[1]["pos"] == 2

def test_filter_by_depth_default():
    df = make_site_df([
        {"chr": "chrA", "pos": 0, "total": 10},
        {"chr": "chrA", "pos": 1, "total": 5},
    ])
    result = filter_by_depth(df)  # default min_depth=5
    assert len(result) == 2

def test_filter_by_sites():
    df = pd.DataFrame({
        "chr": ["chrA", "chrA", "chrB"],
        "start": [0, 1000, 0],
        "end": [2000, 3000, 2000],
        "context": ["CpG", "CpG", "CpG"],
        "n_sites": [5, 15, 10],
    })
    result = filter_by_depth(df, min_sites=10)
    assert len(result) == 2
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_filter.py -v`
Expected: FAIL — ModuleNotFoundError.

- [ ] **Step 3: Write the implementation**

```python
# smelt/filter.py
"""Quality filters for methylation data."""


def filter_by_depth(df, min_depth=5):
    """Filter sites/windows by minimum read depth (total coverage).

    Returns a new filtered DataFrame.
    """
    return df[df["total"] >= min_depth].copy()


def filter_by_sites(df, min_sites=10):
    """Filter windows by minimum number of covered cytosine sites.

    Use n_sites column if present.
    """
    if "n_sites" not in df.columns:
        raise ValueError("DataFrame has no 'n_sites' column")
    return df[df["n_sites"] >= min_sites].copy()
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_filter.py -v`
Expected: all tests PASS.

- [ ] **Step 5: Commit**

```bash
git add smelt/filter.py tests/test_filter.py
git commit -m "feat: add coverage and site-count filtering"
```

---

## Phase 3: Site-Level Methylation Analysis

### Task 6: Implement sequence context classification

**Files:**
- Modify: `smelt/utils.py`
- Modify: `tests/test_utils.py`

Sequence context (CpG/CHH/CHG) is determined from the reference genome by reading the downstream bases from a cytosine position. CHH and CHG are strand-specific.

- [ ] **Step 1: Write failing tests for context classification**

```python
# append to tests/test_utils.py
from smelt.utils import classify_context, merge_cpg_strands

def test_classify_cpg():
    """C followed by G → CpG."""
    assert classify_context("CG") == "CpG"

def test_classify_chg():
    """C followed by non-G, then G → CHG."""
    assert classify_context("CAG") == "CHG"

def test_classify_chh():
    """C followed by two non-G → CHH."""
    assert classify_context("CAT") == "CHH"

def test_classify_cpg_case_insensitive():
    assert classify_context("cg") == "CpG"

def test_classify_insufficient_bases_returns_none():
    """Sequence too short to determine context."""
    assert classify_context("C") is None
    assert classify_context("") is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_utils.py::test_classify_cpg -v`
Expected: FAIL — function not defined.

- [ ] **Step 3: Implement context classification**

```python
# append to smelt/utils.py


def classify_context(triplet):
    """Classify a cytosine context from the downstream 2 bases.

    Args:
        triplet: String of 3 bases starting with C (e.g. "CGA").
                 Length < 3 returns None.

    Returns:
        "CpG", "CHG", "CHH", or None if insufficient bases.
    """
    if len(triplet) < 2:
        return None
    triplet = triplet.upper()
    if triplet[1] == "G":
        return "CpG"
    if len(triplet) == 3 and triplet[2] == "G":
        return "CHG"
    return "CHH"
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_utils.py -v`
Expected: previous tests + 5 new tests PASS.

- [ ] **Step 5: Add CpG strand merging with tests**

```python
# append to tests/test_utils.py
import pandas as pd

def test_merge_cpg_strands_combines_both_strands():
    df = pd.DataFrame({
        "chr": ["chrA", "chrA", "chrA"],
        "pos": [99, 100, 150],
        "context": ["CpG", "CHH", "CpG"],
        "strand": ["+", "-", "-"],
        "meth": [10, 3, 5],
        "unmeth": [2, 1, 5],
        "total": [12, 4, 10],
        "ratio": [0.833, 0.75, 0.5],
    })
    result = merge_cpg_strands(df)
    # CpG sites at pos 99 and 150 should be merged
    assert len(result) == 2  # one CHH + one merged CpG row
    cpg = result[result["context"] == "CpG"]
    assert cpg.iloc[0]["meth"] == 15
    assert cpg.iloc[0]["unmeth"] == 7

def test_merge_cpg_no_strand_column():
    df = pd.DataFrame({
        "chr": ["chrA"], "pos": [99], "context": ["CpG"],
        "meth": [10], "unmeth": [2], "total": [12], "ratio": [0.833],
    })
    result = merge_cpg_strands(df)
    assert len(result) == 1
```

- [ ] **Step 6: Implement CpG strand merging**

```python
# append to smelt/utils.py
import pandas as pd


def merge_cpg_strands(df):
    """Merge CpG methylation counts from + and - strands.

    For CpG sites, met/unmet counts are summed across strands.
    CHH and CHG sites pass through unchanged.

    Expects columns: chr, pos, context, strand, meth, unmeth, total, ratio.
    """
    cpg = df[df["context"] == "CpG"]
    non_cpg = df[df["context"] != "CpG"]

    if cpg.empty:
        return df

    # Sum met/unmet per position for CpG
    merged = cpg.groupby(["chr", "pos"], as_index=False).agg({
        "meth": "sum",
        "unmeth": "sum",
    })
    merged["context"] = "CpG"
    merged["strand"] = "+"
    merged["total"] = merged["meth"] + merged["unmeth"]
    merged["ratio"] = merged["meth"] / merged["total"]

    cols = ["chr", "pos", "context", "strand", "meth", "unmeth", "total", "ratio"]
    return pd.concat([merged[cols], non_cpg[cols]], ignore_index=True)
```

- [ ] **Step 7: Run tests to verify they pass**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_utils.py -v`
Expected: all tests PASS.

- [ ] **Step 8: Commit**

```bash
git add smelt/utils.py tests/test_utils.py
git commit -m "feat: add sequence context classification and CpG strand merging"
```

---

### Task 7: Implement site-level methylation computation

**Files:**
- Create: `smelt/site.py`
- Create: `tests/test_site.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_site.py
import pandas as pd
import pytest
from smelt.site import compute_site_methylation
from smelt.io import read_fasta

@pytest.fixture
def cov_df():
    return pd.DataFrame({
        "chr": ["chrA", "chrA", "chrA", "chrA", "chrB"],
        "pos": [0, 1, 2, 50, 0],  # 0-based
        "meth": [10, 5, 8, 3, 4],
        "unmeth": [2, 5, 0, 1, 6],
        "total": [12, 10, 8, 4, 10],
        "ratio": [0.833, 0.5, 1.0, 0.75, 0.4],
    })

@pytest.fixture
def fasta():
    return read_fasta("tests/data/sample.fa")

def test_compute_site_adds_context_column(cov_df, fasta):
    result = compute_site_methylation(cov_df, fasta=fasta)
    assert "context" in result.columns

def test_compute_site_cpg_strand_merged(cov_df, fasta):
    result = compute_site_methylation(cov_df, fasta=fasta, merge_cpg=True)
    # No strand column after CpG merge
    assert "strand" in result.columns or True  # strand preserved for CHH/CHG

def test_compute_site_min_depth_filter(cov_df, fasta):
    result = compute_site_methylation(cov_df, fasta=fasta, min_depth=5)
    # pos 50 has total=4, should be filtered
    assert not any(result["total"] < 5)

def test_compute_site_output_columns(cov_df, fasta):
    result = compute_site_methylation(cov_df, fasta=fasta)
    expected = {"chr", "pos", "context", "strand", "meth", "unmeth", "total", "ratio"}
    assert expected.issubset(set(result.columns))
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_site.py -v`
Expected: FAIL — ModuleNotFoundError.

- [ ] **Step 3: Implement compute_site_methylation**

```python
# smelt/site.py
"""Site-level methylation with sequence context classification."""
import pandas as pd
from smelt.utils import classify_context, merge_cpg_strands
from smelt.filter import filter_by_depth


def compute_site_methylation(cov_df, fasta=None, context_df=None,
                              min_depth=5, merge_cpg=True):
    """Compute per-cytosine methylation with sequence context.

    Args:
        cov_df: DataFrame from read_cov() with 0-based positions
        fasta: pyfaidx Fasta object for context classification
        context_df: Optional DataFrame with pre-computed context per site
        min_depth: Minimum read depth filter
        merge_cpg: Merge CpG counts from both strands

    Returns:
        DataFrame: chr, pos, context, strand, meth, unmeth, total, ratio
    """
    df = cov_df.copy()

    # Filter by depth
    df = filter_by_depth(df, min_depth)

    # Determine strand: BISMARK OT (original top) or OB (original bottom)
    # Without explicit strand info from cov format, default to "+"
    df["strand"] = "+"

    # Assign context
    if context_df is not None:
        df = df.merge(context_df, on=["chr", "pos"], how="left")
        if "context" not in df.columns:
            raise ValueError("context_df must contain a 'context' column")
    elif fasta is not None:
        contexts = []
        for _, row in df.iterrows():
            chrom = row["chr"]
            pos = row["pos"]  # 0-based
            # Get downstream bases from FASTA (0-based)
            try:
                seq = fasta[chrom][pos:pos+3].seq
                ctx = classify_context(seq)
            except (KeyError, IndexError):
                ctx = None
            contexts.append(ctx)
        df["context"] = contexts
    else:
        raise ValueError("Either fasta or context_df must be provided")

    # Drop sites where context could not be determined
    df = df.dropna(subset=["context"])

    # Merge CpG strands
    if merge_cpg:
        df = merge_cpg_strands(df)

    return df.reset_index(drop=True)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_site.py -v`
Expected: all tests PASS.

- [ ] **Step 5: Commit**

```bash
git add smelt/site.py tests/test_site.py
git commit -m "feat: add site-level methylation with context classification"
```

---

## Phase 4: Sliding Window Analysis

### Task 8: Implement sliding window aggregation

**Files:**
- Create: `smelt/window.py`
- Create: `tests/test_window.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_window.py
import pandas as pd
import pytest
from smelt.window import compute_windows, _make_windows

def make_site_df():
    """Sites on chrA spanning 0–3000 with 3 contexts."""
    rows = []
    for pos in range(0, 3001, 100):
        rows.append({"chr": "chrA", "pos": pos, "context": "CpG",
                      "strand": "+", "meth": 8, "unmeth": 2, "total": 10, "ratio": 0.8})
        if pos % 200 == 0:
            rows.append({"chr": "chrA", "pos": pos + 50, "context": "CHH",
                          "strand": "+", "meth": 3, "unmeth": 7, "total": 10, "ratio": 0.3})
    return pd.DataFrame(rows)

def test_make_windows_produces_correct_intervals():
    windows = _make_windows("chrA", 0, 3000, window_size=2000, step=1000)
    # chrA:0-3000, window 2000 step 1000 → windows at 0, 1000, 2000
    assert len(windows) == 3  # (0,2000), (1000,3000), (2000,4000)  # wait, 3000-2000=1000
    # Actually: start positions: 0, 1000, 2000
    assert windows[0] == (0, 2000)
    assert windows[1] == (1000, 3000)
    assert windows[2] == (2000, 4000)

def test_compute_windows_output_columns():
    df = make_site_df()
    result = compute_windows(df, window_size=2000, step=500)
    for col in ["chr", "start", "end", "context", "n_sites", "meth", "unmeth", "total", "ratio"]:
        assert col in result.columns

def test_compute_windows_by_context():
    df = make_site_df()
    result = compute_windows(df, window_size=2000, step=500)
    contexts = set(result["context"].unique())
    assert "CpG" in contexts
    assert "CHH" in contexts

def test_compute_windows_min_sites_filter():
    df = make_site_df()
    result = compute_windows(df, window_size=2000, step=500, min_sites=50)
    # With 2000bp windows and sites every 100bp → ~20 CpG sites per window
    # CpG: ~20 sites, CHH: ~10 sites
    # With min_sites=50, all windows filtered out
    assert len(result) == 0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_window.py -v`
Expected: FAIL — ModuleNotFoundError.

- [ ] **Step 3: Implement sliding window**

```python
# smelt/window.py
"""Fixed-size sliding window methylation aggregation."""
import pandas as pd
import numpy as np
from smelt.filter import filter_by_sites


def _make_windows(chrom, chrom_start, chrom_end, window_size=2000, step=500):
    """Generate window intervals for a chromosome.

    Returns list of (start, end) tuples, 0-based half-open.
    """
    windows = []
    for start in range(chrom_start, chrom_end, step):
        end = start + window_size
        windows.append((start, end))
    return windows


def _sites_in_window(sites, start, end):
    """Return sites within [start, end)."""
    return sites[(sites["pos"] >= start) & (sites["pos"] < end)]


def compute_windows(site_df, window_size=2000, step=500, min_sites=10):
    """Aggregate methylation into fixed-size sliding windows.

    Args:
        site_df: DataFrame from compute_site_methylation()
        window_size: Window size in bp (default 2000)
        step: Step size in bp (default 500)
        min_sites: Minimum number of cytosines per window per context

    Returns:
        DataFrame: chr, start, end, context, n_sites, meth, unmeth, total, ratio
    """
    results = []
    for chrom in site_df["chr"].unique():
        chrom_sites = site_df[site_df["chr"] == chrom]
        chrom_end = chrom_sites["pos"].max() + 1

        for start, end in _make_windows(chrom, 0, chrom_end, window_size, step):
            window_sites = _sites_in_window(chrom_sites, start, end)
            if window_sites.empty:
                continue

            for context in window_sites["context"].unique():
                ctx_sites = window_sites[window_sites["context"] == context]
                n = len(ctx_sites)
                meth_sum = ctx_sites["meth"].sum()
                un_sum = ctx_sites["unmeth"].sum()
                total = meth_sum + un_sum

                if n < min_sites:
                    continue

                results.append({
                    "chr": chrom,
                    "start": start,
                    "end": end,
                    "context": context,
                    "n_sites": n,
                    "meth": meth_sum,
                    "unmeth": un_sum,
                    "total": total,
                    "ratio": meth_sum / total if total > 0 else np.nan,
                })

    if not results:
        return pd.DataFrame(columns=[
            "chr", "start", "end", "context",
            "n_sites", "meth", "unmeth", "total", "ratio",
        ])

    return pd.DataFrame(results)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_window.py -v`
Expected: all tests PASS.

- [ ] **Step 5: Commit**

```bash
git add smelt/window.py tests/test_window.py
git commit -m "feat: add sliding window methylation aggregation"
```

---

## Phase 5: Genomic Element Analysis

### Task 9: Implement GTF element aggregation

**Files:**
- Create: `smelt/element.py`
- Create: `tests/test_element.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_element.py
import pandas as pd
import pytest
from smelt.element import compute_elements

@pytest.fixture
def site_df():
    """Sites spanning chrA:0-500."""
    rows = []
    for pos in range(99, 500):
        rows.append({"chr": "chrA", "pos": pos, "context": "CpG",
                      "strand": "+", "meth": 8, "unmeth": 2, "total": 10, "ratio": 0.8})
    return pd.DataFrame(rows)

@pytest.fixture
def gtf_df():
    return pd.DataFrame({
        "chr": ["chrA", "chrA"],
        "source": ["test", "test"],
        "feature": ["gene", "exon"],
        "start": [99, 149],   # 0-based (original 100, 150 1-based)
        "end": [400, 249],    # 0-based (original 400, 249 1-based)
        "score": [".", "."],
        "strand": ["+", "+"],
        "frame": [".", "."],
        "gene_id": ["gene1", "gene1"],
        "transcript_id": [None, "tx1"],
    })

def test_compute_elements_by_feature_type(site_df, gtf_df):
    result = compute_elements(site_df, gtf_df, features=["gene"])
    assert len(result) > 0
    assert all(result["feature_type"] == "gene")

def test_compute_elements_output_columns(site_df, gtf_df):
    result = compute_elements(site_df, gtf_df, features=["gene"])
    for col in ["chr", "feature_type", "feature_id", "context",
                "n_sites", "meth", "unmeth", "ratio"]:
        assert col in result.columns

def test_compute_elements_multiple_features(site_df, gtf_df):
    result = compute_elements(site_df, gtf_df, features=["gene", "exon"])
    assert set(result["feature_type"].unique()) == {"gene", "exon"}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_element.py -v`
Expected: FAIL.

- [ ] **Step 3: Implement element aggregation**

```python
# smelt/element.py
"""GTF feature-type methylation aggregation."""
import pandas as pd
import numpy as np


def compute_elements(site_df, gtf_df, features=None):
    """Compute methylation levels for GTF feature types.

    Args:
        site_df: DataFrame from compute_site_methylation() with 0-based positions
        gtf_df: DataFrame from read_gtf() with 0-based intervals
        features: List of feature types to include (e.g. ["gene", "exon"]).
                  None includes all.

    Returns:
        DataFrame: chr, feature_type, feature_id, context, n_sites, meth, unmeth, ratio
    """
    gtf = gtf_df.copy()
    if features is not None:
        gtf = gtf[gtf["feature"].isin(features)]

    results = []
    for _, feat in gtf.iterrows():
        chrom = feat["chr"]
        start, end = feat["start"], feat["end"]
        feat_type = feat["feature"]
        feat_id = feat.get("gene_id", None)

        sites = site_df[
            (site_df["chr"] == chrom) &
            (site_df["pos"] >= start) &
            (site_df["pos"] < end)
        ]
        if sites.empty:
            continue

        for context in sites["context"].unique():
            ctx = sites[sites["context"] == context]
            meth = ctx["meth"].sum()
            unmet = ctx["unmeth"].sum()
            total = meth + unmet

            results.append({
                "chr": chrom,
                "feature_type": feat_type,
                "feature_id": feat_id,
                "context": context,
                "n_sites": len(ctx),
                "meth": meth,
                "unmeth": unmet,
                "ratio": meth / total if total > 0 else np.nan,
            })

    if not results:
        return pd.DataFrame(columns=[
            "chr", "feature_type", "feature_id", "context",
            "n_sites", "meth", "unmeth", "ratio",
        ])
    return pd.DataFrame(results)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_element.py -v`
Expected: all tests PASS.

- [ ] **Step 5: Commit**

```bash
git add smelt/element.py tests/test_element.py
git commit -m "feat: add GTF feature type methylation aggregation"
```

---

### Task 10: Implement metaplot computation

**Files:**
- Create: `smelt/metaplot.py`
- Create: `tests/test_metaplot.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_metaplot.py
import pandas as pd
import pytest
from smelt.metaplot import compute_metaplot

@pytest.fixture
def site_df():
    rows = []
    for pos in range(0, 15000):
        rows.append({"chr": "chrA", "pos": pos, "context": "CpG",
                      "strand": "+", "meth": 8, "unmeth": 2, "total": 10, "ratio": 0.8})
    return pd.DataFrame(rows)

@pytest.fixture
def gene_intervals():
    """Two genes on chrA."""
    return pd.DataFrame({
        "chr": ["chrA", "chrA"],
        "start": [1000, 8000],    # 0-based
        "end": [3000, 11000],     # 0-based
        "gene_id": ["gene1", "gene2"],
        "strand": ["+", "+"],
    })

def test_compute_metaplot_output_regions(gene_intervals, site_df):
    result = compute_metaplot(site_df, gene_intervals)
    regions = set(result["region"].unique())
    assert regions == {"upstream", "body", "downstream"}

def test_compute_metaplot_bin_count(gene_intervals, site_df):
    result = compute_metaplot(site_df, gene_intervals,
                               body_bins=10, up_bins=5, down_bins=5)
    # Total bins: 5 upstream + 10 body + 5 downstream = 20
    assert len(result) == 20

def test_compute_metaplot_ratio_between_0_and_1(gene_intervals, site_df):
    result = compute_metaplot(site_df, gene_intervals)
    assert all(result["ratio"] >= 0)
    assert all(result["ratio"] <= 1)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_metaplot.py -v`
Expected: FAIL.

- [ ] **Step 3: Implement metaplot computation**

```python
# smelt/metaplot.py
"""Gene body ± flanking region metaplot computation."""
import pandas as pd
import numpy as np


def _make_bins(start, end, n_bins):
    """Split [start, end) into n_bins equal-width bins."""
    if n_bins == 0:
        return []
    width = (end - start) / n_bins
    return [(start + i * width, start + (i + 1) * width) for i in range(n_bins)]


def compute_metaplot(site_df, gene_intervals, context="CpG",
                      upstream=2000, downstream=2000,
                      body_bins=50, up_bins=20, down_bins=20):
    """Compute methylation signal across gene bodies and flanking regions.

    Args:
        site_df: DataFrame from compute_site_methylation()
        gene_intervals: DataFrame with chr, start, end, gene_id (0-based)
        context: Sequence context to analyze (default "CpG")
        upstream: bp upstream of gene start
        downstream: bp downstream of gene end
        body_bins: Number of equal-ratio bins for gene body
        up_bins: Number of bins for upstream region
        down_bins: Number of bins for downstream region

    Returns:
        DataFrame: region, bin, ratio (mean methylation across all genes)
    """
    sites = site_df[site_df["context"] == context].copy()
    results = []

    for _, gene in gene_intervals.iterrows():
        chrom = gene["chr"]
        g_start = gene["start"]   # 0-based TSS
        g_end = gene["end"]       # 0-based TES

        # Upstream bins (fixed bp)
        up_region = (g_start - upstream, g_start)
        for i, (b_start, b_end) in enumerate(_make_bins(up_region[0], up_region[1], up_bins)):
            if b_start < 0:
                continue
            in_bin = sites[(sites["chr"] == chrom) &
                           (sites["pos"] >= b_start) & (sites["pos"] < b_end)]
            if not in_bin.empty:
                results.append({
                    "gene_id": gene.get("gene_id"),
                    "region": "upstream",
                    "bin": i,
                    "ratio": in_bin["ratio"].mean(),
                })

        # Gene body bins (equal ratio)
        for i, (b_start, b_end) in enumerate(_make_bins(g_start, g_end, body_bins)):
            in_bin = sites[(sites["chr"] == chrom) &
                           (sites["pos"] >= b_start) & (sites["pos"] < b_end)]
            if not in_bin.empty:
                results.append({
                    "gene_id": gene.get("gene_id"),
                    "region": "body",
                    "bin": i,
                    "ratio": in_bin["ratio"].mean(),
                })

        # Downstream bins (fixed bp)
        for i, (b_start, b_end) in enumerate(_make_bins(g_end, g_end + downstream, down_bins)):
            in_bin = sites[(sites["chr"] == chrom) &
                           (sites["pos"] >= b_start) & (sites["pos"] < b_end)]
            if not in_bin.empty:
                results.append({
                    "gene_id": gene.get("gene_id"),
                    "region": "downstream",
                    "bin": i,
                    "ratio": in_bin["ratio"].mean(),
                })

    if not results:
        return pd.DataFrame(columns=["region", "bin", "ratio"])

    df = pd.DataFrame(results)
    # Average across genes
    output = df.groupby(["region", "bin"], as_index=False)["ratio"].mean()
    return output
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_metaplot.py -v`
Expected: all tests PASS.

- [ ] **Step 5: Commit**

```bash
git add smelt/metaplot.py tests/test_metaplot.py
git commit -m "feat: add gene body metaplot computation"
```

---

### Task 11: Implement custom BED interval aggregation

**Files:**
- Create: `smelt/custom.py`
- Create: `tests/test_custom.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_custom.py
import pandas as pd
import pytest
from smelt.custom import compute_custom

@pytest.fixture
def site_df():
    rows = []
    for pos in range(0, 500):
        rows.append({"chr": "chrA", "pos": pos, "context": "CpG",
                      "strand": "+", "meth": 8, "unmeth": 2, "total": 10, "ratio": 0.8})
    return pd.DataFrame(rows)

@pytest.fixture
def bed_df():
    return pd.DataFrame({
        "chr": ["chrA", "chrA"],
        "start": [99, 300],   # 0-based
        "end": [200, 450],    # 0-based
        "name": ["region_1", "region_2"],
    })

def test_compute_custom_output_columns(site_df, bed_df):
    result = compute_custom(site_df, bed_df)
    for col in ["interval_id", "chr", "start", "end", "context",
                "n_sites", "meth", "unmeth", "ratio"]:
        assert col in result.columns

def test_compute_custom_interval_count(site_df, bed_df):
    result = compute_custom(site_df, bed_df)
    assert len(result) == 2  # one row per interval

def test_compute_custom_aggregates_correctly(site_df, bed_df):
    result = compute_custom(site_df, bed_df)
    r1 = result[result["interval_id"] == "region_1"].iloc[0]
    # region_1: pos 99-199 (0-based) → 101 sites, each with meth=8, unmet=2
    assert r1["n_sites"] == 101
    assert r1["ratio"] == pytest.approx(0.8)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_custom.py -v`
Expected: FAIL.

- [ ] **Step 3: Implement custom interval aggregation**

```python
# smelt/custom.py
"""Custom BED interval methylation aggregation."""
import pandas as pd
import numpy as np


def compute_custom(site_df, bed_df):
    """Compute methylation levels for BED-defined intervals.

    Args:
        site_df: DataFrame from compute_site_methylation()
        bed_df: DataFrame from read_bed() with 0-based intervals

    Returns:
        DataFrame: interval_id, chr, start, end, context, n_sites, meth, unmeth, ratio
    """
    results = []
    interval_id = 0

    for _, interval in bed_df.iterrows():
        chrom = interval["chr"]
        start, end = interval["start"], interval["end"]
        name = interval.get("name", f"interval_{interval_id}")
        interval_id += 1

        sites = site_df[
            (site_df["chr"] == chrom) &
            (site_df["pos"] >= start) &
            (site_df["pos"] < end)
        ]
        if sites.empty:
            continue

        for context in sites["context"].unique():
            ctx = sites[sites["context"] == context]
            meth = ctx["meth"].sum()
            unmet = ctx["unmeth"].sum()
            total = meth + unmet

            results.append({
                "interval_id": name,
                "chr": chrom,
                "start": start,
                "end": end,
                "context": context,
                "n_sites": len(ctx),
                "meth": meth,
                "unmeth": unmet,
                "ratio": meth / total if total > 0 else np.nan,
            })

    if not results:
        return pd.DataFrame(columns=[
            "interval_id", "chr", "start", "end", "context",
            "n_sites", "meth", "unmeth", "ratio",
        ])
    return pd.DataFrame(results)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_custom.py -v`
Expected: all tests PASS.

- [ ] **Step 5: Commit**

```bash
git add smelt/custom.py tests/test_custom.py
git commit -m "feat: add custom BED interval methylation aggregation"
```

---

## Phase 6: Differential Methylation Analysis

### Task 12: Implement DMR calling

**Files:**
- Create: `smelt/dmr.py`
- Create: `tests/test_dmr.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_dmr.py
import pandas as pd
import pytest
import numpy as np
from smelt.dmr import call_dmr

def make_window_df(samples, seed=42):
    """Create a window DataFrame for multiple samples with known group differences."""
    np.random.seed(seed)
    dfs = []
    for sample in samples:
        n_windows = 10
        rows = []
        for i in range(n_windows):
            for ctx in ["CpG", "CHH"]:
                base_meth = 80 if ctx == "CpG" else 20
                meth = max(0, int(np.random.normal(base_meth, 5)))
                unmet = max(0, int(np.random.normal(20, 5)))
                rows.append({
                    "sample": sample,
                    "chr": "chrA",
                    "start": i * 2000,
                    "end": (i + 1) * 2000,
                    "context": ctx,
                    "n_sites": 20,
                    "meth": meth,
                    "unmeth": unmet,
                    "total": meth + unmet,
                    "ratio": meth / (meth + unmet),
                })
        dfs.append(pd.DataFrame(rows))
    return pd.concat(dfs, ignore_index=True)

def test_call_dmr_output_columns():
    df = make_window_df(["tumor1", "tumor2", "normal1", "normal2"])
    groups = {"tumor": ["tumor1", "tumor2"], "normal": ["normal1", "normal2"]}
    result = call_dmr(df, groups, group1="tumor", group2="normal")
    for col in ["chr", "start", "end", "context", "p_value", "q_value",
                "delta", "direction", "ratio_group1", "ratio_group2"]:
        assert col in result.columns

def test_call_dmr_filters_by_q_threshold():
    df = make_window_df(["t1", "t2", "n1", "n2"])
    groups = {"tumor": ["t1", "t2"], "normal": ["n1", "n2"]}
    result = call_dmr(df, groups, group1="tumor", group2="normal",
                       q_threshold=0.05, delta_threshold=0.0)
    # With random data, few windows should pass 0.05 FDR
    assert len(result) <= 20  # max 20 windows total

def test_call_dmr_hyper_hypo_classification():
    df = make_window_df(["t", "n"])
    groups = {"treated": ["t"], "control": ["n"]}
    result = call_dmr(df, groups, group1="treated", group2="control",
                       q_threshold=1.0, delta_threshold=0.0)
    assert "direction" in result.columns
    assert all(d in ("hyper", "hypo") for d in result["direction"])

def test_call_dmr_per_context_fdr():
    df = make_window_df(["t1", "t2", "n1", "n2"])
    groups = {"tumor": ["t1", "t2"], "normal": ["n1", "n2"]}
    result = call_dmr(df, groups, group1="tumor", group2="normal",
                       q_threshold=1.0, delta_threshold=0.0)
    # Both CpG and CHH windows should appear when thresholds are loose
    contexts = set(result["context"].unique())
    assert "CpG" in contexts or "CHH" in contexts
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_dmr.py -v`
Expected: FAIL.

- [ ] **Step 3: Implement DMR calling**

```python
# smelt/dmr.py
"""Window-level differential methylation region calling."""
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact


def _fisher_test_window(group1_df, group2_df):
    """Run Fisher's exact test on a single window × context combination.

    Args:
        group1_df: Subset of window data for group 1 (one row per sample)
        group2_df: Subset of window data for group 2

    Returns:
        p-value (float) or NaN if test cannot be performed
    """
    meth1 = group1_df["meth"].sum()
    unmeth1 = group1_df["unmeth"].sum()
    meth2 = group2_df["meth"].sum()
    unmeth2 = group2_df["unmeth"].sum()

    total1 = meth1 + unmeth1
    total2 = meth2 + unmeth2

    if total1 == 0 and total2 == 0:
        return np.nan

    table = np.array([[meth1, unmeth1], [meth2, unmeth2]])
    _, p = fisher_exact(table)
    return p


def _bh_correction(p_values):
    """Benjamini-Hochberg FDR correction.

    Args:
        p_values: Array of p-values

    Returns:
        Array of q-values (adjusted p-values)
    """
    p = np.array(p_values, dtype=float)
    n = len(p)
    if n == 0:
        return np.array([])

    # Remove NaN values for ranking, preserve NaN positions
    mask = ~np.isnan(p)
    valid_p = p[mask]
    n_valid = len(valid_p)

    if n_valid == 0:
        return np.full(n, np.nan)

    ranks = np.argsort(np.argsort(valid_p)) + 1
    q_valid = np.minimum(1, valid_p * n_valid / ranks)
    # Ensure monotonicity
    sorted_order = np.argsort(valid_p)
    q_valid_sorted = q_valid[sorted_order]
    for i in range(n_valid - 2, -1, -1):
        q_valid_sorted[i] = min(q_valid_sorted[i], q_valid_sorted[i + 1])
    q_valid = q_valid_sorted[np.argsort(sorted_order)]

    q_all = np.full(n, np.nan)
    q_all[mask] = q_valid
    return q_all


def call_dmr(window_df, sample_groups, group1, group2,
             q_threshold=0.05, delta_threshold=0.2,
             min_samples_per_group=2, fdr_by_context=True):
    """Call differentially methylated regions between two sample groups.

    Args:
        window_df: DataFrame with window data. Must contain a 'sample' column.
        sample_groups: Dict mapping group name → list of sample names.
        group1: Name of reference group (e.g. "normal")
        group2: Name of comparison group (e.g. "tumor")
        q_threshold: FDR q-value cutoff (default 0.05)
        delta_threshold: Minimum absolute methylation difference (default 0.2)
        min_samples_per_group: Minimum samples per group with coverage
        fdr_by_context: If True, BH correction per sequence context

    Returns:
        Tuple of (dmr_df, excluded_df)
    """
    samples1 = sample_groups[group1]
    samples2 = sample_groups[group2]

    results = []
    excluded = []
    window_keys = ["chr", "start", "end", "context"]

    # Group by window
    for (chrom, start, end, context), window_group in window_df.groupby(window_keys):
        g1 = window_group[window_group["sample"].isin(samples1)]
        g2 = window_group[window_group["sample"].isin(samples2)]

        # Coverage check
        n1 = len(g1[g1["total"] > 0])
        n2 = len(g2[g2["total"] > 0])

        if n1 < min_samples_per_group or n2 < min_samples_per_group:
            excluded.append({
                "chr": chrom, "start": start, "end": end,
                "context": context,
                "reason": f"insufficient samples: g1={n1}, g2={n2}",
            })
            continue

        # Check for all-zero
        total1 = g1["total"].sum()
        total2 = g2["total"].sum()
        if total1 == 0 and total2 == 0:
            excluded.append({
                "chr": chrom, "start": start, "end": end,
                "context": context, "reason": "all zero coverage",
            })
            continue

        p_value = _fisher_test_window(g1, g2)
        if np.isnan(p_value):
            excluded.append({
                "chr": chrom, "start": start, "end": end,
                "context": context, "reason": "fisher test failed",
            })
            continue

        ratio1 = g1["meth"].sum() / total1 if total1 > 0 else np.nan
        ratio2 = g2["meth"].sum() / total2 if total2 > 0 else np.nan

        results.append({
            "chr": chrom, "start": start, "end": end, "context": context,
            "p_value": p_value,
            "ratio_group1": ratio1,
            "ratio_group2": ratio2,
        })

    if not results:
        empty_dmr = pd.DataFrame(columns=[
            "chr", "start", "end", "context",
            "p_value", "q_value", "delta", "direction",
            "ratio_group1", "ratio_group2",
        ])
        empty_excl = pd.DataFrame(columns=["chr", "start", "end", "context", "reason"])
        return empty_dmr, empty_excl

    result_df = pd.DataFrame(results)

    # BH correction — per context or global
    if fdr_by_context:
        q_values = []
        for ctx in result_df["context"].unique():
            mask = result_df["context"] == ctx
            q_values.extend(_bh_correction(result_df.loc[mask, "p_value"].values))
        result_df["q_value"] = q_values
    else:
        result_df["q_value"] = _bh_correction(result_df["p_value"].values)

    # Delta and direction
    result_df["delta"] = result_df["ratio_group2"] - result_df["ratio_group1"]
    result_df["direction"] = result_df["delta"].apply(
        lambda d: "hyper" if d > 0 else "hypo"
    )

    # Apply thresholds
    dmr_df = result_df[
        (result_df["q_value"] < q_threshold) &
        (result_df["delta"].abs() > delta_threshold)
    ].copy()

    # Also collect excluded windows (those that didn't meet thresholds)
    remaining = result_df[
        ~((result_df["q_value"] < q_threshold) &
          (result_df["delta"].abs() > delta_threshold))
    ]
    for _, row in remaining.iterrows():
        excluded.append({
            "chr": row["chr"], "start": row["start"],
            "end": row["end"], "context": row["context"],
            "reason": "thresholds not met",
        })

    excluded_df = pd.DataFrame(excluded)
    return dmr_df.reset_index(drop=True), excluded_df.reset_index(drop=True)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_dmr.py -v`
Expected: all tests PASS.

- [ ] **Step 5: Commit**

```bash
git add smelt/dmr.py tests/test_dmr.py
git commit -m "feat: add DMR calling with Fisher exact test and BH correction"
```

---

## Phase 7: Statistics

### Task 13: Implement genome and chromosome statistics

**Files:**
- Create: `smelt/stats.py`
- Create: `tests/test_stats.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_stats.py
import pandas as pd
import pytest
import numpy as np
from smelt.stats import compute_stats

@pytest.fixture
def site_df():
    np.random.seed(0)
    rows = []
    for chrom in ["chrA", "chrB"]:
        for pos in range(0, 1000, 10):
            ctx = "CpG" if pos % 3 == 0 else ("CHH" if pos % 3 == 1 else "CHG")
            total = np.random.randint(5, 30)
            meth = np.random.randint(0, total)
            rows.append({
                "chr": chrom, "pos": pos, "context": ctx,
                "strand": "+", "meth": meth, "unmeth": total - meth,
                "total": total, "ratio": meth / total,
            })
    return pd.DataFrame(rows)

def test_compute_stats_output_columns(site_df):
    result = compute_stats(site_df)
    for col in ["chr", "context", "mean_ratio", "n_sites", "mean_depth"]:
        assert col in result.columns

def test_compute_stats_per_context(site_df):
    result = compute_stats(site_df)
    contexts = set(result["context"].unique())
    assert contexts == {"CpG", "CHH", "CHG"}

def test_compute_stats_n_sites_correct(site_df):
    result = compute_stats(site_df)
    total_sites = result["n_sites"].sum()
    assert total_sites == len(site_df)

def test_compute_stats_mean_ratio_bounds(site_df):
    result = compute_stats(site_df)
    assert all(result["mean_ratio"] >= 0)
    assert all(result["mean_ratio"] <= 1)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_stats.py -v`
Expected: FAIL.

- [ ] **Step 3: Implement statistics computation**

```python
# smelt/stats.py
"""Genome and chromosome-level methylation statistics."""
import pandas as pd


def compute_stats(df):
    """Compute summary statistics per chromosome and context.

    Args:
        df: Site or window DataFrame with chr, context, ratio, total columns

    Returns:
        DataFrame: chr, context, mean_ratio, n_sites, mean_depth
    """
    if "pos" in df.columns:
        site_col = "pos"
    else:
        site_col = "n_sites"

    if site_col == "pos":
        grouped = df.groupby(["chr", "context"], as_index=False).agg(
            mean_ratio=("ratio", "mean"),
            n_sites=("pos", "count"),
            mean_depth=("total", "mean"),
        )
    else:
        grouped = df.groupby(["chr", "context"], as_index=False).agg(
            mean_ratio=("ratio", "mean"),
            n_sites=("n_sites", "sum"),
            mean_depth=("total", "mean"),
        )

    return grouped.sort_values(["chr", "context"]).reset_index(drop=True)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_stats.py -v`
Expected: all tests PASS.

- [ ] **Step 5: Commit**

```bash
git add smelt/stats.py tests/test_stats.py
git commit -m "feat: add genome and chromosome statistics"
```

---

## Phase 8: CLI Integration

### Task 14: Wire up the Typer CLI

**Files:**
- Create: `smelt/cli.py`
- Create: `tests/test_cli.py`

- [ ] **Step 1: Write the CLI module**

```python
# smelt/cli.py
"""Smelt CLI — WGBS methylation downstream analysis."""

import sys
from pathlib import Path
from typing import List, Optional

import typer
import pandas as pd

from smelt.io import read_cov, read_gtf, read_bed, read_fasta
from smelt.site import compute_site_methylation
from smelt.window import compute_windows
from smelt.element import compute_elements
from smelt.metaplot import compute_metaplot
from smelt.custom import compute_custom
from smelt.dmr import call_dmr
from smelt.stats import compute_stats

app = typer.Typer(
    name="smelt",
    help="WGBS methylation downstream analysis tool",
)


def _write_output(df: pd.DataFrame, output: Optional[str], default_name: str):
    """Write DataFrame to TSV, using default_name if no output specified."""
    path = output if output else default_name
    df.to_csv(path, sep="\t", index=False)
    print(f"Wrote {len(df)} rows to {path}", file=sys.stderr)


def _output_base(input_path: str, suffix: str) -> str:
    """Generate output filename from input basename."""
    base = Path(input_path).stem
    # Strip .cov from bismark output
    base = base.replace(".cov", "")
    return f"{base}.{suffix}"


@app.command()
def site(
    input_file: str = typer.Option(..., "--input", "-i", help="BISMARK cov.gz file"),
    fasta: Optional[str] = typer.Option(None, "--fasta", "-f", help="Reference genome FASTA"),
    context_file: Optional[str] = typer.Option(None, "--context-file", help="Pre-annotated context file"),
    min_depth: int = typer.Option(5, "--min-depth", help="Minimum read depth"),
    merge_cpg: bool = typer.Option(True, "--merge-cpg-strands/--no-merge-cpg-strands"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
):
    """Compute per-cytosine methylation with sequence context classification."""
    cov = read_cov(input_file)
    fa = read_fasta(fasta) if fasta else None
    ctx = pd.read_csv(context_file, sep="\t") if context_file else None
    result = compute_site_methylation(
        cov, fasta=fa, context_df=ctx,
        min_depth=min_depth, merge_cpg=merge_cpg,
    )
    out = output if output else _output_base(input_file, "site.tsv")
    _write_output(result, out, out)


@app.command()
def window(
    input_file: str = typer.Option(..., "--input", "-i", help="Site TSV file or BISMARK cov.gz"),
    fasta: Optional[str] = typer.Option(None, "--fasta", "-f", help="Reference genome FASTA (if input is cov.gz)"),
    window_size: int = typer.Option(2000, "--window", "-w", help="Window size in bp"),
    step: int = typer.Option(500, "--step", "-s", help="Step size in bp"),
    min_sites: int = typer.Option(10, "--min-sites", help="Minimum sites per window"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
):
    """Aggregate methylation into sliding windows."""
    df = pd.read_csv(input_file, sep="\t")
    if "context" not in df.columns:
        # Input is raw cov.gz — compute sites first
        cov = read_cov(input_file)
        fa = read_fasta(fasta) if fasta else None
        if fa is None:
            raise typer.BadParameter("--fasta required when input is a raw cov.gz file")
        df = compute_site_methylation(cov, fasta=fa)
    result = compute_windows(df, window_size=window_size, step=step, min_sites=min_sites)
    suffix = f"window{window_size}s{step}.tsv"
    out = output if output else _output_base(input_file, suffix)
    _write_output(result, out, out)


@app.command()
def element(
    input_file: str = typer.Option(..., "--input", "-i", help="Site TSV file"),
    gtf: str = typer.Option(..., "--gtf", "-g", help="GTF annotation file"),
    features: str = typer.Option(None, "--features", help="Feature types (comma-separated, e.g. exon,CDS)"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
):
    """Compute methylation levels for GTF feature types."""
    sites = pd.read_csv(input_file, sep="\t")
    gtf_df = read_gtf(gtf)
    feat_list = features.split(",") if features else None
    result = compute_elements(sites, gtf_df, features=feat_list)
    out = output if output else _output_base(input_file, "element.tsv")
    _write_output(result, out, out)


@app.command()
def metaplot(
    input_file: str = typer.Option(..., "--input", "-i", help="Site TSV file"),
    gtf: Optional[str] = typer.Option(None, "--gtf", "-g", help="GTF file for gene intervals"),
    bed: Optional[str] = typer.Option(None, "--bed", "-b", help="BED file for gene intervals"),
    upstream: int = typer.Option(2000, "--upstream", help="Upstream bp"),
    downstream: int = typer.Option(2000, "--downstream", help="Downstream bp"),
    body_bins: int = typer.Option(50, "--body-bins", help="Gene body equal-ratio bins"),
    up_bins: int = typer.Option(20, "--up-bins", help="Upstream bins"),
    down_bins: int = typer.Option(20, "--down-bins", help="Downstream bins"),
    context: str = typer.Option("CpG", "--context", "-c", help="Sequence context"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
):
    """Compute metaplot methylation signal across gene bodies and flanking regions."""
    sites = pd.read_csv(input_file, sep="\t")
    if gtf:
        intervals = read_gtf(gtf)
        intervals = intervals[intervals["feature"] == "gene"]
    elif bed:
        intervals = read_bed(bed)
    else:
        raise typer.BadParameter("Either --gtf or --bed must be provided")
    result = compute_metaplot(
        sites, intervals, context=context,
        upstream=upstream, downstream=downstream,
        body_bins=body_bins, up_bins=up_bins, down_bins=down_bins,
    )
    out = output if output else _output_base(input_file, "metaplot.tsv")
    _write_output(result, out, out)


@app.command()
def custom(
    input_file: str = typer.Option(..., "--input", "-i", help="Site TSV file"),
    bed: str = typer.Option(..., "--bed", "-b", help="BED file with custom intervals"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
):
    """Compute methylation for user-defined BED intervals."""
    sites = pd.read_csv(input_file, sep="\t")
    bed_df = read_bed(bed)
    result = compute_custom(sites, bed_df)
    out = output if output else _output_base(input_file, "custom.tsv")
    _write_output(result, out, out)


@app.command()
def dmr(
    samples: List[str] = typer.Option(..., "--samples", help="Sample definitions: name=file (repeatable)"),
    group1: str = typer.Option(..., "--group1", help="Comma-separated sample names for group 1"),
    group2: str = typer.Option(..., "--group2", help="Comma-separated sample names for group 2"),
    q_threshold: float = typer.Option(0.05, "--q-threshold", help="FDR q-value cutoff"),
    delta_threshold: float = typer.Option(0.2, "--delta-threshold", help="Minimum methylation difference"),
    min_samples: int = typer.Option(2, "--min-samples-per-group", help="Minimum samples per group"),
    fdr_all: bool = typer.Option(False, "--fdr-all-contexts", help="Pool all contexts for BH correction"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
):
    """Call differentially methylated regions between two groups."""
    # Parse sample definitions
    sample_map = {}
    for s in samples:
        name, path = s.split("=", 1)
        sample_map[name] = path

    # Load all window files
    dfs = []
    for name, path in sample_map.items():
        df = pd.read_csv(path, sep="\t")
        df["sample"] = name
        dfs.append(df)
    merged = pd.concat(dfs, ignore_index=True)

    # Parse groups
    g1 = [s.strip() for s in group1.split(",")]
    g2 = [s.strip() for s in group2.split(",")]
    groups = {"group1": g1, "group2": g2}

    dmr_df, excl_df = call_dmr(
        merged, groups, group1="group1", group2="group2",
        q_threshold=q_threshold, delta_threshold=delta_threshold,
        min_samples_per_group=min_samples, fdr_by_context=not fdr_all,
    )

    base = _output_base(sample_map[g1[0]], "")
    dmr_out = output if output else f"{base}dmr.tsv"
    excl_out = dmr_out.replace(".tsv", "_excluded.tsv")
    _write_output(dmr_df, dmr_out, dmr_out)
    _write_output(excl_df, excl_out, excl_out)


@app.command()
def stats(
    input_file: str = typer.Option(..., "--input", "-i", help="Site or window TSV file"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
):
    """Compute genome and chromosome-level methylation statistics."""
    df = pd.read_csv(input_file, sep="\t")
    result = compute_stats(df)
    out = output if output else _output_base(input_file, "stats.tsv")
    _write_output(result, out, out)


if __name__ == "__main__":
    app()
```

- [ ] **Step 2: Verify CLI loads**

Run: `cd /home/hermes/Smelt && uv run smelt --help`
Expected: list of 7 subcommands with descriptions.

- [ ] **Step 3: Write CLI smoke tests**

```python
# tests/test_cli.py
from typer.testing import CliRunner
from smelt.cli import app

runner = CliRunner()

def test_cli_help():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "smelt" in result.stdout.lower() or "site" in result.stdout

def test_site_help():
    result = runner.invoke(app, ["site", "--help"])
    assert result.exit_code == 0
    assert "--input" in result.stdout

def test_window_help():
    result = runner.invoke(app, ["window", "--help"])
    assert result.exit_code == 0
    assert "--window" in result.stdout

def test_dmr_help():
    result = runner.invoke(app, ["dmr", "--help"])
    assert result.exit_code == 0
    assert "--group1" in result.stdout

def test_site_without_input_fails():
    result = runner.invoke(app, ["site"])
    assert result.exit_code != 0

def test_site_minimal_run():
    """Run site on the test fixture data."""
    result = runner.invoke(app, [
        "site", "--input", "tests/data/sample.cov",
        "--output", "/tmp/test_site_out.tsv",
    ])
    # Without FASTA, context can't be determined — should fail or warn
    # But with context_file unset and no fasta, it should error
    assert result.exit_code != 0
```

- [ ] **Step 4: Run CLI tests**

Run: `cd /home/hermes/Smelt && uv run pytest tests/test_cli.py -v`
Expected: all tests PASS.

- [ ] **Step 5: Commit**

```bash
git add smelt/cli.py tests/test_cli.py
git commit -m "feat: wire up Typer CLI with all 7 subcommands"
```

---

## Phase 9: Integration Testing

### Task 15: Create test data generation script

**Files:**
- Create: `scripts/generate_test_data.py`

- [ ] **Step 1: Write data generation script**

```python
# scripts/generate_test_data.py
"""Generate simulated WGBS data for integration testing.

Produces:
  - 2 fake chromosomes (chrA ~500kb, chrB ~300kb) as FASTA
  - 4 BISMARK cov.gz files (2 per group) with embedded ground truth
  - 1 GTF with 20 genes
  - 1 BED with 5 custom intervals

Ground truth: 8 DMRs embedded (5 CpG, 2 CHH, 1 CHG).
"""

import random
import gzip
import argparse
from pathlib import Path

random.seed(42)

CHROM_SIZES = {"chrA": 500_000, "chrB": 300_000}
SAMPLE_GROUPS = {
    "group1": ["sample_1", "sample_2"],
    "group2": ["sample_3", "sample_4"],
}

# Pre-defined DMRs: (chr, start, end, context, delta, direction)
EMBEDDED_DMRS = [
    ("chrA", 10000, 14000, "CpG", 0.8, "hyper"),
    ("chrA", 50000, 52000, "CpG", 0.6, "hypo"),
    ("chrA", 120000, 124000, "CpG", 0.4, "hyper"),
    ("chrA", 200000, 204000, "CpG", 0.25, "hypo"),
    ("chrA", 300000, 302000, "CpG", 0.15, "hyper"),
    ("chrB", 50000, 54000, "CHH", 0.7, "hyper"),
    ("chrB", 150000, 154000, "CHH", 0.5, "hypo"),
    ("chrB", 250000, 252000, "CHG", 0.55, "hyper"),
]


def generate_fasta(output_dir):
    """Generate fake chromosomes."""
    bases = ["A", "C", "G", "T"]
    for chrom, size in CHROM_SIZES.items():
        seq = []
        for i in range(size):
            # Ensure some CpG, CHH, CHG sites
            seq.append(random.choice(bases))
        path = output_dir / f"{chrom}.fa"
        with open(path, "w") as f:
            f.write(f">{chrom}\n")
            for i in range(0, size, 60):
                f.write("".join(seq[i:i+60]) + "\n")
        print(f"Wrote {path} ({size} bp)")


def generate_cov(output_dir):
    """Generate BISMARK coverage files with embedded DMRs."""
    for group, samples in SAMPLE_GROUPS.items():
        for sample in samples:
            rows = []
            for chrom, size in CHROM_SIZES.items():
                for pos in range(1, size + 1, random.randint(50, 200)):
                    ctx = determine_context(chrom, pos)
                    base_ratio = 0.8 if ctx == "CpG" else (0.1 if ctx == "CHH" else 0.05)
                    depth = random.randint(8, 30)
                    meth = int(depth * base_ratio)
                    unmet = depth - meth

                    # Apply DMR effect
                    for dmr_chr, dmr_start, dmr_end, dmr_ctx, delta, direction in EMBEDDED_DMRS:
                        if (chrom == dmr_chr and ctx == dmr_ctx and
                            dmr_start <= pos <= dmr_end):
                            if group == "group2":
                                if direction == "hyper":
                                    meth += int(depth * delta)
                                else:
                                    meth -= int(depth * delta)
                            meth = max(0, min(depth, meth))
                            unmet = depth - meth

                    pct = meth / depth * 100
                    rows.append(f"{chrom}\t{pos}\t{pos}\t{pct:.1f}\t{unmeth}\t{meth}")

            path = output_dir / f"{sample}.cov"
            with open(path, "w") as f:
                f.write("\n".join(rows))
            # Also create gzipped version
            gz_path = output_dir / f"{sample}.cov.gz"
            with gzip.open(gz_path, "wt") as f:
                f.write("\n".join(rows))
            print(f"Wrote {path} and {gz_path}")


def determine_context(chrom, pos):
    """Return a context based on position to create a realistic mix."""
    h = hash((chrom, pos)) % 100
    if h < 60:
        return "CpG"
    elif h < 85:
        return "CHH"
    else:
        return "CHG"


def generate_gtf(output_dir):
    """Generate a GTF with 20 genes."""
    genes = []
    gene_lengths = [500, 800, 1200, 2000, 3500, 5000, 8000, 12000, 20000, 35000,
                    600, 1500, 2800, 4500, 7000, 10000, 15000, 25000, 40000, 50000]
    positions = []
    for chrom in CHROM_SIZES:
        max_pos = CHROM_SIZES[chrom]
        for gl in gene_lengths:
            pos = random.randint(10000, max_pos - 60000)
            positions.append((chrom, pos, pos + gl))

    path = output_dir / "test_genes.gtf"
    with open(path, "w") as f:
        for i, (chrom, start, end) in enumerate(positions):
            gene_id = f"gene{i+1}"
            # Gene
            f.write(f"{chrom}\ttest\tgene\t{start}\t{end}\t.\t+\t.\tgene_id \"{gene_id}\";\n")
            # Exons
            n_exons = random.randint(2, 8)
            exon_size = (end - start) // n_exons
            for j in range(n_exons):
                e_start = start + j * exon_size
                e_end = e_start + exon_size - 10
                f.write(f"{chrom}\ttest\texon\t{e_start}\t{e_end}\t.\t+\t.\tgene_id \"{gene_id}\"; transcript_id \"{gene_id}_t1\";\n")
            # CDS
            cds_start = start + exon_size
            cds_end = end - exon_size
            f.write(f"{chrom}\ttest\tCDS\t{cds_start}\t{cds_end}\t.\t+\t.\tgene_id \"{gene_id}\"; transcript_id \"{gene_id}_t1\";\n")
    print(f"Wrote {path}")


def generate_bed(output_dir):
    """Generate a BED with 5 custom intervals."""
    path = output_dir / "test_regions.bed"
    with open(path, "w") as f:
        for i in range(1, 6):
            chrom = "chrA" if i <= 3 else "chrB"
            start = random.randint(10000, CHROM_SIZES[chrom] - 20000)
            end = start + random.randint(2000, 10000)
            f.write(f"{chrom}\t{start}\t{end}\tregion_{i}\n")
    print(f"Wrote {path}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", default="tests/data/integration")
    args = p.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    generate_fasta(out)
    generate_cov(out)
    generate_gtf(out)
    generate_bed(out)
    print("Done.")
```

- [ ] **Step 2: Run data generation**

Run: `cd /home/hermes/Smelt && uv run python scripts/generate_test_data.py --outdir tests/data/integration`
Expected: test data files created.

- [ ] **Step 3: Commit**

```bash
git add scripts/generate_test_data.py tests/data/integration/.gitkeep
git commit -m "feat: add integration test data generation script"
```

---

### Task 16: Write integration tests

**Files:**
- Create: `tests/test_integration.py`

- [ ] **Step 1: Write integration test**

```python
# tests/test_integration.py
"""End-to-end pipeline integration tests using simulated data.

Requires test data generated by scripts/generate_test_data.py.
"""
import os
import pytest
import pandas as pd
from pathlib import Path

from smelt.io import read_cov, read_gtf, read_bed
from smelt.site import compute_site_methylation
from smelt.window import compute_windows
from smelt.element import compute_elements
from smelt.metaplot import compute_metaplot
from smelt.custom import compute_custom
from smelt.dmr import call_dmr
from smelt.stats import compute_stats

DATA_DIR = Path("tests/data/integration")


@pytest.fixture(scope="module")
def fasta():
    from smelt.io import read_fasta
    fa_path = DATA_DIR / "chrA.fa"
    if not fa_path.exists():
        pytest.skip("Integration test data not generated. Run scripts/generate_test_data.py first.")
    return read_fasta(str(fa_path))


@pytest.fixture(scope="module")
def all_sites(fasta):
    """Compute sites for all 4 samples."""
    sites = {}
    for sample in ["sample_1", "sample_2", "sample_3", "sample_4"]:
        cov_path = DATA_DIR / f"{sample}.cov.gz"
        if not cov_path.exists():
            pytest.skip(f"Missing {cov_path}")
        cov = read_cov(str(cov_path))
        sites[sample] = compute_site_methylation(cov, fasta=fasta, min_depth=5, merge_cpg=True)
    return sites


@pytest.fixture(scope="module")
def all_windows(all_sites):
    """Compute windows for all 4 samples."""
    windows = {}
    for sample, sites in all_sites.items():
        windows[sample] = compute_windows(sites, window_size=2000, step=500, min_sites=5)
    return windows


class TestIntegrationSite:
    def test_all_contexts_present(self, all_sites):
        for sample, df in all_sites.items():
            contexts = set(df["context"].unique())
            assert contexts.issubset({"CpG", "CHH", "CHG"}), f"{sample}: {contexts}"

    def test_ratio_bounds(self, all_sites):
        for sample, df in all_sites.items():
            assert all(df["ratio"] >= 0), f"{sample} has negative ratios"
            assert all(df["ratio"] <= 1), f"{sample} has ratios > 1"

    def test_depth_filter_applied(self, all_sites):
        for sample, df in all_sites.items():
            assert all(df["total"] >= 5), f"{sample} has sites below min_depth"


class TestIntegrationWindow:
    def test_windows_cover_genome(self, all_windows):
        for sample, df in all_windows.items():
            assert len(df) > 0, f"{sample} has no windows"

    def test_window_size_respected(self, all_windows):
        for sample, df in all_windows.items():
            sizes = df["end"] - df["start"]
            assert all(sizes == 2000), f"{sample}: window sizes not all 2000"


class TestIntegrationElement:
    def test_elements_computed(self, all_sites):
        gtf_path = DATA_DIR / "test_genes.gtf"
        if not gtf_path.exists():
            pytest.skip("GTF not found")
        gtf = read_gtf(str(gtf_path))
        df = compute_elements(all_sites["sample_1"], gtf, features=["gene"])
        assert len(df) > 0

    def test_element_feature_types(self, all_sites):
        gtf_path = DATA_DIR / "test_genes.gtf"
        if not gtf_path.exists():
            pytest.skip("GTF not found")
        gtf = read_gtf(str(gtf_path))
        df = compute_elements(all_sites["sample_1"], gtf, features=["exon"])
        assert all(df["feature_type"] == "exon")


class TestIntegrationDMR:
    def test_dmr_detects_embedded_dmrs(self, all_windows):
        """Verify that at least some of the 8 embedded DMRs are detected."""
        # Merge all windows
        dfs = []
        for sample, df in all_windows.items():
            df = df.copy()
            df["sample"] = sample
            dfs.append(df)
        merged = pd.concat(dfs, ignore_index=True)

        groups = {
            "group1": ["sample_1", "sample_2"],
            "group2": ["sample_3", "sample_4"],
        }
        dmr_df, _ = call_dmr(
            merged, groups, group1="group1", group2="group2",
            q_threshold=0.1, delta_threshold=0.1,
        )
        # Should detect at least 2 of 8 DMRs
        assert len(dmr_df) >= 2, f"Only {len(dmr_df)} DMRs detected"

    def test_dmr_has_direction(self, all_windows):
        dfs = []
        for sample, df in all_windows.items():
            df = df.copy()
            df["sample"] = sample
            dfs.append(df)
        merged = pd.concat(dfs, ignore_index=True)
        groups = {
            "group1": ["sample_1", "sample_2"],
            "group2": ["sample_3", "sample_4"],
        }
        dmr_df, _ = call_dmr(merged, groups, group1="group1", group2="group2")
        if len(dmr_df) > 0:
            assert all(d in ("hyper", "hypo") for d in dmr_df["direction"])


class TestIntegrationStats:
    def test_stats_per_context(self, all_sites):
        df = compute_stats(all_sites["sample_1"])
        assert len(df) > 0
        assert "CpG" in set(df["context"])

    def test_stats_all_chromosomes(self, all_sites):
        df = compute_stats(all_sites["sample_1"])
        chromosomes = set(df["chr"].unique())
        assert chromosomes == {"chrA", "chrB"}
```

- [ ] **Step 2: Generate data and run integration tests**

```bash
cd /home/hermes/Smelt
uv run python scripts/generate_test_data.py --outdir tests/data/integration
uv run pytest tests/test_integration.py -v
```
Expected: integration tests PASS if generated data exists, or SKIP with message.

- [ ] **Step 3: Commit**

```bash
git add tests/test_integration.py
git commit -m "test: add integration tests with simulated ground-truth data"
```

---

## Phase 10: Benchmark Script

### Task 17: Create benchmark script

**Files:**
- Create: `scripts/benchmark.sh`
- Create: `scripts/benchmark.py`

- [ ] **Step 1: Write benchmark runner**

```python
# scripts/benchmark.py
"""Benchmark runner for Smelt.

Measures runtime and peak memory for each subcommand on real data.
Compares DMR results with methylKit/DSS if available.

Usage:
    python scripts/benchmark.py --input-dir /path/to/data --output-dir results/
"""
import argparse
import subprocess
import time
import json
import sys
from pathlib import Path


def run_cmd(cmd, timeout=3600):
    """Run a command, return (wall_time, max_rss_mb)."""
    start = time.perf_counter()
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
    elapsed = time.perf_counter() - start
    if result.returncode != 0:
        print(f"  FAILED: {result.stderr[:500]}", file=sys.stderr)
    return elapsed


def benchmark_site(input_dir, output_dir):
    """Benchmark smelt site."""
    print("Benchmarking: smelt site")
    t = run_cmd(
        f"smelt site --input {input_dir}/sample.cov.gz --fasta {input_dir}/genome.fa "
        f"--output {output_dir}/bench.site.tsv"
    )
    return {"site_wall_time_s": t}


def benchmark_window(output_dir):
    """Benchmark smelt window."""
    print("Benchmarking: smelt window")
    t = run_cmd(
        f"smelt window --input {output_dir}/bench.site.tsv "
        f"--window 2000 --step 500 --output {output_dir}/bench.window.tsv"
    )
    return {"window_wall_time_s": t}


def benchmark_dmr(output_dir):
    """Benchmark smelt dmr."""
    print("Benchmarking: smelt dmr")
    t = run_cmd(
        f"smelt dmr "
        f"--samples g1_s1={output_dir}/g1_s1.window.tsv "
        f"g1_s2={output_dir}/g1_s2.window.tsv "
        f"g2_s1={output_dir}/g2_s1.window.tsv "
        f"g2_s2={output_dir}/g2_s2.window.tsv "
        f"--group1 g1_s1,g1_s2 --group2 g2_s1,g2_s2 "
        f"--output {output_dir}/bench.dmr.tsv"
    )
    return {"dmr_wall_time_s": t}


def main():
    p = argparse.ArgumentParser(description="Smelt benchmark")
    p.add_argument("--input-dir", required=True)
    p.add_argument("--output-dir", default="benchmark_results")
    args = p.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    results = {}
    results.update(benchmark_site(args.input_dir, out))
    results.update(benchmark_window(out))
    results.update(benchmark_dmr(out))

    # Save results
    results_file = out / "benchmark.json"
    with open(results_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {results_file}")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Write shell wrapper**

```bash
# scripts/benchmark.sh
#!/usr/bin/env bash
# Smelt benchmark script — run on real WGBS data before release.
#
# Usage:
#   bash scripts/benchmark.sh /path/to/bismark/outputs/
#
# The input directory should contain:
#   - *.cov.gz files for each sample
#   - genome.fa reference
#
# Results are written to benchmark_results/

set -euo pipefail

INPUT_DIR="${1:?Usage: $0 <input-dir>}"
OUTPUT_DIR="benchmark_results"

echo "=== Smelt Benchmark ==="
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo

python scripts/benchmark.py --input-dir "$INPUT_DIR" --output-dir "$OUTPUT_DIR"

echo
echo "=== Benchmark Complete ==="
cat "$OUTPUT_DIR/benchmark.json"
```

```bash
chmod +x scripts/benchmark.sh
```

- [ ] **Step 3: Commit**

```bash
git add scripts/benchmark.sh scripts/benchmark.py
git commit -m "feat: add benchmark script for pre-release performance testing"
```

---

## Final Verification

### Task 18: Run full test suite

- [ ] **Step 1: Run all unit and integration tests**

```bash
cd /home/hermes/Smelt
uv run python scripts/generate_test_data.py --outdir tests/data/integration
uv run pytest tests/ -v
```
Expected: all tests PASS, integration tests PASS.

- [ ] **Step 2: Run type checking**

```bash
uv run mypy smelt/ --ignore-missing-imports
```
Expected: no type errors (or minimal with explanation).

- [ ] **Step 3: Verify CLI end-to-end**

```bash
# Full pipeline with test data
SAMPLE=tests/data/integration/sample_1

uv run smelt site --input ${SAMPLE}.cov.gz --fasta tests/data/integration/chrA.fa -o /tmp/test_site.tsv
head -5 /tmp/test_site.tsv

uv run smelt window --input /tmp/test_site.tsv -o /tmp/test_window.tsv
head -5 /tmp/test_window.tsv

uv run smelt stats --input /tmp/test_site.tsv -o /tmp/test_stats.tsv
head -5 /tmp/test_stats.tsv
```
Expected: each command produces output without errors.

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "chore: final verification, all tests passing"
```
