import pandas as pd
import pytest
from smelt.window import compute_windows


def make_site_df():
    """Sites on chrA spanning 0-3000 with 3 contexts."""
    rows = []
    for pos in range(0, 3001, 100):
        rows.append({"chr": "chrA", "pos": pos, "context": "CpG",
                      "strand": "+", "meth": 8, "unmeth": 2, "total": 10, "ratio": 0.8})
        if pos % 200 == 0:
            rows.append({"chr": "chrA", "pos": pos + 50, "context": "CHH",
                          "strand": "+", "meth": 3, "unmeth": 7, "total": 10, "ratio": 0.3})
    return pd.DataFrame(rows)


def test_window_intervals_correct():
    """Windows cover the genome with correct step intervals."""
    df = pd.DataFrame({
        "chr": ["chrA", "chrA", "chrA", "chrA"],
        "pos": [500, 1500, 2500, 3500],
        "context": ["CpG"] * 4,
        "strand": ["+"] * 4,
        "meth": [5] * 4, "unmeth": [5] * 4,
        "total": [10] * 4, "ratio": [0.5] * 4,
    })
    result = compute_windows(df, window_size=2000, step=1000, min_sites=1)
    starts = sorted(result["start"].unique())
    # Windows should cover all sites: 0, 1000, 2000, 3000
    assert 0 in starts
    assert 1000 in starts
    assert 2000 in starts


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
    assert len(result) == 0


def test_compute_windows_boundary_precision():
    """Half-open [start, end): pos=1999 is inside, pos=2000 is outside."""
    df = pd.DataFrame({
        "chr": ["chrA", "chrA", "chrA"],
        "pos": [1999, 2000, 2001],
        "context": ["CpG", "CpG", "CpG"],
        "strand": ["+", "+", "+"],
        "meth": [5, 5, 5],
        "unmeth": [5, 5, 5],
        "total": [10, 10, 10],
        "ratio": [0.5, 0.5, 0.5],
    })
    result = compute_windows(df, window_size=2000, step=500, min_sites=1)
    w0 = result[result["start"] == 0]
    assert len(w0) == 1
    assert w0.iloc[0]["n_sites"] == 1  # only pos=1999
    # pos=2000 and pos=2001 belong to window [2000, 4000)
    w2000 = result[(result["start"] == 2000) & (result["context"] == "CpG")]
    assert len(w2000) == 1
    assert w2000.iloc[0]["n_sites"] == 2  # pos=2000,2001


def test_compute_windows_step_overlap():
    """Site at pos=1500 appears in exactly 4 windows with W=2000, step=500."""
    df = pd.DataFrame({
        "chr": ["chrA"],
        "pos": [1500],
        "context": ["CpG"],
        "strand": ["+"],
        "meth": [5],
        "unmeth": [5],
        "total": [10],
        "ratio": [0.5],
    })
    result = compute_windows(df, window_size=2000, step=500, min_sites=1)
    # Site at 1500 is in windows: [0,2000), [500,2500), [1000,3000), [1500,3500)
    # NOT in [2000,4000) since 1500 < 2000
    assert len(result) == 4
