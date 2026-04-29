import pandas as pd
import pytest
from smelt.window import compute_windows, _make_windows


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


def test_make_windows_produces_correct_intervals():
    windows = _make_windows("chrA", 0, 3000, window_size=2000, step=1000)
    assert len(windows) == 3
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
    # CpG: ~20 sites per window, CHH: ~10 sites.
    # With min_sites=50, all windows filtered out
    assert len(result) == 0
