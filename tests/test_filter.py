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
    from smelt.filter import filter_by_sites
    df = pd.DataFrame({
        "chr": ["chrA", "chrA", "chrB"],
        "start": [0, 1000, 0],
        "end": [2000, 3000, 2000],
        "context": ["CpG", "CpG", "CpG"],
        "n_sites": [5, 15, 10],
    })
    result = filter_by_sites(df, min_sites=10)
    assert len(result) == 2
