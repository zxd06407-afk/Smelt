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
    assert len(result) == 2


def test_compute_custom_aggregates_correctly(site_df, bed_df):
    result = compute_custom(site_df, bed_df)
    r1 = result[result["interval_id"] == "region_1"].iloc[0]
    # region_1: pos 99-199 (0-based) -> 101 sites, each with meth=8, unmet=2
    assert r1["n_sites"] == 101
    assert r1["ratio"] == pytest.approx(0.8)
