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
        "unmeth": [2, 5, 0, 3, 6],
        "total": [12, 10, 8, 6, 10],
        "ratio": [0.833, 0.5, 1.0, 0.5, 0.4],
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
