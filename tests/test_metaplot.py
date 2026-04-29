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
    assert len(result) == 20


def test_compute_metaplot_ratio_between_0_and_1(gene_intervals, site_df):
    result = compute_metaplot(site_df, gene_intervals)
    assert all(result["ratio"] >= 0)
    assert all(result["ratio"] <= 1)
