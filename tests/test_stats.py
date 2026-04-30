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
    # Per-chromosome rows only (exclude genome row)
    per_chr = result[result["chr"] != "genome"]
    total_sites = per_chr["n_sites"].sum()
    assert total_sites == len(site_df)


def test_compute_stats_mean_ratio_bounds(site_df):
    result = compute_stats(site_df)
    assert all(result["mean_ratio"] >= 0)
    assert all(result["mean_ratio"] <= 1)
