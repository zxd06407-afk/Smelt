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


def test_compute_stats_genome_weighted_mean():
    """Genome mean_ratio is weighted by n_sites across chromosomes."""
    rows = []
    for pos in range(10):
        rows.append({"chr": "chrA", "pos": pos, "context": "CpG",
                       "strand": "+", "meth": 2, "unmeth": 8,
                       "total": 10, "ratio": 0.2})
    for pos in range(100):
        rows.append({"chr": "chrB", "pos": pos, "context": "CpG",
                       "strand": "+", "meth": 8, "unmeth": 2,
                       "total": 10, "ratio": 0.8})
    df = pd.DataFrame(rows)
    result = compute_stats(df)
    genome = result[(result["chr"] == "genome") & (result["context"] == "CpG")]
    assert len(genome) == 1
    assert genome.iloc[0]["n_sites"] == 110
    expected_weighted = (10 * 0.2 + 100 * 0.8) / 110
    assert genome.iloc[0]["mean_ratio"] == pytest.approx(expected_weighted)


def test_compute_stats_sample_name(site_df):
    result = compute_stats(site_df, sample_name="test_sample")
    assert "sample" in result.columns
    assert all(result["sample"] == "test_sample")
