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


def test_compute_site_min_depth_filter(cov_df, fasta):
    result = compute_site_methylation(cov_df, fasta=fasta, min_depth=5)
    # pos 50 has total=4, should be filtered
    assert not any(result["total"] < 5)


def test_compute_site_output_columns(cov_df, fasta):
    result = compute_site_methylation(cov_df, fasta=fasta)
    expected = {"chr", "pos", "context", "strand", "meth", "unmeth", "total", "ratio"}
    assert expected.issubset(set(result.columns))


def test_compute_site_requires_fasta_or_context():
    df = pd.DataFrame({
        "chr": ["chrA"], "pos": [0], "meth": [5], "unmeth": [3],
        "total": [8], "ratio": [0.625],
    })
    with pytest.raises(ValueError, match="Either fasta or context_df"):
        compute_site_methylation(df)
