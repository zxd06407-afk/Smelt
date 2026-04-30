import pandas as pd
import pytest
from smelt.element import compute_elements


@pytest.fixture
def site_df():
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


def test_compute_elements_exact_counts():
    """Known sites inside a feature → exact count verification."""
    site_df = pd.DataFrame({
        "chr": ["chrA", "chrA", "chrA"],
        "pos": [99, 100, 101],
        "context": ["CpG", "CpG", "CpG"],
        "strand": ["+", "+", "+"],
        "meth": [8, 8, 8],
        "unmeth": [2, 2, 2],
        "total": [10, 10, 10],
        "ratio": [0.8, 0.8, 0.8],
    })
    gtf_df = pd.DataFrame({
        "chr": ["chrA"],
        "source": ["test"], "feature": ["gene"],
        "start": [99], "end": [200],
        "score": ["."], "strand": ["+"], "frame": ["."],
        "gene_id": ["gene1"], "transcript_id": [None],
    })
    result = compute_elements(site_df, gtf_df, features=["gene"])
    assert len(result) == 1
    assert result.iloc[0]["n_sites"] == 3
    assert result.iloc[0]["meth"] == 24
    assert result.iloc[0]["unmeth"] == 6
    assert result.iloc[0]["ratio"] == pytest.approx(0.8)


def test_compute_elements_empty_result():
    """Feature types not in data → empty result, not crash."""
    site_df = pd.DataFrame({
        "chr": ["chrA"], "pos": [99], "context": ["CpG"],
        "strand": ["+"], "meth": [8], "unmeth": [2],
        "total": [10], "ratio": [0.8],
    })
    gtf_df = pd.DataFrame({
        "chr": ["chrA"], "source": ["test"], "feature": ["exon"],
        "start": [500], "end": [600],
        "score": ["."], "strand": ["+"], "frame": ["."],
        "gene_id": ["gene1"], "transcript_id": [None],
    })
    result = compute_elements(site_df, gtf_df, features=["exon"])
    assert len(result) == 0
