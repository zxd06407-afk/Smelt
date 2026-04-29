import pytest
import pandas as pd
from smelt.io import read_cov, read_gtf, read_bed, read_fasta


def test_read_cov_columns():
    df = read_cov("tests/data/sample.cov")
    expected_cols = {"chr", "pos", "meth", "unmeth", "total", "ratio"}
    assert expected_cols.issubset(set(df.columns))


def test_read_cov_row_count():
    df = read_cov("tests/data/sample.cov")
    assert len(df) == 6


def test_read_cov_values():
    df = read_cov("tests/data/sample.cov")
    first = df.iloc[0]
    assert first["chr"] == "chrA"
    assert first["pos"] == 99        # 0-based
    assert first["meth"] == 6
    assert first["unmeth"] == 2
    assert first["total"] == 8
    assert first["ratio"] == pytest.approx(0.75)


def test_read_cov_all_zero_unmeth():
    """Site with only methylated reads."""
    df = read_cov("tests/data/sample.cov")
    row = df[df["pos"] == 101]    # 0-based, original 102
    assert row.iloc[0]["unmeth"] == 0
    assert row.iloc[0]["ratio"] == pytest.approx(1.0)


# --- GTF ---

def test_read_gtf_columns():
    df = read_gtf("tests/data/sample.gtf")
    expected = {"chr", "start", "end", "source", "feature", "strand", "gene_id"}
    assert expected.issubset(set(df.columns))


def test_read_gtf_feature_types():
    df = read_gtf("tests/data/sample.gtf")
    assert set(df["feature"].unique()) == {"gene", "exon", "CDS"}


def test_read_gtf_is_zero_based():
    """GTF [100,500] 1-based -> [99,500) 0-based."""
    df = read_gtf("tests/data/sample.gtf")
    gene = df[df["feature"] == "gene"].iloc[0]
    assert gene["start"] == 99
    assert gene["end"] == 500


def test_read_gtf_parses_attributes():
    df = read_gtf("tests/data/sample.gtf")
    assert df.iloc[0]["gene_id"] == "gene1"


# --- BED ---

def test_read_bed_columns():
    df = read_bed("tests/data/sample.bed")
    assert list(df.columns[:4]) == ["chr", "start", "end", "name"]


def test_read_bed_preserves_zero_based():
    """BED is 0-based -- positions unchanged."""
    df = read_bed("tests/data/sample.bed")
    assert df.iloc[0]["start"] == 99
    assert df.iloc[0]["end"] == 200


# --- FASTA ---

def test_read_fasta_returns_faidx():
    fa = read_fasta("tests/data/sample.fa")
    assert "chrA" in fa
    assert len(fa["chrA"]) == 56
