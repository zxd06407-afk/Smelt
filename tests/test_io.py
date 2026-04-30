import pytest
import pandas as pd
from smelt.io import read_gtf, read_bed, read_fasta, read_cx_report


# --- CX_report.txt ---

def test_read_cx_report_columns():
    df = read_cx_report("tests/data/test.CX_report.txt")
    expected = {"chr", "pos", "strand", "meth", "unmeth", "context", "total", "ratio"}
    assert expected.issubset(set(df.columns))


def test_read_cx_report_strand():
    df = read_cx_report("tests/data/test.CX_report.txt")
    assert set(df["strand"].unique()) == {"+", "-"}


def test_read_cx_report_context():
    df = read_cx_report("tests/data/test.CX_report.txt")
    assert set(df["context"].unique()) == {"CpG", "CHH", "CHG"}


def test_read_cx_report_zero_based():
    """CX_report 1-based pos 100 -> 0-based pos 99."""
    df = read_cx_report("tests/data/test.CX_report.txt")
    assert df.iloc[0]["pos"] == 99


def test_read_cx_report_counts():
    df = read_cx_report("tests/data/test.CX_report.txt")
    first = df.iloc[0]
    assert first["meth"] == 6
    assert first["unmeth"] == 2
    assert first["total"] == 8


# --- GTF ---

def test_read_gtf_columns():
    df = read_gtf("tests/data/sample.gtf")
    expected = {"chr", "start", "end", "source", "feature", "strand", "gene_id"}
    assert expected.issubset(set(df.columns))


def test_read_gtf_feature_types():
    df = read_gtf("tests/data/sample.gtf")
    assert set(df["feature"].unique()) == {"gene", "exon", "CDS"}


def test_read_gtf_is_zero_based():
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
    df = read_bed("tests/data/sample.bed")
    assert df.iloc[0]["start"] == 99
    assert df.iloc[0]["end"] == 200


# --- FASTA ---

def test_read_fasta_returns_faidx():
    fa = read_fasta("tests/data/sample.fa")
    assert "chrA" in fa
    assert len(fa["chrA"]) == 56
