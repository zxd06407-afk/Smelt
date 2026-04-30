import pytest
from smelt.utils import to_zero_based, to_one_based, classify_context, merge_cpg_strands
import pandas as pd


def test_to_zero_based_single_position():
    """1-based pos 100 -> 0-based pos 99."""
    assert to_zero_based(100) == 99


def test_to_zero_based_gtf_interval():
    """GTF 1-based closed [100, 200] -> 0-based half-open [99, 200)."""
    assert to_zero_based(100, is_interval=True, end=200) == (99, 200)


def test_to_zero_based_bed_already_zero():
    """BED 0-based [99, 200) -- should be a no-op."""
    assert to_zero_based(99, source="bed") == 99


def test_to_one_based_single_position():
    """0-based pos 99 -> 1-based pos 100."""
    assert to_one_based(99) == 100


def test_roundtrip():
    """to_zero_based then to_one_based returns original."""
    for pos in [1, 100, 1000]:
        assert to_one_based(to_zero_based(pos)) == pos


# --- context classification ---

def test_classify_cpg():
    assert classify_context("CG") == "CpG"


def test_classify_chg():
    assert classify_context("CAG") == "CHG"


def test_classify_chh():
    assert classify_context("CAT") == "CHH"


def test_classify_cpg_case_insensitive():
    assert classify_context("cg") == "CpG"


def test_classify_insufficient_bases_returns_none():
    assert classify_context("C") is None
    assert classify_context("") is None


# --- CpG strand merging ---

def test_merge_cpg_strands_combines_both_strands():
    df = pd.DataFrame({
        "chr": ["chrA", "chrA", "chrA"],
        "pos": [99, 100, 150],
        "context": ["CpG", "CHH", "CpG"],
        "strand": ["+", "-", "-"],
        "meth": [10, 3, 5],
        "unmeth": [2, 1, 5],
        "total": [12, 4, 10],
        "ratio": [0.833, 0.75, 0.5],
    })
    result = merge_cpg_strands(df)
    assert len(result) == 3  # one CHH + two CpG at different positions (not merged)
    cpg = result[result["context"] == "CpG"]
    assert len(cpg) == 2


def test_merge_cpg_dyad_adjacent_positions():
    """CpG at pos 99/+ and pos 100/- merge into dyad at pos 99."""
    df = pd.DataFrame({
        "chr": ["chrA", "chrA", "chrA"],
        "pos": [99, 100, 200],
        "context": ["CpG", "CpG", "CpG"],
        "strand": ["+", "-", "+"],
        "meth": [10, 5, 8],
        "unmeth": [2, 1, 3],
        "total": [12, 6, 11],
        "ratio": [0.833, 0.833, 0.727],
    })
    result = merge_cpg_strands(df)
    dyad = result[result["pos"] == 99]
    assert len(dyad) == 1
    assert dyad.iloc[0]["meth"] == 15
    assert dyad.iloc[0]["unmeth"] == 3
    lone = result[result["pos"] == 200]
    assert lone.iloc[0]["meth"] == 8


def test_merge_cpg_dyad_non_adjacent_not_merged():
    """CpG at pos 99/+ and pos 101/- should NOT merge."""
    df = pd.DataFrame({
        "chr": ["chrA", "chrA"],
        "pos": [99, 101],
        "context": ["CpG", "CpG"],
        "strand": ["+", "-"],
        "meth": [10, 5],
        "unmeth": [2, 1],
        "total": [12, 6],
        "ratio": [0.833, 0.833],
    })
    result = merge_cpg_strands(df)
    assert len(result) == 2


def test_merge_cpg_no_strand_column():
    df = pd.DataFrame({
        "chr": ["chrA"], "pos": [99], "context": ["CpG"],
        "meth": [10], "unmeth": [2], "total": [12], "ratio": [0.833],
    })
    result = merge_cpg_strands(df)
    assert len(result) == 1
