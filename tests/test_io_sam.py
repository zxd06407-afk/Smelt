import pytest
from smelt.io import read_bismark_sam


def test_read_bismark_sam_returns_dataframe():
    df = read_bismark_sam("tests/data/test.sam")
    assert len(df) > 0


def test_read_bismark_sam_columns():
    df = read_bismark_sam("tests/data/test.sam")
    expected = {"chr", "pos", "meth", "unmeth", "total", "ratio"}
    assert expected.issubset(set(df.columns))


def test_read_bismark_sam_counts():
    """read1 and read2 align at chrA:100 (1-based, 99 0-based).
    XM has 'Z' (meth CpG) at i%4==0 positions.
    pos=99 gets 2 meth from read1+read2."""
    df = read_bismark_sam("tests/data/test.sam")
    row = df[df["pos"] == 99]
    assert len(row) == 1
    assert row.iloc[0]["meth"] == 2
    assert row.iloc[0]["unmeth"] == 0


def test_read_bismark_sam_unmeth_counts():
    """XM has 'z' (unmeth CpG) at i%4==2 positions.
    pos=101 (0-based) gets 2 unmeth from read1+read2."""
    df = read_bismark_sam("tests/data/test.sam")
    row = df[df["pos"] == 101]
    assert len(row) == 1
    assert row.iloc[0]["unmeth"] == 2


def test_read_bismark_sam_mixed_strands():
    """read3 has 'H' at i%4==0 (meth CHH) and 'h' at i%4==2 (unmeth CHH).
    read4 has 'x' at i%4==2 (unmeth CHG).
    pos=149: read3 XM[0]='H' -> 1 meth"""
    df = read_bismark_sam("tests/data/test.sam")
    row = df[df["pos"] == 149]
    assert len(row) == 1
    assert row.iloc[0]["meth"] == 1
    assert row.iloc[0]["unmeth"] == 0


def test_read_bismark_sam_bam_input():
    """Verify BAM input also works."""
    df = read_bismark_sam("tests/data/test.bam")
    assert len(df) == 36


def test_read_bismark_sam_position_range():
    """Positions with data from read1+read2: 99-133 step 2.
    Positions from read3+read4: 149-183 step 2."""
    df = read_bismark_sam("tests/data/test.sam")
    positions = set(df["pos"])
    assert 99 in positions
    assert 133 in positions
    assert 149 in positions
    assert 183 in positions
