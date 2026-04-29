"""File I/O for BISMARK, GTF, BED, and FASTA formats."""
import re

import pandas as pd
import pyfaidx

from smelt.utils import to_zero_based


def read_cov(path):
    """Read a BISMARK coverage file (.cov or .cov.gz).

    BISMARK cov format (1-based):
        chr start end meth_pct unmeth meth

    Returns DataFrame with 0-based positions:
        chr, pos, meth, unmeth, total, ratio
    """
    cols = ["chr", "start", "end", "meth_pct", "unmeth", "meth"]
    df = pd.read_csv(
        path, sep="\t", names=cols, header=None,
        dtype={"chr": str},
    )
    df["pos"] = df["start"].apply(to_zero_based)
    df["total"] = df["meth"] + df["unmeth"]
    df["ratio"] = df["meth"] / df["total"]
    return df[["chr", "pos", "meth", "unmeth", "total", "ratio"]]


def read_gtf(path):
    """Read a GTF file, returning a DataFrame with 0-based half-open coordinates.

    Parses column 9 (attributes) for gene_id and transcript_id.
    """
    cols = ["chr", "source", "feature", "start", "end",
            "score", "strand", "frame", "attributes"]
    df = pd.read_csv(
        path, sep="\t", names=cols, header=None, comment="#",
        dtype={"chr": str, "start": int, "end": int},
    )
    df["start"] = df["start"] - 1

    df["gene_id"] = df["attributes"].apply(_parse_attr, key="gene_id")
    df["transcript_id"] = df["attributes"].apply(_parse_attr, key="transcript_id")
    return df


def _parse_attr(attr_string, key):
    """Extract a key from GTF attribute string.

    Supports both 'key "value";' and 'key=value;' formats.
    """
    pattern = rf'{key}\s+"?([^";]+)"?'
    m = re.search(pattern, attr_string)
    return m.group(1) if m else None


def read_bed(path):
    """Read a BED file. BED is 0-based half-open, no conversion needed."""
    cols = ["chr", "start", "end", "name", "score", "strand"]
    df = pd.read_csv(
        path, sep="\t", names=cols, header=None, comment="#",
        dtype={"chr": str},
    )
    n_cols = min(len(df.columns), len(cols))
    df.columns = cols[:n_cols]
    return df


def read_fasta(path):
    """Open a FASTA file return a pyfaidx.Fasta object for random access."""
    return pyfaidx.Fasta(path)
