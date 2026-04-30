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


def read_bismark_sam(path):
    """Read BISMARK SAM/BAM file, extracting per-cytosine methylation counts.

    Parses the XM tag (methylation call string) from each alignment.
    Uppercase characters (Z/H/X) count as methylated, lowercase (z/h/x) as
    unmethylated.  Positions with '.' or missing XM are skipped.

    Returns a DataFrame compatible with read_cov() output:
        chr, pos (0-based), meth, unmeth, total, ratio
    """
    import pysam
    from collections import defaultdict

    mode = "rb" if str(path).endswith(".bam") else "r"
    counts = defaultdict(lambda: [0, 0])  # (meth, unmeth)

    with pysam.AlignmentFile(str(path), mode) as sam:
        for read in sam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if not read.has_tag("XM"):
                continue

            xm = read.get_tag("XM")
            chrom = sam.get_reference_name(read.reference_id)
            ref_positions = read.get_reference_positions(full_length=True)

            for read_idx, ref_pos in enumerate(ref_positions):
                if ref_pos is None or read_idx >= len(xm):
                    continue
                xm_char = xm[read_idx]
                if xm_char == ".":
                    continue

                if xm_char.isupper():
                    counts[(chrom, ref_pos)][0] += 1
                elif xm_char.islower():
                    counts[(chrom, ref_pos)][1] += 1

    rows = []
    for (chrom, pos), (meth, unmeth) in counts.items():
        total = meth + unmeth
        rows.append({
            "chr": chrom,
            "pos": pos,  # pysam returns 0-based positions
            "meth": meth,
            "unmeth": unmeth,
            "total": total,
            "ratio": meth / total if total > 0 else 0.0,
        })

    return pd.DataFrame(rows, columns=["chr", "pos", "meth", "unmeth", "total", "ratio"])
