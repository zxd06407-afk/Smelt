"""File I/O for BISMARK, GTF, BED, and FASTA formats."""
import re
from collections import defaultdict

import pandas as pd
import pyfaidx

from smelt.utils import to_zero_based


def read_cx_report(path):
    """Read BISMARK coverage2cytosine CX_report.txt file.

    Format (1-based, tab-separated):
        chr pos strand meth unmeth context trinucleotide

    Context values from BISMARK: CG, CHH, CHG.

    Returns DataFrame with 0-based positions:
        chr, pos, strand, meth, unmeth, context, total, ratio
    """
    cols = ["chr", "pos_raw", "strand", "meth", "unmeth", "context", "trinucleotide"]
    df = pd.read_csv(
        path, sep="\t", names=cols, header=None,
        dtype={"chr": str, "strand": str, "context": str},
    )
    df["pos"] = df["pos_raw"].apply(to_zero_based)
    ctx_map = {"CG": "CpG"}
    df["context"] = df["context"].replace(ctx_map)
    df["total"] = df["meth"] + df["unmeth"]
    df["ratio"] = df["meth"] / df["total"]
    return df[["chr", "pos", "strand", "meth", "unmeth", "context", "total", "ratio"]]


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


XM_CONTEXT = {"Z": "CpG", "z": "CpG",
               "H": "CHH", "h": "CHH",
               "X": "CHG", "x": "CHG",
               "U": None, "u": None}


def read_bismark_sam(path):
    """Read BISMARK SAM/BAM file, extracting per-cytosine methylation counts.

    Parses the XM tag for context (Z=CpG, H=CHH, X=CHG) and case for
    methylation status (upper=meth, lower=unmeth).  Extracts strand from the
    alignment flag.

    Returns DataFrame:
        chr, pos (0-based), strand, meth, unmeth, context, total, ratio
    """
    import pysam

    # (chrom, pos, strand, context) -> [meth, unmeth]
    counts = defaultdict(lambda: [0, 0])

    mode = "rb" if str(path).endswith(".bam") else "r"
    with pysam.AlignmentFile(str(path), mode) as sam:
        for read in sam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if not read.has_tag("XM"):
                continue

            xm = read.get_tag("XM")
            chrom = sam.get_reference_name(read.reference_id)
            strand = "-" if read.is_reverse else "+"
            ref_positions = read.get_reference_positions(full_length=True)

            for read_idx, ref_pos in enumerate(ref_positions):
                if ref_pos is None or read_idx >= len(xm):
                    continue
                xm_char = xm[read_idx]
                if xm_char == ".":
                    continue

                ctx = XM_CONTEXT.get(xm_char.upper())
                if ctx is None:
                    continue

                is_meth = xm_char.isupper()
                if is_meth:
                    counts[(chrom, ref_pos, strand, ctx)][0] += 1
                else:
                    counts[(chrom, ref_pos, strand, ctx)][1] += 1

    rows = []
    for (chrom, pos, strand, context), (meth, unmeth) in counts.items():
        total = meth + unmeth
        rows.append({
            "chr": chrom,
            "pos": pos,
            "strand": strand,
            "context": context,
            "meth": meth,
            "unmeth": unmeth,
            "total": total,
            "ratio": meth / total if total > 0 else 0.0,
        })

    return pd.DataFrame(rows, columns=[
        "chr", "pos", "strand", "context", "meth", "unmeth", "total", "ratio",
    ])
