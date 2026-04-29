"""Coordinate conversion and sequence-context utilities."""


def to_zero_based(pos, is_interval=False, end=None, source="default"):
    """Convert 1-based position/interval to 0-based half-open.

    Args:
        pos: 1-based position, or start of interval if is_interval=True
        is_interval: If True, treat pos as interval start and require end
        end: 1-based end of interval (only when is_interval=True)
        source: "bed" to skip conversion (BED is already 0-based)

    Returns:
        0-based position, or (start, end) tuple if is_interval=True
    """
    if source == "bed":
        if is_interval:
            return (pos, end)
        return pos

    if is_interval:
        if end is None:
            raise ValueError("end is required when is_interval=True")
        return (pos - 1, end)

    return pos - 1


def to_one_based(pos):
    """Convert 0-based position to 1-based."""
    return pos + 1


def classify_context(triplet):
    """Classify a cytosine context from the downstream 2 bases.

    Args:
        triplet: String of 3 bases starting with C (e.g. "CGA").
                 Length < 3 returns None.

    Returns:
        "CpG", "CHG", "CHH", or None if insufficient bases.
    """
    if len(triplet) < 2:
        return None
    triplet = triplet.upper()
    if triplet[1] == "G":
        return "CpG"
    if len(triplet) == 3 and triplet[2] == "G":
        return "CHG"
    return "CHH"


def merge_cpg_strands(df):
    """Merge CpG methylation counts from + and - strands.

    For CpG sites, met/unmet counts are summed across strands.
    CHH and CHG sites pass through unchanged.

    Expects columns: chr, pos, context, strand, meth, unmeth, total, ratio.
    """
    import pandas as pd

    df = df.copy()
    if "strand" not in df.columns:
        df["strand"] = "+"

    cpg = df[df["context"] == "CpG"]
    non_cpg = df[df["context"] != "CpG"]

    if cpg.empty:
        return df

    merged = cpg.groupby(["chr", "pos"], as_index=False).agg({
        "meth": "sum",
        "unmeth": "sum",
    })
    merged["context"] = "CpG"
    merged["strand"] = "+"
    merged["total"] = merged["meth"] + merged["unmeth"]
    merged["ratio"] = merged["meth"] / merged["total"]

    cols = ["chr", "pos", "context", "strand", "meth", "unmeth", "total", "ratio"]
    return pd.concat([merged[cols], non_cpg[cols]], ignore_index=True)
