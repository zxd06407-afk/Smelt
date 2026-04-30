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
    """Merge CpG dyad partners on opposite strands.

    A CpG dyad consists of C on the + strand at position N and C on the
    - strand at position N+1.  These measure the same methylation event
    and are merged into one record at position N with combined counts.

    CHH and CHG sites pass through unchanged.

    Expects columns: chr, pos, context, strand, meth, unmeth, total, ratio.
    """
    import pandas as pd

    df = df.copy()
    if "strand" not in df.columns or "context" not in df.columns:
        return df

    cpg = df[df["context"] == "CpG"].copy()
    non_cpg = df[df["context"] != "CpG"].copy()

    if cpg.empty:
        return df

    cpg_plus = cpg[cpg["strand"] == "+"].copy()
    cpg_minus = cpg[cpg["strand"] == "-"].copy()

    if cpg_minus.empty:
        return pd.concat([cpg_plus, non_cpg], ignore_index=True)

    # Build index of minus-strand CpG by (chr, pos-1) for dyad lookup
    cpg_minus_idx = cpg_minus.set_index(["chr", "pos"])
    matched_minus_positions = set()

    rows = []
    for _, plus_row in cpg_plus.iterrows():
        chrom = plus_row["chr"]
        pos = plus_row["pos"]
        # Look for minus-strand partner at pos+1
        partner_key = (chrom, pos + 1)
        if partner_key in cpg_minus_idx.index:
            minus_row = cpg_minus_idx.loc[partner_key]
            # Handle case where multiple minus rows match
            if isinstance(minus_row, pd.DataFrame):
                minus_row = minus_row.iloc[0]
            rows.append({
                "chr": chrom, "pos": pos,
                "context": "CpG", "strand": "+",
                "meth": plus_row["meth"] + minus_row["meth"],
                "unmeth": plus_row["unmeth"] + minus_row["unmeth"],
            })
            matched_minus_positions.add(partner_key)
        else:
            rows.append({
                "chr": chrom, "pos": pos,
                "context": "CpG", "strand": "+",
                "meth": plus_row["meth"],
                "unmeth": plus_row["unmeth"],
            })

    # Add unmatched minus-strand CpG sites
    for idx, minus_row in cpg_minus.iterrows():
        if (minus_row["chr"], minus_row["pos"]) not in matched_minus_positions:
            rows.append({
                "chr": minus_row["chr"], "pos": minus_row["pos"],
                "context": "CpG", "strand": "+",
                "meth": minus_row["meth"],
                "unmeth": minus_row["unmeth"],
            })

    merged = pd.DataFrame(rows)
    merged["total"] = merged["meth"] + merged["unmeth"]
    merged["ratio"] = merged["meth"] / merged["total"]

    cols = ["chr", "pos", "context", "strand", "meth", "unmeth", "total", "ratio"]
    return pd.concat([merged[cols], non_cpg[cols]], ignore_index=True)


def parallel_chromosomes(df, func, threads=1, **kwargs):
    """Process each chromosome in parallel using multiprocessing.

    Args:
        df: DataFrame with a 'chr' column.
        func: Function to apply per chromosome, called as func(chrom_df, **kwargs).
        threads: Number of worker processes (1 = sequential).
        **kwargs: Additional arguments passed to func.

    Returns:
        Concatenated DataFrame of results from all chromosomes.
    """
    import pandas as pd
    from concurrent.futures import ProcessPoolExecutor

    chromosomes = sorted(df["chr"].unique())
    if threads <= 1 or len(chromosomes) <= 1:
        results = []
        for chrom in chromosomes:
            chrom_df = df[df["chr"] == chrom].copy()
            results.append(func(chrom_df, **kwargs))
        if results:
            return pd.concat(results, ignore_index=True)
        return pd.DataFrame()

    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {}
        for chrom in chromosomes:
            chrom_df = df[df["chr"] == chrom].copy()
            futures[executor.submit(func, chrom_df, **kwargs)] = chrom

        results = []
        for future in futures:
            results.append(future.result())

    if results:
        return pd.concat(results, ignore_index=True)
    return pd.DataFrame()
