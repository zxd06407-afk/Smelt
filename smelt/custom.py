"""Custom BED interval methylation aggregation."""
import pandas as pd
import numpy as np


def compute_custom(site_df, bed_df):
    """Compute methylation levels for BED-defined intervals.

    Args:
        site_df: DataFrame from compute_site_methylation()
        bed_df: DataFrame from read_bed() with 0-based intervals

    Returns:
        DataFrame: interval_id, chr, start, end, context, n_sites, meth, unmeth, ratio
    """
    results = []
    interval_id = 0

    for _, interval in bed_df.iterrows():
        chrom = interval["chr"]
        start, end = interval["start"], interval["end"]
        name = interval.get("name", f"interval_{interval_id}")
        interval_id += 1

        sites = site_df[
            (site_df["chr"] == chrom) &
            (site_df["pos"] >= start) &
            (site_df["pos"] < end)
        ]
        if sites.empty:
            continue

        for context in sites["context"].unique():
            ctx = sites[sites["context"] == context]
            meth = ctx["meth"].sum()
            unmet = ctx["unmeth"].sum()
            total = meth + unmet

            results.append({
                "interval_id": name,
                "chr": chrom,
                "start": start,
                "end": end,
                "context": context,
                "n_sites": len(ctx),
                "meth": meth,
                "unmeth": unmet,
                "ratio": meth / total if total > 0 else np.nan,
            })

    if not results:
        return pd.DataFrame(columns=[
            "interval_id", "chr", "start", "end", "context",
            "n_sites", "meth", "unmeth", "ratio",
        ])
    return pd.DataFrame(results)
