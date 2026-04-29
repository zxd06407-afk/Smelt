"""Fixed-size sliding window methylation aggregation."""
import pandas as pd
import numpy as np


def _make_windows(chrom, chrom_start, chrom_end, window_size=2000, step=500):
    """Generate window intervals for a chromosome.

    Returns list of (start, end) tuples, 0-based half-open.
    """
    windows = []
    for start in range(chrom_start, chrom_end, step):
        end = start + window_size
        windows.append((start, end))
    return windows


def _sites_in_window(sites, start, end):
    """Return sites within [start, end)."""
    return sites[(sites["pos"] >= start) & (sites["pos"] < end)]


def compute_windows(site_df, window_size=2000, step=500, min_sites=10):
    """Aggregate methylation into fixed-size sliding windows.

    Args:
        site_df: DataFrame from compute_site_methylation()
        window_size: Window size in bp (default 2000)
        step: Step size in bp (default 500)
        min_sites: Minimum number of cytosines per window per context

    Returns:
        DataFrame: chr, start, end, context, n_sites, meth, unmeth, total, ratio
    """
    results = []
    for chrom in site_df["chr"].unique():
        chrom_sites = site_df[site_df["chr"] == chrom]
        chrom_end = chrom_sites["pos"].max() + 1

        for start, end in _make_windows(chrom, 0, chrom_end, window_size, step):
            window_sites = _sites_in_window(chrom_sites, start, end)
            if window_sites.empty:
                continue

            for context in window_sites["context"].unique():
                ctx_sites = window_sites[window_sites["context"] == context]
                n = len(ctx_sites)
                meth_sum = ctx_sites["meth"].sum()
                un_sum = ctx_sites["unmeth"].sum()
                total = meth_sum + un_sum

                if n < min_sites:
                    continue

                results.append({
                    "chr": chrom,
                    "start": start,
                    "end": end,
                    "context": context,
                    "n_sites": n,
                    "meth": meth_sum,
                    "unmeth": un_sum,
                    "total": total,
                    "ratio": meth_sum / total if total > 0 else np.nan,
                })

    if not results:
        return pd.DataFrame(columns=[
            "chr", "start", "end", "context",
            "n_sites", "meth", "unmeth", "total", "ratio",
        ])

    return pd.DataFrame(results)
