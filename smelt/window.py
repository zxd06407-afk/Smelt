"""Fixed-size sliding window methylation aggregation."""
import pandas as pd
import numpy as np
from smelt.utils import parallel_chromosomes


def _window_chromosome(chrom_sites, window_size=2000, step=500, min_sites=10):
    """Compute sliding windows for a single chromosome using binary search.

    Uses numpy.searchsorted to find window boundaries in O(log n) per window
    instead of O(n) boolean indexing.
    """
    results = []
    chrom = chrom_sites["chr"].iloc[0]
    chrom_end = chrom_sites["pos"].max() + 1

    for context in chrom_sites["context"].unique():
        ctx_sites = chrom_sites[chrom_sites["context"] == context].sort_values("pos")
        if ctx_sites.empty:
            continue
        positions = ctx_sites["pos"].values

        for start in range(0, chrom_end, step):
            end = start + window_size
            left = np.searchsorted(positions, start, side="left")
            right = np.searchsorted(positions, end, side="left")
            n = right - left
            if n < min_sites:
                continue

            in_window = ctx_sites.iloc[left:right]
            meth_sum = in_window["meth"].sum()
            un_sum = in_window["unmeth"].sum()
            total = meth_sum + un_sum
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

    return pd.DataFrame(results)


def compute_windows(site_df, window_size=2000, step=500, min_sites=10,
                     sample_name=None, threads=1):
    """Aggregate methylation into fixed-size sliding windows.

    Args:
        site_df: DataFrame from compute_site_methylation()
        window_size: Window size in bp (default 2000)
        step: Step size in bp (default 500)
        min_sites: Minimum number of cytosines per window per context
        sample_name: Optional sample identifier added to output
        threads: Number of worker processes for per-chromosome parallelism

    Returns:
        DataFrame: chr, start, end, context, n_sites, meth, unmeth, total, ratio
    """
    result = parallel_chromosomes(
        site_df, _window_chromosome, threads=threads,
        window_size=window_size, step=step, min_sites=min_sites,
    )
    if result.empty:
        result = pd.DataFrame(columns=[
            "chr", "start", "end", "context",
            "n_sites", "meth", "unmeth", "total", "ratio",
        ])
    if sample_name is not None:
        result["sample"] = sample_name
    return result
