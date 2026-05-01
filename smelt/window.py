"""Fixed-size sliding window methylation aggregation."""
import sys
import pandas as pd
import numpy as np
from smelt.utils import parallel_chromosomes


def _window_chromosome(chrom_sites, window_size=2000, step=500, min_sites=10):
    """Compute sliding windows for a single chromosome.

    Uses two-pointer O(n+m) algorithm: sites and windows are sorted,
    pointers only advance forward.
    """
    chrom = str(chrom_sites["chr"].iloc[0])
    sites = chrom_sites.sort_values("pos").reset_index(drop=True)
    n_sites = len(sites)
    chrom_end = sites["pos"].max() + 1

    results = []
    left = 0   # first site index >= window start

    for start in range(0, chrom_end, step):
        end = start + window_size

        # Advance left pointer to first site >= start
        while left < n_sites and sites["pos"].iloc[left] < start:
            left += 1

        # Advance right pointer to first site >= end
        right = left
        while right < n_sites and sites["pos"].iloc[right] < end:
            right += 1

        if left == right:
            continue  # empty window

        window_sites = sites.iloc[left:right]

        for context in ("CpG", "CHH", "CHG"):
            ctx_sites = window_sites[window_sites["context"] == context]
            n = len(ctx_sites)
            if n < min_sites:
                continue

            meth_sum = ctx_sites["meth"].sum()
            un_sum = ctx_sites["unmeth"].sum()
            total = meth_sum + un_sum

            results.append({
                "chr": chrom, "start": start, "end": end,
                "context": context, "n_sites": n,
                "meth": meth_sum, "unmeth": un_sum,
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
