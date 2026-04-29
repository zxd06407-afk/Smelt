"""Genome and chromosome-level methylation statistics."""
import pandas as pd


def compute_stats(df):
    """Compute summary statistics per chromosome and context.

    Args:
        df: Site or window DataFrame with chr, context, ratio, total columns

    Returns:
        DataFrame: chr, context, mean_ratio, n_sites, mean_depth
    """
    if "pos" in df.columns:
        site_col = "pos"
    else:
        site_col = "n_sites"

    if site_col == "pos":
        grouped = df.groupby(["chr", "context"], as_index=False).agg(
            mean_ratio=("ratio", "mean"),
            n_sites=("pos", "count"),
            mean_depth=("total", "mean"),
        )
    else:
        grouped = df.groupby(["chr", "context"], as_index=False).agg(
            mean_ratio=("ratio", "mean"),
            n_sites=("n_sites", "sum"),
            mean_depth=("total", "mean"),
        )

    return grouped.sort_values(["chr", "context"]).reset_index(drop=True)
