"""Genome and chromosome-level methylation statistics."""
import pandas as pd


def compute_stats(df, sample_name=None):
    """Compute summary statistics per chromosome and context, plus genome-wide.

    Args:
        df: Site or window DataFrame with chr, context, ratio, total columns
        sample_name: Optional sample identifier added to output

    Returns:
        DataFrame: chr, context, mean_ratio, n_sites, mean_depth
                   Includes a chr="genome" row per context.
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

    # Append genome-wide rows (weighted by n_sites)
    genome_rows = []
    for ctx in grouped["context"].unique():
        ctx_df = grouped[grouped["context"] == ctx]
        total_sites = ctx_df["n_sites"].sum()
        if total_sites > 0:
            w_mean_ratio = (ctx_df["mean_ratio"] * ctx_df["n_sites"]).sum() / total_sites
            w_mean_depth = (ctx_df["mean_depth"] * ctx_df["n_sites"]).sum() / total_sites
        else:
            w_mean_ratio = 0.0
            w_mean_depth = 0.0
        genome_rows.append({
            "chr": "genome",
            "context": ctx,
            "mean_ratio": w_mean_ratio,
            "n_sites": total_sites,
            "mean_depth": w_mean_depth,
        })
    genome_df = pd.DataFrame(genome_rows)
    grouped = pd.concat([grouped, genome_df], ignore_index=True)

    grouped = grouped.sort_values(["chr", "context"]).reset_index(drop=True)

    if sample_name is not None:
        grouped["sample"] = sample_name

    return grouped
