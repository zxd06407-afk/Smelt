"""Site-level methylation computation."""
import pandas as pd
from smelt.utils import merge_cpg_strands
from smelt.filter import filter_by_depth


def compute_site_methylation(site_df, min_depth=5, merge_cpg=True,
                              sample_name=None, threads=1):
    """Compute per-cytosine methylation.

    Context and strand must already be present in the DataFrame
    (set by read_bismark_sam or read_cx_report).

    Args:
        site_df: DataFrame with chr, pos, strand, meth, unmeth, context columns
        min_depth: Minimum read depth filter
        merge_cpg: Merge CpG dyad counts across strands
        sample_name: Optional sample identifier added to output
        threads: Number of workers (unused, accepted for API compatibility)

    Returns:
        DataFrame: chr, pos, context, strand, meth, unmeth, total, ratio[, sample]
    """
    expected_cols = {"chr", "pos", "meth", "unmeth"}
    missing = expected_cols - set(site_df.columns)
    if missing:
        raise ValueError(f"Input missing required columns: {missing}")

    df = filter_by_depth(site_df, min_depth)

    if merge_cpg:
        df = merge_cpg_strands(df)

    if sample_name is not None:
        df["sample"] = sample_name

    return df.reset_index(drop=True)
