"""Gene body +/- flanking region metaplot computation."""
import pandas as pd
import numpy as np


def _make_bins(start, end, n_bins):
    """Split [start, end) into n_bins equal-width bins."""
    if n_bins == 0:
        return []
    width = (end - start) / n_bins
    return [(start + i * width, start + (i + 1) * width) for i in range(n_bins)]


def compute_metaplot(site_df, gene_intervals, context="CpG",
                      upstream=2000, downstream=2000,
                      body_bins=50, up_bins=20, down_bins=20):
    """Compute methylation signal across gene bodies and flanking regions.

    Args:
        site_df: DataFrame from compute_site_methylation()
        gene_intervals: DataFrame with chr, start, end, gene_id (0-based)
        context: Sequence context to analyze (default "CpG")
        upstream: bp upstream of gene start
        downstream: bp downstream of gene end
        body_bins: Number of equal-ratio bins for gene body
        up_bins: Number of bins for upstream region
        down_bins: Number of bins for downstream region

    Returns:
        DataFrame: region, bin, ratio (mean methylation across all genes)
    """
    sites = site_df[site_df["context"] == context].copy()
    results = []

    for _, gene in gene_intervals.iterrows():
        chrom = gene["chr"]
        g_start = gene["start"]   # 0-based TSS
        g_end = gene["end"]       # 0-based TES

        # Upstream bins (fixed bp)
        up_region = (g_start - upstream, g_start)
        for i, (b_start, b_end) in enumerate(_make_bins(up_region[0], up_region[1], up_bins)):
            if b_start < 0:
                continue
            in_bin = sites[(sites["chr"] == chrom) &
                           (sites["pos"] >= b_start) & (sites["pos"] < b_end)]
            if not in_bin.empty:
                results.append({
                    "gene_id": gene.get("gene_id"),
                    "region": "upstream",
                    "bin": i,
                    "ratio": in_bin["ratio"].mean(),
                })

        # Gene body bins (equal ratio)
        for i, (b_start, b_end) in enumerate(_make_bins(g_start, g_end, body_bins)):
            in_bin = sites[(sites["chr"] == chrom) &
                           (sites["pos"] >= b_start) & (sites["pos"] < b_end)]
            if not in_bin.empty:
                results.append({
                    "gene_id": gene.get("gene_id"),
                    "region": "body",
                    "bin": i,
                    "ratio": in_bin["ratio"].mean(),
                })

        # Downstream bins (fixed bp)
        for i, (b_start, b_end) in enumerate(_make_bins(g_end, g_end + downstream, down_bins)):
            in_bin = sites[(sites["chr"] == chrom) &
                           (sites["pos"] >= b_start) & (sites["pos"] < b_end)]
            if not in_bin.empty:
                results.append({
                    "gene_id": gene.get("gene_id"),
                    "region": "downstream",
                    "bin": i,
                    "ratio": in_bin["ratio"].mean(),
                })

    if not results:
        return pd.DataFrame(columns=["region", "bin", "ratio"])

    df = pd.DataFrame(results)
    output = df.groupby(["region", "bin"], as_index=False)["ratio"].mean()
    return output
