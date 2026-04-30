"""GTF feature-type methylation aggregation."""
import pandas as pd
import numpy as np


def compute_elements(site_df, gtf_df, features=None, sample_name=None):
    """Compute methylation levels for GTF feature types.

    Args:
        site_df: DataFrame from compute_site_methylation() with 0-based positions
        gtf_df: DataFrame from read_gtf() with 0-based intervals
        features: List of feature types to include (e.g. ["gene", "exon"]).
                  None includes all.
        sample_name: Optional sample identifier added to output

    Returns:
        DataFrame: chr, feature_type, feature_id, context, n_sites, meth, unmeth, ratio
    """
    gtf = gtf_df.copy()
    if features is not None:
        gtf = gtf[gtf["feature"].isin(features)]

    results = []
    for _, feat in gtf.iterrows():
        chrom = feat["chr"]
        start, end = feat["start"], feat["end"]
        feat_type = feat["feature"]
        feat_id = feat.get("gene_id", None)

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
                "chr": chrom,
                "feature_type": feat_type,
                "feature_id": feat_id,
                "context": context,
                "n_sites": len(ctx),
                "meth": meth,
                "unmeth": unmet,
                "ratio": meth / total if total > 0 else np.nan,
            })

    if not results:
        result = pd.DataFrame(columns=[
            "chr", "feature_type", "feature_id", "context",
            "n_sites", "meth", "unmeth", "ratio",
        ])
    else:
        result = pd.DataFrame(results)
    if sample_name is not None:
        result["sample"] = sample_name
    return result
