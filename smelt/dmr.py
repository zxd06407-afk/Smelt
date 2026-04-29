"""Window-level differential methylation region calling."""
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact


def _fisher_test_window(group1_df, group2_df):
    """Run Fisher's exact test on a single window x context combination.

    Args:
        group1_df: Subset of window data for group 1 (one row per sample)
        group2_df: Subset of window data for group 2

    Returns:
        p-value (float) or NaN if test cannot be performed
    """
    meth1 = group1_df["meth"].sum()
    unmeth1 = group1_df["unmeth"].sum()
    meth2 = group2_df["meth"].sum()
    unmeth2 = group2_df["unmeth"].sum()

    total1 = meth1 + unmeth1
    total2 = meth2 + unmeth2

    if total1 == 0 and total2 == 0:
        return np.nan

    table = np.array([[meth1, unmeth1], [meth2, unmeth2]])
    _, p = fisher_exact(table)
    return p


def _bh_correction(p_values):
    """Benjamini-Hochberg FDR correction.

    Args:
        p_values: Array of p-values

    Returns:
        Array of q-values (adjusted p-values)
    """
    p = np.array(p_values, dtype=float)
    n = len(p)
    if n == 0:
        return np.array([])

    mask = ~np.isnan(p)
    valid_p = p[mask]
    n_valid = len(valid_p)

    if n_valid == 0:
        return np.full(n, np.nan)

    ranks = np.argsort(np.argsort(valid_p)) + 1
    q_valid = np.minimum(1, valid_p * n_valid / ranks)
    sorted_order = np.argsort(valid_p)
    q_valid_sorted = q_valid[sorted_order]
    for i in range(n_valid - 2, -1, -1):
        q_valid_sorted[i] = min(q_valid_sorted[i], q_valid_sorted[i + 1])
    q_valid = q_valid_sorted[np.argsort(sorted_order)]

    q_all = np.full(n, np.nan)
    q_all[mask] = q_valid
    return q_all


def call_dmr(window_df, sample_groups, group1, group2,
             q_threshold=0.05, delta_threshold=0.2,
             min_samples_per_group=2, fdr_by_context=True):
    """Call differentially methylated regions between two sample groups.

    Args:
        window_df: DataFrame with window data. Must contain a 'sample' column.
        sample_groups: Dict mapping group name -> list of sample names.
        group1: Name of reference group (e.g. "normal")
        group2: Name of comparison group (e.g. "tumor")
        q_threshold: FDR q-value cutoff (default 0.05)
        delta_threshold: Minimum absolute methylation difference (default 0.2)
        min_samples_per_group: Minimum samples per group with coverage
        fdr_by_context: If True, BH correction per sequence context

    Returns:
        Tuple of (dmr_df, excluded_df)
    """
    samples1 = sample_groups[group1]
    samples2 = sample_groups[group2]

    results = []
    excluded = []
    window_keys = ["chr", "start", "end", "context"]

    for (chrom, start, end, context), window_group in window_df.groupby(window_keys):
        g1 = window_group[window_group["sample"].isin(samples1)]
        g2 = window_group[window_group["sample"].isin(samples2)]

        n1 = len(g1[g1["total"] > 0])
        n2 = len(g2[g2["total"] > 0])

        if n1 < min_samples_per_group or n2 < min_samples_per_group:
            excluded.append({
                "chr": chrom, "start": start, "end": end,
                "context": context,
                "reason": f"insufficient samples: g1={n1}, g2={n2}",
            })
            continue

        total1 = g1["total"].sum()
        total2 = g2["total"].sum()
        if total1 == 0 and total2 == 0:
            excluded.append({
                "chr": chrom, "start": start, "end": end,
                "context": context, "reason": "all zero coverage",
            })
            continue

        p_value = _fisher_test_window(g1, g2)
        if np.isnan(p_value):
            excluded.append({
                "chr": chrom, "start": start, "end": end,
                "context": context, "reason": "fisher test failed",
            })
            continue

        ratio1 = g1["meth"].sum() / total1 if total1 > 0 else np.nan
        ratio2 = g2["meth"].sum() / total2 if total2 > 0 else np.nan

        results.append({
            "chr": chrom, "start": start, "end": end, "context": context,
            "p_value": p_value,
            "ratio_group1": ratio1,
            "ratio_group2": ratio2,
        })

    if not results:
        empty_dmr = pd.DataFrame(columns=[
            "chr", "start", "end", "context",
            "p_value", "q_value", "delta", "direction",
            "ratio_group1", "ratio_group2",
        ])
        empty_excl = pd.DataFrame(columns=["chr", "start", "end", "context", "reason"])
        return empty_dmr, empty_excl

    result_df = pd.DataFrame(results)

    if fdr_by_context:
        q_values = []
        for ctx in result_df["context"].unique():
            mask = result_df["context"] == ctx
            q_values.extend(_bh_correction(result_df.loc[mask, "p_value"].values))
        result_df["q_value"] = q_values
    else:
        result_df["q_value"] = _bh_correction(result_df["p_value"].values)

    result_df["delta"] = result_df["ratio_group2"] - result_df["ratio_group1"]
    result_df["direction"] = result_df["delta"].apply(
        lambda d: "hyper" if d > 0 else "hypo"
    )

    dmr_df = result_df[
        (result_df["q_value"] < q_threshold) &
        (result_df["delta"].abs() > delta_threshold)
    ].copy()

    remaining = result_df[
        ~((result_df["q_value"] < q_threshold) &
          (result_df["delta"].abs() > delta_threshold))
    ]
    for _, row in remaining.iterrows():
        excluded.append({
            "chr": row["chr"], "start": row["start"],
            "end": row["end"], "context": row["context"],
            "reason": "thresholds not met",
        })

    excluded_df = pd.DataFrame(excluded)
    return dmr_df.reset_index(drop=True), excluded_df.reset_index(drop=True)
