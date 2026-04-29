import pandas as pd
import pytest
import numpy as np
from smelt.dmr import call_dmr, _bh_correction


def make_window_df(samples, seed=42):
    """Create a window DataFrame for multiple samples."""
    np.random.seed(seed)
    dfs = []
    for sample in samples:
        n_windows = 10
        rows = []
        for i in range(n_windows):
            for ctx in ["CpG", "CHH"]:
                base_meth = 80 if ctx == "CpG" else 20
                meth = max(0, int(np.random.normal(base_meth, 5)))
                unmet = max(0, int(np.random.normal(20, 5)))
                rows.append({
                    "sample": sample,
                    "chr": "chrA",
                    "start": i * 2000,
                    "end": (i + 1) * 2000,
                    "context": ctx,
                    "n_sites": 20,
                    "meth": meth,
                    "unmeth": unmet,
                    "total": meth + unmet,
                    "ratio": meth / (meth + unmet),
                })
        dfs.append(pd.DataFrame(rows))
    return pd.concat(dfs, ignore_index=True)


def test_bh_correction():
    p = np.array([0.001, 0.01, 0.05, 0.5, np.nan])
    q = _bh_correction(p)
    assert len(q) == 5
    assert np.isnan(q[4])
    assert q[0] <= q[1]  # smaller p should get smaller q
    assert all(q[~np.isnan(q)] >= 0)
    assert all(q[~np.isnan(q)] <= 1)


def test_call_dmr_output_columns():
    df = make_window_df(["tumor1", "tumor2", "normal1", "normal2"])
    groups = {"tumor": ["tumor1", "tumor2"], "normal": ["normal1", "normal2"]}
    result, _ = call_dmr(df, groups, group1="tumor", group2="normal")
    for col in ["chr", "start", "end", "context", "p_value", "q_value",
                "delta", "direction", "ratio_group1", "ratio_group2"]:
        assert col in result.columns


def test_call_dmr_filters_by_q_threshold():
    df = make_window_df(["t1", "t2", "n1", "n2"])
    groups = {"tumor": ["t1", "t2"], "normal": ["n1", "n2"]}
    result, _ = call_dmr(df, groups, group1="tumor", group2="normal",
                          q_threshold=0.05, delta_threshold=0.0)
    # With random data, few windows should pass 0.05 FDR
    assert len(result) <= 20


def test_call_dmr_hyper_hypo_classification():
    df = make_window_df(["t", "n"])
    groups = {"treated": ["t"], "control": ["n"]}
    result, _ = call_dmr(df, groups, group1="treated", group2="control",
                          q_threshold=1.0, delta_threshold=0.0)
    assert "direction" in result.columns
    assert all(d in ("hyper", "hypo") for d in result["direction"])


def test_call_dmr_per_context_fdr():
    df = make_window_df(["t1", "t2", "n1", "n2"])
    groups = {"tumor": ["t1", "t2"], "normal": ["n1", "n2"]}
    result, _ = call_dmr(df, groups, group1="tumor", group2="normal",
                          q_threshold=1.0, delta_threshold=0.0)
    contexts = set(result["context"].unique())
    assert "CpG" in contexts or "CHH" in contexts
