import pandas as pd
import pytest
import numpy as np
from smelt.dmr import call_dmr, _bh_correction


def make_known_df():
    """Constructed data with known meth/unmeth counts for precise assertions.

    Window 1: group1=(80,20), group2=(20,80) → hyper DMR (delta=-0.6)
    Window 2: group1=(50,50), group2=(50,50) → no difference (delta=0)
    Window 3: group1=(10,90), group2=(90,10) → hyper DMR (delta=+0.8)
    Window 4: group1=(40,60), group2=(60,40) → weak signal (delta=+0.2)
    Each window has 2 samples per group with equal split of counts.
    """
    rows = []
    for sample_idx, (sample, g1_meth, g1_unmeth, g2_meth, g2_unmeth) in enumerate([
        ("ctrl_1", 80, 20, 20, 80),
        ("ctrl_2", 80, 20, 20, 80),
        ("case_1", 20, 80, 80, 20),
        ("case_2", 20, 80, 80, 20),
    ]):
        for i, (meth_base, unmeth_base) in enumerate([
            (g1_meth, g1_unmeth),  # window 1: strong hypo
            (50, 50),               # window 2: no diff
            (g2_meth, g2_unmeth),  # window 3: strong hyper
            (40, 60),               # window 4: weak (delta=0.2 for case vs ctrl)
        ]):
            meth = meth_base if sample.startswith("ctrl") else (
                20 if i == 0 else (50 if i == 1 else (90 if i == 2 else 60))
            )
            unmeth = unmeth_base if sample.startswith("ctrl") else (
                80 if i == 0 else (50 if i == 1 else (10 if i == 2 else 40))
            )
            rows.append({
                "sample": sample,
                "chr": "chrA",
                "start": i * 2000,
                "end": (i + 1) * 2000,
                "context": "CpG",
                "n_sites": 10,
                "meth": meth,
                "unmeth": unmeth,
                "total": meth + unmeth,
                "ratio": meth / (meth + unmeth),
            })
    return pd.DataFrame(rows)


def test_bh_correction():
    p = np.array([0.001, 0.01, 0.05, 0.5, np.nan])
    q = _bh_correction(p)
    assert len(q) == 5
    assert np.isnan(q[4])
    assert q[0] <= q[1]
    assert all(q[~np.isnan(q)] >= 0)
    assert all(q[~np.isnan(q)] <= 1)


def test_bh_correction_monotonic():
    """BH q-values should be monotonic with p-values after sorting."""
    p = np.array([0.001, 0.05, 0.01, 0.5, 0.1])
    q = _bh_correction(p)
    # Sort both; q should be in same order
    order = np.argsort(p)
    q_sorted = q[order]
    for i in range(len(q_sorted) - 1):
        assert q_sorted[i] <= q_sorted[i + 1] + 1e-10


def test_call_dmr_output_columns():
    df = make_known_df()
    groups = {"ctrl": ["ctrl_1", "ctrl_2"], "case": ["case_1", "case_2"]}
    result, _ = call_dmr(df, groups, group1="ctrl", group2="case")
    for col in ["chr", "start", "end", "context", "p_value", "q_value",
                "delta", "direction", "ratio_group1", "ratio_group2"]:
        assert col in result.columns


def test_call_dmr_with_known_data():
    """Precise assertions on constructed data with known counts."""
    df = make_known_df()
    groups = {"ctrl": ["ctrl_1", "ctrl_2"], "case": ["case_1", "case_2"]}

    # Include all windows including noise (q_threshold=1, delta_threshold=-1)
    result, excl = call_dmr(df, groups, group1="ctrl", group2="case",
                             q_threshold=1.0, delta_threshold=-1.0)

    # Window 1: ctrl(160,40) vs case(40,160) → strong hypo, delta ≈ -0.6
    w1 = result[result["start"] == 0]
    assert len(w1) == 1
    assert w1.iloc[0]["direction"] == "hypo"
    assert w1.iloc[0]["delta"] < -0.5
    assert w1.iloc[0]["p_value"] < 1e-10

    # Window 3: ctrl(20,180) vs case(180,20) → strong hyper, delta ≈ +0.8
    w3 = result[result["start"] == 4000]
    assert len(w3) == 1
    assert w3.iloc[0]["direction"] == "hyper"
    assert w3.iloc[0]["delta"] >= 0.7

    # Window 2: delta=0, q=1.0, excluded with no_signal
    w2_excl = excl[(excl["start"] == 2000) & (excl["context"] == "CpG")]
    assert len(w2_excl) == 1
    assert w2_excl.iloc[0]["category"] == "no_signal"

    # Window 4: weak signal (delta=0.2) passes threshold
    w4 = result[result["start"] == 6000]
    assert len(w4) == 1
    assert abs(w4.iloc[0]["delta"]) > 0.15


def test_call_dmr_filters_by_thresholds():
    """With strict thresholds, only strong DMRs pass."""
    df = make_known_df()
    groups = {"ctrl": ["ctrl_1", "ctrl_2"], "case": ["case_1", "case_2"]}
    result, _ = call_dmr(df, groups, group1="ctrl", group2="case",
                          q_threshold=0.05, delta_threshold=0.3)
    # Windows 1 (hypo, delta≈-0.6) and 3 (hyper, delta≈+0.8) pass
    # Windows 2 (no diff) and 4 (weak, delta≈+0.2) excluded
    starts = set(result["start"])
    assert 0 in starts       # window 1
    assert 2000 not in starts  # window 2 excluded
    assert 4000 in starts    # window 3
    assert 6000 not in starts  # window 4 excluded


def test_call_dmr_excluded_categories():
    """Verify each exclusion category appears correctly."""
    df = make_known_df()
    groups = {"ctrl": ["ctrl_1", "ctrl_2"], "case": ["case_1", "case_2"]}

    _, excl = call_dmr(df, groups, group1="ctrl", group2="case",
                        q_threshold=0.05, delta_threshold=0.3)
    assert "category" in excl.columns
    if len(excl) > 0:
        valid = {"low_coverage", "all_zero", "fisher_failed", "no_signal"}
        assert all(c in valid for c in excl["category"].unique())
        # Windows 2 and 4 should be "no_signal"
        ns = excl[excl["category"] == "no_signal"]
        assert len(ns) >= 2


def test_call_dmr_per_context_fdr_independent():
    """BH correction runs independently per context."""
    rows = []
    # 3 CpG windows with various p-values (will get BH corrected together)
    # 2 CHH windows (will get their own BH correction)
    for ctx, base in [("CpG", 80), ("CHH", 20)]:
        for i in range(3 if ctx == "CpG" else 2):
            rows.append({
                "sample": "ctrl", "chr": "chrA",
                "start": i * 2000, "end": (i + 1) * 2000,
                "context": ctx, "n_sites": 10,
                "meth": base, "unmeth": 20,
                "total": base + 20, "ratio": base / (base + 20),
            })
            rows.append({
                "sample": "case", "chr": "chrA",
                "start": i * 2000, "end": (i + 1) * 2000,
                "context": ctx, "n_sites": 10,
                "meth": 20 + i * 10, "unmeth": 80 - i * 10,
                "total": 100, "ratio": (20 + i * 10) / 100,
            })
    df = pd.DataFrame(rows)
    groups = {"ctrl": ["ctrl"], "case": ["case"]}
    result, _ = call_dmr(df, groups, group1="ctrl", group2="case",
                          q_threshold=1.0, delta_threshold=0.0,
                          min_samples_per_group=1)
    # Both contexts should appear (each independently corrected)
    assert "CpG" in set(result["context"])
    assert "CHH" in set(result["context"])
    # CpG windows should have valid q-values (not all 1.0)
    cpg = result[result["context"] == "CpG"]
    assert any(cpg["q_value"] < 1.0)
