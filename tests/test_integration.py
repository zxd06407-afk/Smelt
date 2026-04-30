"""End-to-end pipeline integration tests using simulated CX_report data."""

import pytest
import pandas as pd
from pathlib import Path

from smelt.io import read_cx_report, read_gtf, read_bed
from smelt.site import compute_site_methylation
from smelt.window import compute_windows
from smelt.element import compute_elements
from smelt.metaplot import compute_metaplot
from smelt.custom import compute_custom
from smelt.dmr import call_dmr
from smelt.stats import compute_stats

DATA_DIR = Path("tests/data/integration")


@pytest.fixture(scope="module")
def all_sites():
    """Load pre-computed CX_report files for all 4 samples."""
    sites = {}
    for sample in ["sample_1", "sample_2", "sample_3", "sample_4"]:
        path = DATA_DIR / f"{sample}.CX_report.txt"
        if not path.exists():
            pytest.skip(f"Missing {path}. Run scripts/generate_test_data.py first.")
        df = read_cx_report(str(path))
        sites[sample] = compute_site_methylation(df, min_depth=5, merge_cpg=True,
                                                  sample_name=sample)
    return sites


@pytest.fixture(scope="module")
def all_windows(all_sites):
    windows = {}
    for sample, sites in all_sites.items():
        windows[sample] = compute_windows(sites, window_size=2000, step=500, min_sites=5,
                                           sample_name=sample)
    return windows


class TestIntegrationSite:
    def test_all_contexts_present(self, all_sites):
        for sample, df in all_sites.items():
            contexts = set(df["context"].unique())
            assert contexts.issubset({"CpG", "CHH", "CHG"}), f"{sample}: {contexts}"
            assert len(contexts) > 0

    def test_ratio_bounds(self, all_sites):
        for sample, df in all_sites.items():
            assert all(df["ratio"] >= 0), f"{sample} has negative ratios"
            assert all(df["ratio"] <= 1), f"{sample} has ratios > 1"

    def test_depth_filter_applied(self, all_sites):
        for sample, df in all_sites.items():
            assert all(df["total"] >= 5), f"{sample} has sites below min_depth"

    def test_sample_column_present(self, all_sites):
        for sample, df in all_sites.items():
            assert "sample" in df.columns
            assert all(df["sample"] == sample)


class TestIntegrationWindow:
    def test_windows_cover_genome(self, all_windows):
        for sample, df in all_windows.items():
            assert len(df) > 0, f"{sample} has no windows"

    def test_window_size_respected(self, all_windows):
        for sample, df in all_windows.items():
            sizes = (df["end"] - df["start"]).unique()
            assert set(sizes) == {2000}, f"{sample}: window sizes not all 2000"


class TestIntegrationElement:
    def test_elements_computed(self, all_sites):
        gtf_path = DATA_DIR / "test_genes.gtf"
        if not gtf_path.exists():
            pytest.skip("GTF not found")
        gtf = read_gtf(str(gtf_path))
        df = compute_elements(all_sites["sample_1"], gtf, features=["gene"])
        assert len(df) > 0

    def test_element_feature_types(self, all_sites):
        gtf_path = DATA_DIR / "test_genes.gtf"
        if not gtf_path.exists():
            pytest.skip("GTF not found")
        gtf = read_gtf(str(gtf_path))
        df = compute_elements(all_sites["sample_1"], gtf, features=["exon"])
        assert all(df["feature_type"] == "exon")


class TestIntegrationDMR:
    def test_dmr_output_has_required_columns(self, all_windows):
        dfs = []
        for sample, df in all_windows.items():
            df = df.copy()
            df["sample"] = sample
            dfs.append(df)
        merged = pd.concat(dfs, ignore_index=True)
        groups = {"group1": ["sample_1", "sample_2"],
                   "group2": ["sample_3", "sample_4"]}
        dmr_df, excl_df = call_dmr(
            merged, groups, group1="group1", group2="group2",
            q_threshold=0.05, delta_threshold=0.15,
        )
        assert "q_value" in dmr_df.columns
        assert "delta" in dmr_df.columns

    def test_dmr_has_direction(self, all_windows):
        dfs = []
        for sample, df in all_windows.items():
            df = df.copy()
            df["sample"] = sample
            dfs.append(df)
        merged = pd.concat(dfs, ignore_index=True)
        groups = {"group1": ["sample_1", "sample_2"],
                   "group2": ["sample_3", "sample_4"]}
        dmr_df, _ = call_dmr(merged, groups, group1="group1", group2="group2")
        if len(dmr_df) > 0:
            assert all(d in ("hyper", "hypo") for d in dmr_df["direction"])

    def test_dmr_excluded_has_category(self, all_windows):
        dfs = []
        for sample, df in all_windows.items():
            df = df.copy()
            df["sample"] = sample
            dfs.append(df)
        merged = pd.concat(dfs, ignore_index=True)
        groups = {"group1": ["sample_1", "sample_2"],
                   "group2": ["sample_3", "sample_4"]}
        _, excl_df = call_dmr(merged, groups, group1="group1", group2="group2")
        assert "category" in excl_df.columns
        if len(excl_df) > 0:
            valid = {"low_coverage", "all_zero", "fisher_failed", "no_signal"}
            assert all(c in valid for c in excl_df["category"].unique())


class TestIntegrationStats:
    def test_stats_per_context(self, all_sites):
        df = compute_stats(all_sites["sample_1"])
        assert len(df) > 0
        assert "CpG" in set(df["context"])

    def test_stats_has_genome_row(self, all_sites):
        df = compute_stats(all_sites["sample_1"])
        assert "genome" in set(df["chr"].unique())

    def test_stats_all_chromosomes(self, all_sites):
        df = compute_stats(all_sites["sample_1"])
        chromosomes = set(df[df["chr"] != "genome"]["chr"].unique())
        assert chromosomes == {"chrA", "chrB"}


EMBEDDED_DMRS = [
    ("chrA", 10000, 18000, "CpG", 0.7, "hyper"),
    ("chrA", 30000, 38000, "CpG", 0.5, "hypo"),
    ("chrA", 50000, 58000, "CpG", 0.4, "hyper"),
    ("chrA", 70000, 78000, "CpG", 0.3, "hypo"),
    ("chrA", 90000, 98000, "CpG", 0.2, "hyper"),
    ("chrB", 10000, 18000, "CHH", 0.6, "hyper"),
    ("chrB", 30000, 38000, "CHH", 0.5, "hypo"),
    ("chrB", 50000, 58000, "CHG", 0.5, "hyper"),
]


class TestIntegrationAccuracy:
    """Ground-truth assertions using 8 embedded DMRs in simulated data."""

    def test_all_eight_dmrs_detected(self, all_windows):
        """At q<0.05, |delta|>0.15, all 8 embedded DMRs are detected."""
        dfs = []
        for sample, df in all_windows.items():
            df = df.copy()
            df["sample"] = sample
            dfs.append(df)
        merged = pd.concat(dfs, ignore_index=True)
        groups = {"group1": ["sample_1", "sample_2"],
                   "group2": ["sample_3", "sample_4"]}
        dmr_df, _ = call_dmr(merged, groups, group1="group1", group2="group2",
                              q_threshold=0.05, delta_threshold=0.15)
        hits = set()
        for _, dmr in dmr_df.iterrows():
            for i, ed in enumerate(EMBEDDED_DMRS):
                if (dmr["chr"] == ed[0] and dmr["context"] == ed[3] and
                    dmr["start"] < ed[2] and dmr["end"] > ed[1]):
                    hits.add(i)
        assert len(hits) >= 8, f"Only {len(hits)}/8 DMRs detected"

    def test_dmr_direction_matches_preset(self, all_windows):
        """Each detected DMR region has correct hyper/hypo direction."""
        dfs = []
        for sample, df in all_windows.items():
            df = df.copy()
            df["sample"] = sample
            dfs.append(df)
        merged = pd.concat(dfs, ignore_index=True)
        groups = {"group1": ["sample_1", "sample_2"],
                   "group2": ["sample_3", "sample_4"]}
        dmr_df, _ = call_dmr(merged, groups, group1="group1", group2="group2",
                              q_threshold=0.05, delta_threshold=0.15)
        for ed in EMBEDDED_DMRS:
            matches = dmr_df[
                (dmr_df["chr"] == ed[0]) & (dmr_df["context"] == ed[3]) &
                (dmr_df["start"] < ed[2]) & (dmr_df["end"] > ed[1])
            ]
            if len(matches) > 0:
                expected_dir = ed[5]
                assert all(d == expected_dir for d in matches["direction"]), \
                    f"DMR {ed[0]}:{ed[1]}-{ed[2]} direction mismatch"

    def test_dmr_q_values_below_threshold(self, all_windows):
        """All detected DMRs have q-value below threshold."""
        dfs = []
        for sample, df in all_windows.items():
            df = df.copy()
            df["sample"] = sample
            dfs.append(df)
        merged = pd.concat(dfs, ignore_index=True)
        groups = {"group1": ["sample_1", "sample_2"],
                   "group2": ["sample_3", "sample_4"]}
        dmr_df, _ = call_dmr(merged, groups, group1="group1", group2="group2",
                              q_threshold=0.05, delta_threshold=0.15)
        if len(dmr_df) > 0:
            assert all(dmr_df["q_value"] < 0.05), "Some DMRs have q >= 0.05"

    def test_stats_genome_weighted(self, all_sites):
        """Genome row exists and n_sites equals sum of per-chr n_sites."""
        df = compute_stats(all_sites["sample_1"])
        genome = df[df["chr"] == "genome"]
        assert len(genome) > 0
        per_chr = df[df["chr"] != "genome"]
        for ctx in genome["context"].unique():
            g_sites = genome[genome["context"] == ctx]["n_sites"].sum()
            p_sites = per_chr[per_chr["context"] == ctx]["n_sites"].sum()
            assert g_sites == pytest.approx(p_sites), \
                f"Genome n_sites mismatch for {ctx}"
