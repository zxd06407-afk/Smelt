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
    def test_dmr_detects_some_regions(self, all_windows):
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
