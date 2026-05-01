"""Smelt CLI -- WGBS methylation downstream analysis."""

import sys
from pathlib import Path
from typing import List, Optional

import typer
import pandas as pd

from smelt.io import (read_cx_report, read_gtf, read_bed, read_bismark_sam)
from smelt.site import compute_site_methylation
from smelt.window import compute_windows
from smelt.element import compute_elements
from smelt.metaplot import compute_metaplot
from smelt.custom import compute_custom
from smelt.dmr import call_dmr
from smelt.stats import compute_stats

app = typer.Typer(
    name="smelt",
    help="WGBS methylation downstream analysis tool",
)


def _write_output(df: pd.DataFrame, output: Optional[str], default_name: str):
    """Write DataFrame, using default_name if no output specified.
    Defaults to Parquet for intermediate files; use .parquet suffix for TSV."""
    path = output if output else default_name
    if "chr" in df.columns:
        df = df.copy()
        df["chr"] = df["chr"].astype(str)
    ext = Path(path).suffix.lower()
    if ext == ".tsv":
        df.to_csv(path, sep="\t", index=False)
    else:
        if not ext:
            path = path + ".parquet"
        df.to_parquet(path, index=False)
    print(f"Wrote {len(df)} rows to {path}", file=sys.stderr)


def _read_input(path: str) -> pd.DataFrame:
    """Read a DataFrame from TSV or Parquet based on extension."""
    ext = Path(path).suffix.lower()
    if ext == ".parquet":
        return pd.read_parquet(path)
    else:
        return pd.read_csv(path, sep="\t")


def _output_base(input_path: str, suffix: str) -> str:
    """Generate output filename from input basename."""
    base = Path(input_path).stem
    return f"{base}.{suffix}"


def _infer_sample_name(input_path: str) -> str:
    """Infer sample name from input filename."""
    return Path(input_path).stem


def _validate_chromosomes(site_df, annot_df, annot_label):
    """Check that chromosome names overlap between site and annotation files."""
    site_chroms = set(site_df["chr"].unique())
    annot_chroms = set(annot_df["chr"].unique())
    common = site_chroms & annot_chroms
    if not common:
        raise typer.BadParameter(
            f"No chromosomes in common between site file and {annot_label}.\n"
            f"  Site file chromosomes: {sorted(site_chroms)}\n"
            f"  {annot_label} chromosomes: {sorted(annot_chroms)}"
        )


@app.command()
def site(
    input_file: str = typer.Option(..., "--input", "-i",
        help="BISMARK SAM/BAM or CX_report.txt file"),
    min_depth: int = typer.Option(5, "--min-depth", help="Minimum read depth"),
    merge_cpg: bool = typer.Option(True, "--merge-cpg-strands/--no-merge-cpg-strands"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name",
        help="Sample identifier (default: inferred from filename)"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
):
    """Compute per-cytosine methylation with sequence context classification."""
    input_lower = input_file.lower()
    if input_lower.endswith((".bam", ".sam")):
        df = read_bismark_sam(input_file)
    elif "cx_report" in input_lower:
        df = read_cx_report(input_file)
    else:
        raise typer.BadParameter(
            "Unsupported input format. Use BISMARK SAM/BAM or CX_report.txt."
        )
    name = sample_name if sample_name else _infer_sample_name(input_file)
    result = compute_site_methylation(
        df, min_depth=min_depth, merge_cpg=merge_cpg,
        sample_name=name, threads=threads,
    )
    out = output if output else f"{name}.site.parquet"
    _write_output(result, out, out)


@app.command()
def window(
    input_file: str = typer.Option(..., "--input", "-i",
        help="Site TSV file or BISMARK SAM/BAM/CX_report.txt"),
    window_size: int = typer.Option(2000, "--window", "-w", help="Window size in bp"),
    step: int = typer.Option(500, "--step", "-s", help="Step size in bp"),
    min_sites: int = typer.Option(10, "--min-sites", help="Minimum sites per window"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name",
        help="Sample identifier (default: inferred from filename)"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
):
    """Aggregate methylation into sliding windows."""
    input_lower = input_file.lower()
    if input_lower.endswith((".bam", ".sam")):
        df = read_bismark_sam(input_file)
    elif "cx_report" in input_lower:
        df = read_cx_report(input_file)
    else:
        df = _read_input(input_file)

    if "context" not in df.columns:
        name = sample_name if sample_name else _infer_sample_name(input_file)
        df = compute_site_methylation(df, min_depth=5, sample_name=name, threads=threads)
    elif sample_name and "sample" not in df.columns:
        df["sample"] = sample_name

    result = compute_windows(df, window_size=window_size, step=step,
                              min_sites=min_sites, threads=threads)
    suffix = f"window{window_size}s{step}.parquet"
    name = sample_name if sample_name else _infer_sample_name(input_file)
    out = output if output else f"{name}.{suffix}"
    _write_output(result, out, out)


@app.command()
def element(
    input_file: str = typer.Option(..., "--input", "-i", help="Site TSV file"),
    gtf: str = typer.Option(..., "--gtf", "-g", help="GTF annotation file"),
    features: Optional[str] = typer.Option(None, "--features",
        help="Feature types (comma-separated, e.g. exon,CDS)"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name",
        help="Sample identifier (default: inferred from site filename)"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
):
    """Compute methylation levels for GTF feature types."""
    sites = _read_input(input_file)
    gtf_df = read_gtf(gtf)
    _validate_chromosomes(sites, gtf_df, "GTF file")
    feat_list = features.split(",") if features else None
    name = sample_name if sample_name else _infer_sample_name(input_file)
    result = compute_elements(sites, gtf_df, features=feat_list, sample_name=name)
    out = output if output else f"{name}.element.parquet"
    _write_output(result, out, out)


@app.command()
def metaplot(
    input_file: str = typer.Option(..., "--input", "-i", help="Site TSV file"),
    gtf: Optional[str] = typer.Option(None, "--gtf", "-g",
        help="GTF file for gene intervals"),
    bed: Optional[str] = typer.Option(None, "--bed", "-b",
        help="BED file for gene intervals"),
    upstream: int = typer.Option(2000, "--upstream", help="Upstream bp"),
    downstream: int = typer.Option(2000, "--downstream", help="Downstream bp"),
    body_bins: int = typer.Option(50, "--body-bins", help="Gene body equal-ratio bins"),
    up_bins: int = typer.Option(20, "--up-bins", help="Upstream bins"),
    down_bins: int = typer.Option(20, "--down-bins", help="Downstream bins"),
    context: str = typer.Option("CpG", "--context", "-c", help="Sequence context"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name",
        help="Sample identifier (default: inferred from site filename)"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
):
    """Compute metaplot methylation signal across gene bodies and flanking regions."""
    sites = _read_input(input_file)
    if gtf:
        intervals = read_gtf(gtf)
        intervals = intervals[intervals["feature"] == "gene"]
    elif bed:
        intervals = read_bed(bed)
    else:
        raise typer.BadParameter("Either --gtf or --bed must be provided")
    _validate_chromosomes(sites, intervals, "GTF/BED file")
    name = sample_name if sample_name else _infer_sample_name(input_file)
    result = compute_metaplot(
        sites, intervals, context=context,
        upstream=upstream, downstream=downstream,
        body_bins=body_bins, up_bins=up_bins, down_bins=down_bins,
        sample_name=name,
    )
    out = output if output else f"{name}.metaplot.parquet"
    _write_output(result, out, out)


@app.command()
def custom(
    input_file: str = typer.Option(..., "--input", "-i", help="Site TSV file"),
    bed: str = typer.Option(..., "--bed", "-b",
        help="BED file with custom intervals"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name",
        help="Sample identifier (default: inferred from site filename)"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
):
    """Compute methylation for user-defined BED intervals."""
    sites = _read_input(input_file)
    bed_df = read_bed(bed)
    _validate_chromosomes(sites, bed_df, "BED file")
    name = sample_name if sample_name else _infer_sample_name(input_file)
    result = compute_custom(sites, bed_df, sample_name=name)
    out = output if output else f"{name}.custom.parquet"
    _write_output(result, out, out)


@app.command()
def dmr(
    samples: Optional[List[str]] = typer.Option(None, "--samples",
        help="Sample definitions: name=file (repeatable). "
             "May be omitted if --group1/--group2 contain name=file entries."),
    group1: str = typer.Option(..., "--group1",
        help="Sample names (comma-separated), or name=file entries"),
    group2: str = typer.Option(..., "--group2",
        help="Sample names (comma-separated), or name=file entries"),
    q_threshold: float = typer.Option(0.05, "--q-threshold",
        help="FDR q-value cutoff"),
    delta_threshold: float = typer.Option(0.2, "--delta-threshold",
        help="Minimum methylation difference"),
    min_samples: int = typer.Option(2, "--min-samples-per-group",
        help="Minimum samples per group"),
    fdr_all: bool = typer.Option(False, "--fdr-all-contexts",
        help="Pool all contexts for BH correction"),
    output: Optional[str] = typer.Option(None, "--output", "-o",
        help="Output file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
):
    """Call differentially methylated regions between two groups.

    Two syntaxes are supported:

    \b
    1. Explicit --samples:
       smelt dmr --samples a=file1 --samples b=file2 \\
                 --group1 a --group2 b

    \b
    2. Inline name=file in group args:
       smelt dmr --group1 a=file1,b=file2 --group2 c=file3,d=file4
    """
    sample_map = {}
    if samples:
        for s in samples:
            name, path = s.split("=", 1)
            sample_map[name] = path
    else:
        for entry in group1.split(","):
            entry = entry.strip()
            if "=" in entry:
                name, path = entry.split("=", 1)
                sample_map[name.strip()] = path.strip()
        for entry in group2.split(","):
            entry = entry.strip()
            if "=" in entry:
                name, path = entry.split("=", 1)
                sample_map[name.strip()] = path.strip()
    if not sample_map:
        raise typer.BadParameter(
            "No sample files specified. Use --samples name=file or "
            "provide name=file entries in --group1/--group2."
        )

    dfs = []
    for name, path in sample_map.items():
        df = _read_input(path)
        df["sample"] = name
        dfs.append(df)
    merged = pd.concat(dfs, ignore_index=True)

    # Extract sample names from group args (strip file paths if present)
    g1 = []
    for s in group1.split(","):
        s = s.strip()
        g1.append(s.split("=")[0].strip() if "=" in s else s)
    g2 = []
    for s in group2.split(","):
        s = s.strip()
        g2.append(s.split("=")[0].strip() if "=" in s else s)
    groups = {"group1": g1, "group2": g2}

    dmr_df, excl_df = call_dmr(
        merged, groups, group1="group1", group2="group2",
        q_threshold=q_threshold, delta_threshold=delta_threshold,
        min_samples_per_group=min_samples, fdr_by_context=not fdr_all,
        threads=threads,
    )

    base = Path(sample_map[g1[0]]).stem
    dmr_out = output if output else f"{base}.dmr.parquet"
    excl_out = dmr_out.replace(".parquet", "_excluded.parquet")
    _write_output(dmr_df, dmr_out, dmr_out)
    _write_output(excl_df, excl_out, excl_out)


@app.command()
def stats(
    input_file: str = typer.Option(..., "--input", "-i",
        help="Site or window TSV file"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name",
        help="Sample identifier (default: inferred from filename)"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
):
    """Compute genome and chromosome-level methylation statistics."""
    df = _read_input(input_file)
    name = sample_name if sample_name else _infer_sample_name(input_file)
    result = compute_stats(df, sample_name=name)
    out = output if output else f"{name}.stats.parquet"
    _write_output(result, out, out)


if __name__ == "__main__":
    app()
