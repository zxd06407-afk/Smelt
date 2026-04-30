"""Smelt CLI -- WGBS methylation downstream analysis."""

import sys
from pathlib import Path
from typing import List, Optional

import typer
import pandas as pd

from smelt.io import read_cov, read_gtf, read_bed, read_fasta
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
    """Write DataFrame to TSV, using default_name if no output specified."""
    path = output if output else default_name
    df.to_csv(path, sep="\t", index=False)
    print(f"Wrote {len(df)} rows to {path}", file=sys.stderr)


def _output_base(input_path: str, suffix: str) -> str:
    """Generate output filename from input basename."""
    base = Path(input_path).stem
    base = base.replace(".cov", "")
    return f"{base}.{suffix}"


@app.command()
def site(
    input_file: str = typer.Option(..., "--input", "-i", help="BISMARK cov.gz file"),
    fasta: Optional[str] = typer.Option(None, "--fasta", "-f", help="Reference genome FASTA"),
    context_file: Optional[str] = typer.Option(None, "--context-file", help="Pre-annotated context file"),
    min_depth: int = typer.Option(5, "--min-depth", help="Minimum read depth"),
    merge_cpg: bool = typer.Option(True, "--merge-cpg-strands/--no-merge-cpg-strands"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
):
    """Compute per-cytosine methylation with sequence context classification."""
    cov = read_cov(input_file)
    fa = read_fasta(fasta) if fasta else None
    ctx = pd.read_csv(context_file, sep="\t") if context_file else None
    result = compute_site_methylation(
        cov, fasta=fa, context_df=ctx,
        min_depth=min_depth, merge_cpg=merge_cpg, threads=threads,
    )
    out = output if output else _output_base(input_file, "site.tsv")
    _write_output(result, out, out)


@app.command()
def window(
    input_file: str = typer.Option(..., "--input", "-i", help="Site TSV file or BISMARK cov.gz"),
    fasta: Optional[str] = typer.Option(None, "--fasta", "-f", help="Reference genome FASTA (if input is cov.gz)"),
    window_size: int = typer.Option(2000, "--window", "-w", help="Window size in bp"),
    step: int = typer.Option(500, "--step", "-s", help="Step size in bp"),
    min_sites: int = typer.Option(10, "--min-sites", help="Minimum sites per window"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
):
    """Aggregate methylation into sliding windows."""
    df = pd.read_csv(input_file, sep="\t")
    if "context" not in df.columns:
        cov = read_cov(input_file)
        fa = read_fasta(fasta) if fasta else None
        if fa is None:
            raise typer.BadParameter("--fasta required when input is a raw cov.gz file")
        df = compute_site_methylation(cov, fasta=fa, threads=threads)
    result = compute_windows(df, window_size=window_size, step=step, min_sites=min_sites, threads=threads)
    suffix = f"window{window_size}s{step}.tsv"
    out = output if output else _output_base(input_file, suffix)
    _write_output(result, out, out)


@app.command()
def element(
    input_file: str = typer.Option(..., "--input", "-i", help="Site TSV file"),
    gtf: str = typer.Option(..., "--gtf", "-g", help="GTF annotation file"),
    features: Optional[str] = typer.Option(None, "--features", help="Feature types (comma-separated, e.g. exon,CDS)"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
):
    """Compute methylation levels for GTF feature types."""
    sites = pd.read_csv(input_file, sep="\t")
    gtf_df = read_gtf(gtf)
    feat_list = features.split(",") if features else None
    result = compute_elements(sites, gtf_df, features=feat_list)
    out = output if output else _output_base(input_file, "element.tsv")
    _write_output(result, out, out)


@app.command()
def metaplot(
    input_file: str = typer.Option(..., "--input", "-i", help="Site TSV file"),
    gtf: Optional[str] = typer.Option(None, "--gtf", "-g", help="GTF file for gene intervals"),
    bed: Optional[str] = typer.Option(None, "--bed", "-b", help="BED file for gene intervals"),
    upstream: int = typer.Option(2000, "--upstream", help="Upstream bp"),
    downstream: int = typer.Option(2000, "--downstream", help="Downstream bp"),
    body_bins: int = typer.Option(50, "--body-bins", help="Gene body equal-ratio bins"),
    up_bins: int = typer.Option(20, "--up-bins", help="Upstream bins"),
    down_bins: int = typer.Option(20, "--down-bins", help="Downstream bins"),
    context: str = typer.Option("CpG", "--context", "-c", help="Sequence context"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
):
    """Compute metaplot methylation signal across gene bodies and flanking regions."""
    sites = pd.read_csv(input_file, sep="\t")
    if gtf:
        intervals = read_gtf(gtf)
        intervals = intervals[intervals["feature"] == "gene"]
    elif bed:
        intervals = read_bed(bed)
    else:
        raise typer.BadParameter("Either --gtf or --bed must be provided")
    result = compute_metaplot(
        sites, intervals, context=context,
        upstream=upstream, downstream=downstream,
        body_bins=body_bins, up_bins=up_bins, down_bins=down_bins,
    )
    out = output if output else _output_base(input_file, "metaplot.tsv")
    _write_output(result, out, out)


@app.command()
def custom(
    input_file: str = typer.Option(..., "--input", "-i", help="Site TSV file"),
    bed: str = typer.Option(..., "--bed", "-b", help="BED file with custom intervals"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
):
    """Compute methylation for user-defined BED intervals."""
    sites = pd.read_csv(input_file, sep="\t")
    bed_df = read_bed(bed)
    result = compute_custom(sites, bed_df)
    out = output if output else _output_base(input_file, "custom.tsv")
    _write_output(result, out, out)


@app.command()
def dmr(
    samples: List[str] = typer.Option(..., "--samples", help="Sample definitions: name=file (repeatable)"),
    group1: str = typer.Option(..., "--group1", help="Comma-separated sample names for group 1"),
    group2: str = typer.Option(..., "--group2", help="Comma-separated sample names for group 2"),
    q_threshold: float = typer.Option(0.05, "--q-threshold", help="FDR q-value cutoff"),
    delta_threshold: float = typer.Option(0.2, "--delta-threshold", help="Minimum methylation difference"),
    min_samples: int = typer.Option(2, "--min-samples-per-group", help="Minimum samples per group"),
    fdr_all: bool = typer.Option(False, "--fdr-all-contexts", help="Pool all contexts for BH correction"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
):
    """Call differentially methylated regions between two groups."""
    sample_map = {}
    for s in samples:
        name, path = s.split("=", 1)
        sample_map[name] = path

    dfs = []
    for name, path in sample_map.items():
        df = pd.read_csv(path, sep="\t")
        df["sample"] = name
        dfs.append(df)
    merged = pd.concat(dfs, ignore_index=True)

    g1 = [s.strip() for s in group1.split(",")]
    g2 = [s.strip() for s in group2.split(",")]
    groups = {"group1": g1, "group2": g2}

    dmr_df, excl_df = call_dmr(
        merged, groups, group1="group1", group2="group2",
        q_threshold=q_threshold, delta_threshold=delta_threshold,
        min_samples_per_group=min_samples, fdr_by_context=not fdr_all,
    )

    base = _output_base(sample_map[g1[0]], "")
    dmr_out = output if output else f"{base}dmr.tsv"
    excl_out = dmr_out.replace(".tsv", "_excluded.tsv")
    _write_output(dmr_df, dmr_out, dmr_out)
    _write_output(excl_df, excl_out, excl_out)


@app.command()
def stats(
    input_file: str = typer.Option(..., "--input", "-i", help="Site or window TSV file"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
):
    """Compute genome and chromosome-level methylation statistics."""
    df = pd.read_csv(input_file, sep="\t")
    result = compute_stats(df)
    out = output if output else _output_base(input_file, "stats.tsv")
    _write_output(result, out, out)


if __name__ == "__main__":
    app()
