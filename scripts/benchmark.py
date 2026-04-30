#!/usr/bin/env python3
"""Smelt vs BatMeth2 benchmark runner.

Compares methylation analysis results on Arabidopsis met1 vs WT WGBS data.
Outputs a JSON report with per-module metrics.
"""

import argparse
import csv
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr, spearmanr


def run_cmd(cmd, timeout=7200):
    """Run a shell command, return (exit_code, wall_time)."""
    start = time.perf_counter()
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True,
                            timeout=timeout)
    elapsed = time.perf_counter() - start
    if result.returncode != 0:
        print(f"  WARNING: command failed (exit={result.returncode})")
        print(f"  stderr: {result.stderr[:300]}", file=sys.stderr)
    return result.returncode, elapsed


def load_tsv(path):
    """Load a TSV file into list of dicts."""
    p = Path(path)
    if not p.exists():
        return None
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def compare_ratios(r1, r2):
    """Compute Pearson/Spearman correlation between two ratio arrays."""
    a = np.array(r1, dtype=float)
    b = np.array(r2, dtype=float)
    mask = ~(np.isnan(a) | np.isnan(b))
    a, b = a[mask], b[mask]
    if len(a) < 3:
        return {"pearson_r": None, "spearman_rho": None, "n": len(a)}
    r, p_r = pearsonr(a, b)
    rho, p_rho = spearmanr(a, b)
    return {"pearson_r": round(float(r), 6),
            "pearson_p": round(float(p_r), 8),
            "spearman_rho": round(float(rho), 6),
            "spearman_p": round(float(p_rho), 8),
            "n": len(a)}


def compare_dmr(smelt_data, batmeth_data):
    """Compare DMR overlap between two tools."""
    if not smelt_data or not batmeth_data:
        return {"error": "missing DMR data"}
    smelt_set = set()
    for r in smelt_data:
        smelt_set.add((r["chr"], int(r["start"]), int(r["end"]),
                        r.get("context", "")))
    batmeth_set = set()
    for r in batmeth_data:
        batmeth_set.add((r.get("chrom", r.get("chr", "")),
                          int(r.get("start", 0)),
                          int(r.get("end", 0)),
                          r.get("context", "")))
    inter = smelt_set & batmeth_set
    union = smelt_set | batmeth_set
    jaccard = len(inter) / len(union) if union else 0
    return {"smelt_count": len(smelt_set),
            "batmeth_count": len(batmeth_set),
            "overlap_count": len(inter),
            "jaccard_index": round(jaccard, 4)}


def run_site_benchmark(bams, out_dir, threads):
    """Compare site-level methylation."""
    print("\n=== Site-Level Benchmark ===")
    result = {}
    for name, bam in bams.items():
        out = out_dir / f"{name}.site.tsv"
        cmd = (f"smelt site --input {bam} --sample-name {name} "
               f"--min-depth 5 -o {out} -t {threads}")
        rc, t = run_cmd(cmd)
        data = load_tsv(out)
        result[name] = {"wall_time_s": round(t, 1),
                        "n_sites": len(data) if data else 0}
        print(f"  {name}: {result[name]['n_sites']} sites, {t:.0f}s")

    # Cross-tool comparison: site-level methylation ratio correlation
    # (BatMeth2 calmeth comparison done after BatMeth2 runs)
    return result


def run_window_benchmark(site_tsvs, out_dir, threads):
    """Generate unified windows and run custom interval methylation."""
    print("\n=== Window Benchmark ===")
    result = {}

    # Generate unified window BED from first sample's window output
    first_sample = list(site_tsvs.keys())[0]
    window_tsv = out_dir / f"{first_sample}.window.tsv"
    window_bed = out_dir / "unified_windows.bed"

    cmd = (f"smelt window --input {site_tsvs[first_sample]} "
           f"-w 2000 -s 1000 --min-sites 5 -o {window_tsv} -t {threads}")
    run_cmd(cmd)

    # Extract unique windows as BED
    data = load_tsv(window_tsv)
    if data:
        seen = set()
        with open(window_bed, "w") as f:
            for row in data:
                key = (row["chr"], row["start"], row["end"], row["context"])
                if key not in seen:
                    seen.add(key)
                    f.write(f"{row['chr']}\t{row['start']}\t{row['end']}\t"
                            f"{row['context']}_{row['start']}_{row['end']}\n")
        result["n_windows"] = len(seen)
        print(f"  Unified windows: {len(seen)}")

    # Smelt custom on unified windows per sample
    for name, site_tsv in site_tsvs.items():
        out = out_dir / f"{name}.custom.tsv"
        cmd = (f"smelt custom --input {site_tsv} --bed {window_bed} "
               f"--sample-name {name} -o {out}")
        rc, t = run_cmd(cmd)
        data = load_tsv(out)
        result[name] = {"wall_time_s": round(t, 1),
                        "n_windows": len(data) if data else 0}
        if data:
            result[name]["mean_ratio"] = round(
                np.mean([float(r["ratio"]) for r in data]), 4)

    return result


def run_dmr_benchmark(group1_windows, group2_windows, out_dir, threads):
    """Run Smelt DMR on unified windows."""
    print("\n=== DMR Benchmark ===")
    result = {}

    # Smelt DMR
    samples_args = " ".join(f"{k}={v}" for k, v in
                            {**group1_windows, **group2_windows}.items())
    g1_names = ",".join(group1_windows.keys())
    g2_names = ",".join(group2_windows.keys())

    out_dmr = out_dir / "smelt.dmr.tsv"
    cmd = (f"smelt dmr --samples {samples_args} "
           f"--group1 {g1_names} --group2 {g2_names} "
           f"--q-threshold 0.05 --delta-threshold 0.2 "
           f"-o {out_dmr} -t {threads}")
    rc, t = run_cmd(cmd)
    data = load_tsv(out_dmr)

    result["smelt"] = {"wall_time_s": round(t, 1),
                       "n_dmr": len(data) if data else 0}
    if data:
        for ctx in ["CpG", "CHH", "CHG"]:
            ctx_data = [r for r in data if r.get("context") == ctx]
            if ctx_data:
                result["smelt"][f"n_{ctx}"] = len(ctx_data)
                result["smelt"][f"mean_delta_{ctx}"] = round(
                    np.mean([abs(float(r["delta"])) for r in ctx_data]), 4)
        print(f"  Smelt DMRs: {result['smelt']['n_dmr']} total "
              f"(CpG={result['smelt'].get('n_CpG',0)}, "
              f"CHH={result['smelt'].get('n_CHH',0)}, "
              f"CHG={result['smelt'].get('n_CHG',0)})")

    return result


def run_stats_benchmark(site_tsvs, out_dir, threads):
    """Compare genome/chromosome statistics."""
    print("\n=== Stats Benchmark ===")
    result = {}
    for name, site_tsv in site_tsvs.items():
        out = out_dir / f"{name}.stats.tsv"
        cmd = (f"smelt stats --input {site_tsv} "
               f"--sample-name {name} -o {out}")
        rc, t = run_cmd(cmd)
        data = load_tsv(out)
        if data:
            genome_rows = [r for r in data if r["chr"] == "genome"]
            result[name] = {"wall_time_s": round(t, 1)}
            for gr in genome_rows:
                ctx = gr["context"]
                result[name][f"genome_{ctx}_ratio"] = round(
                    float(gr["mean_ratio"]), 4)
                result[name][f"genome_{ctx}_sites"] = int(gr["n_sites"])
            print(f"  {name}: {result[name]}")

    return result


def main():
    p = argparse.ArgumentParser(description="Smelt Benchmark Runner")
    p.add_argument("--bam-dir", help="Directory with BISMARK BAM files")
    p.add_argument("--fastq-dir", help="Directory with FASTQ files "
                   "(for running BISMARK)")
    p.add_argument("--genome", help="Reference genome FASTA for BISMARK")
    p.add_argument("--output-dir", default="benchmark_results")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--module",
                   choices=["site", "window", "dmr", "stats", "all"],
                   default="all")
    p.add_argument("--skip-bismark", action="store_true",
                   help="Skip BISMARK alignment (BAM files already exist)")
    args = p.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Sample definitions
    samples = {
        "WT_rep1": "SRR921028",
        "WT_rep2": "SRR921029",
        "met1_rep1": "SRR921034",
        "met1_rep2": "SRR921035",
    }
    group1_windows = {"WT_rep1": None, "WT_rep2": None}
    group2_windows = {"met1_rep1": None, "met1_rep2": None}

    # Step 1: BISMARK alignment (if needed)
    bam_dir = args.bam_dir
    if not args.skip_bismark and args.fastq_dir:
        print("=== BISMARK Alignment ===")
        fastq_dir = Path(args.fastq_dir)
        bam_dir = out_dir / "bam"
        bam_dir.mkdir(exist_ok=True)
        for name, srr in samples.items():
            r1 = fastq_dir / f"{srr}_1.fastq.gz"
            r2 = fastq_dir / f"{srr}_2.fastq.gz"
            out_bam = bam_dir / f"{name}.bam"
            if out_bam.exists():
                print(f"  {name}: BAM already exists, skipping")
                continue
            cmd = (f"bismark --bowtie2 -p {args.threads} "
                   f"--genome {Path(args.genome).parent} "
                   f"-1 {r1} -2 {r2} -o {bam_dir}/")
            rc, t = run_cmd(cmd)
            print(f"  {name}: aligned in {t:.0f}s")
    elif not bam_dir:
        print("ERROR: --bam-dir or --fastq-dir required")
        sys.exit(1)

    # Build BAM file map
    bams = {}
    for name in samples:
        bam_path = Path(bam_dir) / f"{name}.bam"
        if bam_path.exists():
            bams[name] = str(bam_path)
        else:
            bam_path = Path(bam_dir) / f"{name}_pe.bam"
            if bam_path.exists():
                bams[name] = str(bam_path)

    if len(bams) < 4:
        print(f"WARNING: Only {len(bams)}/4 BAM files found")
        print(f"  Found: {list(bams.keys())}")
        if len(bams) < 2:
            print("Need at least 2 BAM files for comparison, exiting")
            sys.exit(1)

    all_results = {}

    # Step 2: Site-level benchmark
    if args.module in ("site", "all"):
        all_results["site"] = run_site_benchmark(bams, out_dir, args.threads)

    # Step 3: Window benchmark
    if args.module in ("window", "all"):
        site_tsvs = {name: str(out_dir / f"{name}.site.tsv") for name in bams}
        # Verify site TSVs exist
        for name, path in list(site_tsvs.items()):
            if not Path(path).exists():
                print(f"  Site TSV for {name} not found, removing from window analysis")
                del site_tsvs[name]
        if site_tsvs:
            all_results["window"] = run_window_benchmark(
                site_tsvs, out_dir, args.threads)

            # Build window file paths for DMR
            for name in site_tsvs:
                custom_tsv = out_dir / f"{name}.custom.tsv"
                if custom_tsv.exists():
                    if name in group1_windows:
                        group1_windows[name] = str(custom_tsv)
                    else:
                        group2_windows[name] = str(custom_tsv)

    # Step 4: DMR benchmark
    if args.module in ("dmr", "all"):
        g1 = {k: v for k, v in group1_windows.items() if v}
        g2 = {k: v for k, v in group2_windows.items() if v}
        if len(g1) >= 1 and len(g2) >= 1:
            all_results["dmr"] = run_dmr_benchmark(g1, g2, out_dir, args.threads)
        else:
            print("  Skipping DMR: need at least 1 sample per group")

    # Step 5: Stats benchmark
    if args.module in ("stats", "all"):
        site_tsvs = {name: str(out_dir / f"{name}.site.tsv") for name in bams}
        existing = {k: v for k, v in site_tsvs.items() if Path(v).exists()}
        if existing:
            all_results["stats"] = run_stats_benchmark(
                existing, out_dir, args.threads)

    # Save report
    report_path = out_dir / "benchmark_report.json"
    with open(report_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\n=== Report saved to {report_path} ===")
    print(json.dumps(all_results, indent=2))


if __name__ == "__main__":
    main()
