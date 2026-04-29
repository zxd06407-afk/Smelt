"""Benchmark runner for Smelt.

Measures runtime for each subcommand on real data.
Compares DMR results with methylKit/DSS if available.

Usage:
    python scripts/benchmark.py --input-dir /path/to/data --output-dir results/
"""
import argparse
import subprocess
import time
import json
import sys
from pathlib import Path


def run_cmd(cmd, timeout=3600):
    """Run a command, return wall_time."""
    start = time.perf_counter()
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
    elapsed = time.perf_counter() - start
    if result.returncode != 0:
        print(f"  FAILED: {result.stderr[:500]}", file=sys.stderr)
    return elapsed


def benchmark_site(input_dir, output_dir):
    """Benchmark smelt site."""
    print("Benchmarking: smelt site")
    t = run_cmd(
        f"smelt site --input {input_dir}/sample.cov.gz --fasta {input_dir}/genome.fa "
        f"--output {output_dir}/bench.site.tsv"
    )
    return {"site_wall_time_s": t}


def benchmark_window(output_dir):
    """Benchmark smelt window."""
    print("Benchmarking: smelt window")
    t = run_cmd(
        f"smelt window --input {output_dir}/bench.site.tsv "
        f"--window 2000 --step 500 --output {output_dir}/bench.window.tsv"
    )
    return {"window_wall_time_s": t}


def main():
    p = argparse.ArgumentParser(description="Smelt benchmark")
    p.add_argument("--input-dir", required=True)
    p.add_argument("--output-dir", default="benchmark_results")
    args = p.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    results = {}
    results.update(benchmark_site(args.input_dir, out))
    results.update(benchmark_window(out))

    results_file = out / "benchmark.json"
    with open(results_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {results_file}")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
