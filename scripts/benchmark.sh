#!/usr/bin/env bash
# Smelt benchmark script — run on real WGBS data before release.
#
# Usage:
#   bash scripts/benchmark.sh /path/to/bismark/outputs/
#
# The input directory should contain:
#   - *.cov.gz files for each sample
#   - genome.fa reference
#
# Results are written to benchmark_results/

set -euo pipefail

INPUT_DIR="${1:?Usage: $0 <input-dir>}"
OUTPUT_DIR="benchmark_results"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "=== Smelt Benchmark ==="
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo

python "$SCRIPT_DIR/benchmark.py" --input-dir "$INPUT_DIR" --output-dir "$OUTPUT_DIR"

echo
echo "=== Benchmark Complete ==="
cat "$OUTPUT_DIR/benchmark.json"
