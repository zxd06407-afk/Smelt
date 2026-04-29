# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Smelt — WGBS methylation downstream analysis CLI tool. Reads BISMARK output, computes methylation at sites / sliding windows / genomic elements, DMR calling, and summary statistics. Outputs TSV for external plotting.

- **Language:** Python 3
- **CLI:** Typer
- **Dependency management:** uv
- **Logging:** structlog
- **Testing:** pytest

## Design spec

Full design document at `docs/superpowers/specs/2026-04-29-smelt-design.md`. Read it before any implementation work.

## Commands

```bash
# Build / install
uv sync                          # Install dependencies
uv pip install -e .              # Install smelt in editable mode

# Tests
pytest                           # All tests
pytest tests/test_site.py        # Single module
pytest -k "test_cpg_merge"       # Single test by name

# Linting (once configured)
ruff check smelt/ tests/
mypy smelt/

# Run smelt during development
python -m smelt.cli --help
python -m smelt.cli site --input data/sample.cov.gz
```

## Architecture

7 subcommands under `smelt/` as one module per subcommand:

| Module | Subcommand | Purpose |
|--------|-----------|---------|
| `site.py` | `smelt site` | Per-cytosine methylation + context classification |
| `window.py` | `smelt window` | Fixed-window sliding-window aggregation |
| `element.py` | `smelt element` | GTF feature-type methylation |
| `metaplot.py` | `smelt metaplot` | Gene body ± flanking equal-ratio binning |
| `custom.py` | `smelt custom` | BED custom intervals |
| `dmr.py` | `smelt dmr` | Window-level Fisher + BH, per-context FDR |
| `stats.py` | `smelt stats` | Genome/chromosome summary by context |

Shared: `io.py` (file I/O), `filter.py` (coverage filtering), `utils.py` (coordinate conversion, context helpers), `cli.py` (Typer entry point).

## Coordinate conventions

**Internal: 0-based half-open `[start, end)`.** Always convert at I/O boundaries:

| Input format | Convention | Conversion |
|-------------|-----------|------------|
| BISMARK cov.gz | 1-based pos | `pos - 1` |
| GTF/GFF | 1-based closed | `start-1`, `end` unchanged |
| BED | 0-based half-open | None |

Output is always 0-based half-open. Coordinate helpers live in `utils.py`.

## Key design rules

- **Sequence contexts:** CpG, CHH, CHG. Context determined from FASTA reference or pre-annotated BISMARK context file.
- **CpG strand merging:** Default on (`--merge-cpg-strands`). CHH/CHG always strand-separated.
- **Chromosome naming:** Strict exact match. Mismatch → error, not silent data loss.
- **FDR correction:** Per-context (CpG/CHH/CHG separately). `--fdr-all-contexts` to pool.
- **Per-chromosome processing:** Default single-threaded, `--threads N` for multi-process.
- **DMR double threshold:** Both FDR < 0.05 and |delta| > 0.2 required (configurable).
- **Output naming:** Auto-generated from input basename + parameters, `--output` overrides.

## Testing tiers

1. **Unit tests** — constructed minimal data, CI, validate individual functions
2. **Integration tests** — simulated data (2 fake chr, 4 samples, 8 embedded DMRs), CI < 2 min
3. **Benchmark** — real public data, cross-tool comparison, manual pre-release
