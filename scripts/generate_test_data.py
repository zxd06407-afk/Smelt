"""Generate simulated WGBS data for integration testing.

Produces:
  - 2 fake chromosomes (chrA ~100kb, chrB ~50kb) as FASTA
  - 4 BISMARK cov files (2 per group) with embedded ground truth
  - 1 GTF with 20 genes
  - 1 BED with 5 custom intervals

Ground truth: 8 DMRs embedded (5 CpG, 2 CHH, 1 CHG).
"""

import random
import gzip
import argparse
from pathlib import Path

random.seed(42)

CHROM_SIZES = {"chrA": 100_000, "chrB": 50_000}
SAMPLE_GROUPS = {
    "group1": ["sample_1", "sample_2"],
    "group2": ["sample_3", "sample_4"],
}

# Pre-defined DMRs: (chr, start, end, context, delta, direction)
EMBEDDED_DMRS = [
    ("chrA", 10000, 14000, "CpG", 0.8, "hyper"),
    ("chrA", 50000, 52000, "CpG", 0.6, "hypo"),
    ("chrA", 70000, 72000, "CpG", 0.4, "hyper"),
    ("chrA", 85000, 87000, "CpG", 0.25, "hypo"),
    ("chrA", 90000, 92000, "CpG", 0.15, "hyper"),
    ("chrB", 10000, 14000, "CHH", 0.7, "hyper"),
    ("chrB", 25000, 29000, "CHH", 0.5, "hypo"),
    ("chrB", 40000, 42000, "CHG", 0.55, "hyper"),
]

BASES = ["A", "C", "G", "T"]


def _is_in_dmr(chrom, pos, context):
    """Check if a position falls within any embedded DMR."""
    for dmr_chr, dmr_start, dmr_end, dmr_ctx, _delta, _dir in EMBEDDED_DMRS:
        if chrom == dmr_chr and context == dmr_ctx and dmr_start <= pos <= dmr_end:
            return _delta, _dir
    return None


def generate_fasta(output_dir):
    """Generate a single multi-chromosome FASTA."""
    path = output_dir / "genome.fa"
    with open(path, "w") as f:
        for chrom, size in CHROM_SIZES.items():
            f.write(f">{chrom}\n")
            seq = []
            for i in range(size):
                # Inject CpG, CHH, CHG patterns
                if i % 30 < 10:
                    seq.append("C")
                elif i % 30 < 15:
                    seq.append("G")
                else:
                    seq.append(random.choice(BASES))
            for i in range(0, size, 60):
                f.write("".join(seq[i:i+60]) + "\n")
    print(f"Wrote {path} ({sum(CHROM_SIZES.values())} bp)")


def _determine_context(chrom, pos):
    """Return a context based on position hash for realistic mix."""
    h = hash((chrom, pos)) % 100
    if h < 60:
        return "CpG"
    elif h < 85:
        return "CHH"
    else:
        return "CHG"


def generate_cov_files(output_dir):
    """Generate BISMARK coverage files with embedded DMRs.

    All samples share the same cytosine positions (fixed step per chromosome).
    Variation comes from depth, methylation ratio, and group-specific DMR effects.
    """
    # Build shared position list with pre-determined contexts
    shared_sites = []
    for chrom, size in CHROM_SIZES.items():
        step = 100  # fixed step so all samples have sites at identical positions
        for pos in range(1, size + 1, step):
            ctx = _determine_context(chrom, pos)
            shared_sites.append((chrom, pos, ctx))
    print(f"Shared site template: {len(shared_sites)} positions across all chromosomes")

    for group, samples in SAMPLE_GROUPS.items():
        for sample in samples:
            rows = []
            for chrom, pos, ctx in shared_sites:
                base_ratio = 0.8 if ctx == "CpG" else (0.1 if ctx == "CHH" else 0.05)
                # Add biological variation: random perturbation around base_ratio
                base_ratio += random.gauss(0, 0.05)
                base_ratio = max(0, min(1, base_ratio))

                depth = random.randint(10, 40)
                meth = int(depth * base_ratio)

                # Apply DMR effect for group2 samples
                dmr = _is_in_dmr(chrom, pos, ctx)
                if dmr and group == "group2":
                    delta, direction = dmr
                    if direction == "hyper":
                        meth += int(depth * delta * 0.8)
                    else:
                        meth -= int(depth * delta * 0.8)
                meth = max(0, min(depth, meth))
                unmet = depth - meth
                pct = meth / depth * 100
                rows.append(f"{chrom}\t{pos}\t{pos}\t{pct:.1f}\t{unmet}\t{meth}")

            path = output_dir / f"{sample}.cov"
            with open(path, "w") as f:
                f.write("\n".join(rows))
            gz_path = output_dir / f"{sample}.cov.gz"
            with gzip.open(gz_path, "wt") as f:
                f.write("\n".join(rows))
            print(f"Wrote {path} ({len(rows)} sites) and .gz")


def generate_gtf(output_dir):
    """Generate a GTF with 20 genes per chromosome."""
    genes = []
    gene_lengths = [500, 800, 1200, 2000, 3500, 5000, 8000,
                    600, 1500, 2800, 4500, 7000, 10000, 15000,
                    25000, 40000, 550, 1800, 3300, 6200]
    path = output_dir / "test_genes.gtf"
    with open(path, "w") as f:
        gene_idx = 0
        for chrom in CHROM_SIZES:
            max_end = CHROM_SIZES[chrom] - 10000
            for gl in [gl for gl in gene_lengths if gl < max_end - 5000]:
                pos = random.randint(5000, max_end - gl)
                gene_idx += 1
                gene_id = f"gene{gene_idx}"
                f.write(f"{chrom}\ttest\tgene\t{pos}\t{pos+gl}\t.\t+\t.\tgene_id \"{gene_id}\";\n")
                n_exons = random.randint(2, 5)
                exon_size = max(50, gl // n_exons)
                for j in range(n_exons):
                    e_start = pos + j * exon_size
                    e_end = min(e_start + exon_size - 10, pos + gl)
                    f.write(f"{chrom}\ttest\texon\t{e_start}\t{e_end}\t.\t+\t.\tgene_id \"{gene_id}\"; transcript_id \"{gene_id}_t1\";\n")
                if gl > 300:
                    cds_start = pos + min(50, gl // 4)
                    cds_end = pos + gl - min(50, gl // 4)
                    f.write(f"{chrom}\ttest\tCDS\t{cds_start}\t{cds_end}\t.\t+\t.\tgene_id \"{gene_id}\"; transcript_id \"{gene_id}_t1\";\n")
    print(f"Wrote {path} ({gene_idx} genes)")


def generate_bed(output_dir):
    """Generate a BED with 5 custom intervals."""
    path = output_dir / "test_regions.bed"
    intervals = [(1, "chrA", 15000, 25000), (2, "chrA", 40000, 48000),
                 (3, "chrA", 60000, 68000), (4, "chrB", 5000, 15000),
                 (5, "chrB", 30000, 38000)]
    with open(path, "w") as f:
        for idx, chrom, start, end in intervals:
            f.write(f"{chrom}\t{start}\t{end}\tregion_{idx}\n")
    print(f"Wrote {path}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", default="tests/data/integration")
    args = p.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    generate_fasta(out)
    generate_cov_files(out)
    generate_gtf(out)
    generate_bed(out)
    print("Done.")
