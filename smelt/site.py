"""Site-level methylation with sequence context classification."""
import pandas as pd
from smelt.utils import classify_context, merge_cpg_strands, parallel_chromosomes
from smelt.filter import filter_by_depth


def _classify_chromosome(chrom_df, fasta=None, context_df=None):
    """Assign context to sites on a single chromosome."""
    df = chrom_df.copy()

    if context_df is not None:
        chrom_name = df["chr"].iloc[0]
        ctx_chrom = context_df[context_df["chr"] == chrom_name]
        if not ctx_chrom.empty:
            df = df.merge(ctx_chrom, on=["chr", "pos"], how="left")
        if "context" not in df.columns:
            df["context"] = None
    elif fasta is not None:
        contexts = []
        chrom_name = df["chr"].iloc[0]
        for _, row in df.iterrows():
            pos = row["pos"]
            try:
                seq = fasta[chrom_name][pos:pos+3].seq
                ctx = classify_context(seq)
            except (KeyError, IndexError):
                ctx = None
            contexts.append(ctx)
        df["context"] = contexts

    return df.dropna(subset=["context"])


def compute_site_methylation(cov_df, fasta=None, context_df=None,
                              min_depth=5, merge_cpg=True, threads=1):
    """Compute per-cytosine methylation with sequence context.

    Args:
        cov_df: DataFrame from read_cov() with 0-based positions
        fasta: pyfaidx Fasta object for context classification
        context_df: Optional DataFrame with pre-computed context per site
        min_depth: Minimum read depth filter
        merge_cpg: Merge CpG counts from both strands
        threads: Number of worker processes for per-chromosome parallelism

    Returns:
        DataFrame: chr, pos, context, strand, meth, unmeth, total, ratio
    """
    if fasta is None and context_df is None:
        raise ValueError("Either fasta or context_df must be provided")

    df = filter_by_depth(cov_df, min_depth)
    df["strand"] = "+"

    df = parallel_chromosomes(
        df, _classify_chromosome, threads=threads,
        fasta=fasta, context_df=context_df,
    )

    if merge_cpg:
        df = merge_cpg_strands(df)

    return df.reset_index(drop=True)
