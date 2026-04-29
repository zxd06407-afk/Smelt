"""Site-level methylation with sequence context classification."""
import pandas as pd
from smelt.utils import classify_context, merge_cpg_strands
from smelt.filter import filter_by_depth


def compute_site_methylation(cov_df, fasta=None, context_df=None,
                              min_depth=5, merge_cpg=True):
    """Compute per-cytosine methylation with sequence context.

    Args:
        cov_df: DataFrame from read_cov() with 0-based positions
        fasta: pyfaidx Fasta object for context classification
        context_df: Optional DataFrame with pre-computed context per site
        min_depth: Minimum read depth filter
        merge_cpg: Merge CpG counts from both strands

    Returns:
        DataFrame: chr, pos, context, strand, meth, unmeth, total, ratio
    """
    df = cov_df.copy()

    # Filter by depth
    df = filter_by_depth(df, min_depth)

    # Determine strand: BISMARK OT (original top) or OB (original bottom)
    # Without explicit strand info from cov format, default to "+"
    df["strand"] = "+"

    # Assign context
    if context_df is not None:
        df = df.merge(context_df, on=["chr", "pos"], how="left")
        if "context" not in df.columns:
            raise ValueError("context_df must contain a 'context' column")
    elif fasta is not None:
        contexts = []
        for _, row in df.iterrows():
            chrom = row["chr"]
            pos = row["pos"]  # 0-based
            try:
                seq = fasta[chrom][pos:pos+3].seq
                ctx = classify_context(seq)
            except (KeyError, IndexError):
                ctx = None
            contexts.append(ctx)
        df["context"] = contexts
    else:
        raise ValueError("Either fasta or context_df must be provided")

    # Drop sites where context could not be determined
    df = df.dropna(subset=["context"])

    # Merge CpG strands
    if merge_cpg:
        df = merge_cpg_strands(df)

    return df.reset_index(drop=True)
