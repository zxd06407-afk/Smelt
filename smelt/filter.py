"""Quality filters for methylation data."""


def filter_by_depth(df, min_depth=5):
    """Filter sites/windows by minimum read depth (total coverage).

    Returns a new filtered DataFrame.
    """
    return df[df["total"] >= min_depth].copy()


def filter_by_sites(df, min_sites=10):
    """Filter windows by minimum number of covered cytosine sites.

    Use n_sites column if present.
    """
    if "n_sites" not in df.columns:
        raise ValueError("DataFrame has no 'n_sites' column")
    return df[df["n_sites"] >= min_sites].copy()
