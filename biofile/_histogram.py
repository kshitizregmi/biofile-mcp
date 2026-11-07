"""Helpers for turning count tables into compact histograms."""
import numpy as np


def bin_counts(values, counts, max_bins: int = 100) -> dict:
    """Pool a value -> count table into at most ``max_bins`` equal-width bins.

    ``values`` are numbers (integers or floats) and ``counts`` how often each
    occurs. Returns ``start``, ``end`` and ``count`` lists. When the values are
    integers spanning no more than ``max_bins``, every value keeps its own bin.
    """
    values = np.asarray(values, dtype=float)
    counts = np.asarray(counts, dtype=np.int64)
    if values.size == 0:
        return {"start": [], "end": [], "count": []}

    lo, hi = float(values.min()), float(values.max())
    integers = bool(np.all(values == np.round(values)))
    if integers:
        # Whole-number values: bins of a whole-number width, one per value when they fit
        width = float(-(-(hi - lo + 1) // max_bins))  # ceiling division
        n_bins = int((hi - lo) // width) + 1
        starts = lo + width * np.arange(n_bins)
        ends = np.minimum(starts + width - 1, hi)
    else:
        n_bins = max_bins if hi > lo else 1
        width = (hi - lo) / n_bins or 1.0
        starts = lo + width * np.arange(n_bins)
        ends = starts + width

    index = np.minimum(((values - lo) // width).astype(int), n_bins - 1)
    binned = np.bincount(index, weights=counts, minlength=n_bins).astype(np.int64)

    def tidy(a):
        return [int(v) if integers else round(float(v), 4) for v in a]

    return {"start": tidy(starts), "end": tidy(ends), "count": binned.tolist()}
