"""Sequencing read (FASTQ) summaries."""
import os
from collections import Counter

import numpy as np
import pyfastx

# FASTQ quality characters are Phred scores shifted by 33 ("!" = Q0, "I" = Q40)
_PHRED_OFFSET = 33


def fastq_stats(path: str, max_reads: int | None = None) -> dict:
    """Summarize reads, base composition and base quality in a FASTQ file.

    Reads are streamed one at a time, so memory use stays constant regardless
    of file size. No index file is written next to the input.

    Parameters
    ----------
    path : str
        Path to a ``.fastq``/``.fq`` file, optionally gzip-compressed.
    max_reads : int, optional
        Stop after this many reads. Use for a quick estimate on large files;
        the result is then marked ``sampled``. By default every read is read.

    Returns
    -------
    dict
        ``file``
            The path that was read.
        ``reads``, ``sampled``
            Number of reads examined, and whether reading stopped early at
            ``max_reads``.
        ``total_bases``, ``min_len``, ``avg_len``, ``max_len``
            Read length statistics.
        ``mean_base_quality``
            Arithmetic mean of the Phred scores of all bases. This is higher
            than seqkit's ``AvgQual``, which averages error probabilities.
        ``q20_pct``, ``q30_pct``
            Percentage of bases with Phred quality of at least 20 and 30.
        ``gc_pct``
            Percentage of G and C bases (either case).
        ``n_bases``
            Number of ambiguous ``N``/``n`` bases.
        ``lowercase_bases``
            Number of lowercase bases. Raw sequencer output is uppercase, so a
            non-zero value usually means the reads were soft-masked or edited.
        ``other_characters``
            Characters other than A, C, G, T and N, such as IUPAC codes.
        ``base_counts``
            Case-sensitive count of every character seen.

        An empty file returns ``reads=0`` and a ``warning`` instead.

    Raises
    ------
    OSError
        If the file does not exist or cannot be read.
    RuntimeError
        If the file is not valid FASTQ.
    """
    # pyfastx refuses to open zero-byte files, so report them directly
    if os.path.getsize(path) == 0:
        return {"file": path, "reads": 0, "warning": "file is empty"}

    reads = 0
    total_bases = 0
    min_len = None
    max_len = 0
    base_counts = Counter()
    quality_sum = 0
    q20 = 0
    q30 = 0

    for _name, seq, qual in pyfastx.Fastx(path):
        reads += 1
        length = len(seq)
        total_bases += length
        min_len = length if min_len is None else min(min_len, length)
        max_len = max(max_len, length)
        base_counts.update(seq)

        # Convert the whole quality string to Phred scores in one step
        scores = np.frombuffer(qual.encode(), dtype=np.uint8) - _PHRED_OFFSET
        quality_sum += int(scores.sum())
        q20 += int((scores >= 20).sum())
        q30 += int((scores >= 30).sum())

        if max_reads and reads >= max_reads:
            break

    if reads == 0:
        return {"file": path, "reads": 0, "warning": "file has no reads"}

    gc = sum(base_counts[b] for b in "GCgc")
    lowercase = sum(n for b, n in base_counts.items() if b.islower())
    other = sum(n for b, n in base_counts.items() if b.upper() not in "ACGTN")

    return {
        "file": path,
        "reads": reads,
        "sampled": bool(max_reads and reads >= max_reads),
        "total_bases": total_bases,
        "min_len": min_len,
        "avg_len": round(total_bases / reads, 1),
        "max_len": max_len,
        "mean_base_quality": round(quality_sum / total_bases, 2),
        "q20_pct": round(100 * q20 / total_bases, 2),
        "q30_pct": round(100 * q30 / total_bases, 2),
        "gc_pct": round(100 * gc / total_bases, 2),
        "n_bases": base_counts["N"] + base_counts["n"],
        "lowercase_bases": lowercase,
        "other_characters": other,
        "base_counts": dict(base_counts),
    }


# Highest Phred score that can be written in Phred+33 ("~" = Q93)
_MAX_PHRED = 93

# Base composition columns, and a lookup from byte value to column (other = N)
_BASES = "ACGTN"
_BASE_INDEX = np.full(256, 4, dtype=np.int64)
for _i, _b in enumerate("ACGT"):
    _BASE_INDEX[ord(_b)] = _i
    _BASE_INDEX[ord(_b.lower())] = _i


def _bin_edges(n: int, max_bins: int) -> np.ndarray:
    """Start index of each bin when splitting n items into equal-width bins.

    Bins are one item wide when n <= max_bins; otherwise the smallest width
    that keeps the number of bins at or below max_bins.
    """
    width = -(-n // max_bins)  # ceiling division
    return np.arange(0, n, width)


def _quantile(cumulative: np.ndarray, totals: np.ndarray, q: float) -> np.ndarray:
    """Per-row quantile of a score histogram, given its row-wise cumulative sum."""
    return np.argmax(cumulative >= (q * totals)[:, None], axis=1)


def fastq_quality_profile(
    path: str, max_reads: int | None = 100_000, max_points: int = 300
) -> dict:
    """Base quality along the read, and the distribution of read lengths.

    This is the data behind a FastQC-style "per base sequence quality" plot.
    Use it to decide whether reads need 3' trimming and to what length.

    Parameters
    ----------
    path : str
        Path to a ``.fastq``/``.fq`` file, optionally gzip-compressed.
    max_reads : int or None, default 100000
        Number of reads to examine, taken from the start of the file. A sample
        of this size gives a stable profile; pass ``None`` to read every read.
    max_points : int, default 300
        Maximum number of positions returned. Reads up to this length are
        reported base by base; for longer reads, neighbouring positions are
        pooled into bins of equal width.

    Returns
    -------
    dict
        ``file``, ``reads``, ``sampled``
            The path, number of reads examined and whether that was a sample.
        ``position_start``, ``position_end``
            1-based first and last read position of each bin (equal when
            bins are single positions).
        ``mean``, ``median``, ``q1``, ``q3``, ``p10``, ``p90``
            Phred quality at each bin: mean, median, quartiles, and 10th/90th
            percentiles across all bases in that bin.
        ``reads_covering``
            Number of bases contributing to each bin. This falls toward the
            3' end when reads have different lengths.
        ``median_below_q30_from``, ``median_below_q20_from``
            First position at which the median falls below Q30 / Q20, or
            ``None`` if it never does. A natural trimming point.
        ``length_start``, ``length_end``, ``length_reads``
            Read-length histogram, at most 50 bins.
        ``base_pct``
            Percentage of each base at each bin, as lists keyed ``"A"``,
            ``"C"``, ``"G"``, ``"T"`` and ``"N"`` (any other character counts
            as N). Lines that diverge at the start of reads usually mean
            primers or adapters; for amplicon data such as 16S this is expected.
        ``read_quality``
            Reads per mean Phred quality (rounded), as ``quality`` and
            ``reads`` lists: FastQC's "per sequence quality scores".
        ``read_gc``
            Reads per GC percentage (rounded), as ``gc_pct`` and ``reads``
            lists. A second peak can point to contamination.

        An empty file returns ``reads=0`` and a ``warning`` instead.

    Raises
    ------
    OSError
        If the file does not exist or cannot be read.
    RuntimeError
        If the file is not valid FASTQ.
    """
    if os.path.getsize(path) == 0:
        return {"file": path, "reads": 0, "warning": "file is empty"}

    # hist[position, score] = number of bases with that score at that position
    hist = np.zeros((0, _MAX_PHRED + 1), dtype=np.int64)
    # bases[position, A/C/G/T/N] = number of each base at that position
    bases = np.zeros((0, len(_BASES)), dtype=np.int64)
    read_quality = np.zeros(_MAX_PHRED + 1, dtype=np.int64)
    read_gc = np.zeros(101, dtype=np.int64)
    lengths = Counter()
    reads = 0

    for _name, seq, qual in pyfastx.Fastx(path):
        scores = np.frombuffer(qual.encode(), dtype=np.uint8) - _PHRED_OFFSET
        scores = np.minimum(scores, _MAX_PHRED)
        codes = _BASE_INDEX[np.frombuffer(seq.encode(), dtype=np.uint8)]
        length = scores.size
        if length > hist.shape[0]:
            extra = length - hist.shape[0]
            hist = np.vstack([hist, np.zeros((extra, _MAX_PHRED + 1), np.int64)])
            bases = np.vstack([bases, np.zeros((extra, len(_BASES)), np.int64)])
        positions = np.arange(length)
        hist[positions, scores] += 1  # each position appears once per read
        bases[positions, codes] += 1
        read_quality[int(round(scores.mean()))] += 1
        gc = np.count_nonzero((codes == 1) | (codes == 2))
        read_gc[int(round(100 * gc / length))] += 1
        lengths[length] += 1
        reads += 1
        if max_reads and reads >= max_reads:
            break

    if reads == 0:
        return {"file": path, "reads": 0, "warning": "file has no reads"}

    # Pool positions into bins, then summarise each bin's score histogram
    starts = _bin_edges(hist.shape[0], max_points)
    ends = np.append(starts[1:], hist.shape[0])
    binned = np.add.reduceat(hist, starts, axis=0)
    base_bins = np.add.reduceat(bases, starts, axis=0)
    base_pct = 100 * base_bins / np.maximum(base_bins.sum(axis=1, keepdims=True), 1)
    totals = binned.sum(axis=1)
    cumulative = binned.cumsum(axis=1)
    mean = (binned * np.arange(_MAX_PHRED + 1)).sum(axis=1) / totals
    median = _quantile(cumulative, totals, 0.5)

    def first_below(threshold: int):
        below = np.flatnonzero(median < threshold)
        return int(starts[below[0]] + 1) if below.size else None

    # Read-length histogram
    lo, hi = min(lengths), max(lengths)
    edges = np.arange(lo, hi + 1) if hi - lo < 50 else np.linspace(lo, hi + 1, 51).astype(int)
    edges = np.unique(edges)
    bounds = np.append(edges, hi + 1)
    length_reads = [
        sum(n for L, n in lengths.items() if bounds[i] <= L < bounds[i + 1])
        for i in range(len(edges))
    ]

    return {
        "file": path,
        "reads": reads,
        "sampled": bool(max_reads and reads >= max_reads),
        "position_start": (starts + 1).tolist(),
        "position_end": ends.tolist(),
        "mean": np.round(mean, 2).tolist(),
        "median": median.tolist(),
        "q1": _quantile(cumulative, totals, 0.25).tolist(),
        "q3": _quantile(cumulative, totals, 0.75).tolist(),
        "p10": _quantile(cumulative, totals, 0.10).tolist(),
        "p90": _quantile(cumulative, totals, 0.90).tolist(),
        "reads_covering": totals.tolist(),
        "median_below_q30_from": first_below(30),
        "median_below_q20_from": first_below(20),
        "length_start": edges.tolist(),
        "length_end": (bounds[1:] - 1).tolist(),
        "length_reads": length_reads,
        "base_pct": {b: np.round(base_pct[:, i], 2).tolist() for i, b in enumerate(_BASES)},
        "read_quality": {
            "quality": np.flatnonzero(read_quality).tolist(),
            "reads": read_quality[read_quality > 0].tolist(),
        },
        "read_gc": {
            "gc_pct": np.flatnonzero(read_gc).tolist(),
            "reads": read_gc[read_gc > 0].tolist(),
        },
    }
