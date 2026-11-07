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
