"""FASTQ tools."""
import os
from collections import Counter

import numpy as np
import pyfastx


def fastq_stats(path: str, max_reads: int | None = None) -> dict:
    """Summarize a FASTQ or FASTQ.gz file: number of reads, total bases, read length
    (min/avg/max), mean base quality, % bases >= Q20 and Q30, GC %, N count,
    lowercase count and per-base counts. If max_reads is set, only the first
    max_reads reads are checked (faster for big files)."""
    if os.path.getsize(path) == 0:  # pyfastx can't open empty files
        return {"file": path, "reads": 0, "warning": "file is empty"}

    reads = 0
    total_bases = 0
    min_len = None
    max_len = 0
    base_counts = Counter()
    qual_sum = 0
    q20 = 0
    q30 = 0

    for name, seq, qual in pyfastx.Fastx(path):
        reads += 1
        length = len(seq)
        total_bases += length
        min_len = length if min_len is None else min(min_len, length)
        max_len = max(max_len, length)
        base_counts.update(seq)

        q = np.frombuffer(qual.encode(), dtype=np.uint8) - 33
        qual_sum += int(q.sum())
        q20 += int((q >= 20).sum())
        q30 += int((q >= 30).sum())

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
        "mean_base_quality": round(qual_sum / total_bases, 2),
        "q20_pct": round(100 * q20 / total_bases, 2),
        "q30_pct": round(100 * q30 / total_bases, 2),
        "gc_pct": round(100 * gc / total_bases, 2),
        "n_bases": base_counts["N"] + base_counts["n"],
        "lowercase_bases": lowercase,
        "other_characters": other,
        "base_counts": dict(base_counts),
    }
