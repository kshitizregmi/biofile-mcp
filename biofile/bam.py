"""Sequence alignment (BAM) summaries."""
import json
import re

import numpy as np
import pysam
from pysam import samtools

from biofile._histogram import bin_counts

# How many contig names to show when a region names a contig that isn't there
_CONTIGS_TO_SHOW = 10


def _region_contig(region: str, contigs) -> str:
    """Return the contig part of a samtools region string.

    Follows samtools: if the whole string is a contig name in the header, it is
    taken as the contig. Otherwise a trailing ``:start-end`` is removed.
    """
    # Some contig names (e.g. HLA alleles) contain colons themselves,
    # so check the header before splitting
    if region in contigs:
        return region
    contig, sep, span = region.rpartition(":")
    if sep and re.fullmatch(r"[\d,]+(-[\d,]*)?", span):
        return contig
    return region


def _require_contig(contig: str, contigs) -> None:
    """Raise a ValueError naming the header's contigs if contig is not one of them."""
    if contig not in contigs:
        shown = list(contigs[:_CONTIGS_TO_SHOW])
        more = f" and {len(contigs) - len(shown)} more" if len(contigs) > len(shown) else ""
        raise ValueError(
            f"Contig '{contig}' is not in the BAM header. Contigs in this file: {shown}{more}"
        )


def bam_stats(path: str, region: str | None = None) -> dict:
    """Summarize read counts in a BAM file, and optionally coverage of one region.

    Read counts come from ``samtools flagstat``, which reads the whole file.
    Region coverage comes from ``samtools coverage`` and uses the BAM index,
    so it is fast even for large files.

    Parameters
    ----------
    path : str
        Path to a ``.bam`` file.
    region : str, optional
        Region in samtools format, e.g. ``"chr17:43044295-43125483"`` (1-based,
        inclusive). The contig name must match the BAM header exactly
        (``"chr1"`` and ``"1"`` are different). Requires a ``.bai`` or ``.csi``
        index next to the BAM.

    Returns
    -------
    dict
        ``file``
            The path that was read.
        ``total_reads``, ``mapped_reads``, ``mapped_pct``, ``duplicates``
            Counts for QC-passed reads.
        ``region``, ``reads_in_region``, ``pct_bases_covered``, ``mean_depth``
            Present only when ``region`` is given. Duplicate, secondary,
            QC-failed and unmapped reads are excluded from these values.

    Raises
    ------
    ValueError
        If the region's contig is not in the BAM header. The message lists
        contig names from the header, which helps spot ``chr1`` vs ``1``
        naming mismatches.
    pysam.utils.SamtoolsError
        If the file is missing or invalid, the index is missing, or the region
        cannot be parsed.
    """
    flagstat = json.loads(samtools.flagstat("-O", "json", path))
    counts = flagstat["QC-passed reads"]

    summary = {
        "file": path,
        "total_reads": counts["total"],
        "mapped_reads": counts["mapped"],
        "mapped_pct": counts["mapped %"],
        "duplicates": counts["duplicates"],
    }

    if region:
        with pysam.AlignmentFile(path) as bam:
            contigs = bam.references
        _require_contig(_region_contig(region, contigs), contigs)

        # samtools coverage prints a header row and one row of values
        table = samtools.coverage("-r", region, path).strip().splitlines()
        header = table[0].lstrip("#").split("\t")
        values = table[1].split("\t")
        coverage = dict(zip(header, values))

        summary["region"] = region
        summary["reads_in_region"] = int(coverage["numreads"])
        summary["pct_bases_covered"] = float(coverage["coverage"])
        summary["mean_depth"] = float(coverage["meandepth"])

    return summary


# Longest region bam_depth_profile will read; per-base depth is held in memory
_MAX_PROFILE_BP = 10_000_000


def _parse_region(region: str, contigs, lengths) -> tuple[str, int, int]:
    """Split a samtools region into (contig, start, end), 1-based inclusive."""
    contig = _region_contig(region, contigs)
    _require_contig(contig, contigs)
    contig_length = lengths[contigs.index(contig)]
    span = region[len(contig) + 1:].replace(",", "")
    if not span:
        return contig, 1, contig_length

    first, _, last = span.partition("-")
    start = int(first)
    end = int(last) if last else contig_length
    if not 1 <= start <= end:
        raise ValueError(f"Invalid region '{region}': start must be at least 1 and not after end")
    return contig, start, min(end, contig_length)


def bam_depth_profile(
    path: str, region: str, min_depth: int = 20, max_points: int = 500
) -> dict:
    """Read depth along a region, for plotting coverage and finding gaps.

    Depth is counted per base with the same read filters as
    ``samtools coverage``: unmapped, secondary, QC-failed and duplicate reads
    are excluded, and no base-quality cut-off is applied.

    Parameters
    ----------
    path : str
        Path to an indexed ``.bam`` file.
    region : str
        Region in samtools format, e.g. ``"chr17:43044295-43125483"`` (1-based,
        inclusive), or a contig name for the whole contig. At most 10 Mb.
    min_depth : int, default 20
        Depth regarded as adequate; used for ``pct_bases_at_min_depth``.
    max_points : int, default 500
        Maximum number of points returned. Longer regions are split into
        bins of equal width and the mean depth of each bin is reported.

    Returns
    -------
    dict
        ``file``, ``region``
            The path and the region actually read, as ``contig:start-end``.
        ``position_start``, ``position_end``, ``depth``
            1-based bin bounds and mean depth per bin.
        ``mean_depth``, ``median_depth``, ``min_depth``, ``max_depth``
            Per-base depth statistics over the whole region.
        ``min_depth_threshold``, ``pct_bases_at_min_depth``
            The threshold used and the percentage of bases at or above it.

    Raises
    ------
    ValueError
        If the contig is not in the BAM header, the coordinates are invalid,
        the region is longer than 10 Mb, or the BAM has no index.
    """
    with pysam.AlignmentFile(path) as bam:
        contig, start, end = _parse_region(region, bam.references, bam.lengths)
        size = end - start + 1
        if size > _MAX_PROFILE_BP:
            raise ValueError(
                f"Region is {size:,} bp; the limit is {_MAX_PROFILE_BP:,} bp. Choose a smaller region."
            )
        # count_coverage returns four arrays (A, C, G, T); their sum is the depth
        acgt = bam.count_coverage(contig, start - 1, end, quality_threshold=0, read_callback="all")
    depth = np.sum(np.array(acgt, dtype=np.int64), axis=0)

    width = -(-size // max_points)  # ceiling division
    starts = np.arange(0, size, width)
    binned = np.add.reduceat(depth, starts) / np.diff(np.append(starts, size))

    return {
        "file": path,
        "region": f"{contig}:{start}-{end}",
        "position_start": (start + starts).tolist(),
        "position_end": (start + np.append(starts[1:], size) - 1).tolist(),
        "depth": np.round(binned, 2).tolist(),
        "mean_depth": round(float(depth.mean()), 2),
        "median_depth": float(np.median(depth)),
        "min_depth": int(depth.min()),
        "max_depth": int(depth.max()),
        "min_depth_threshold": min_depth,
        "pct_bases_at_min_depth": round(100 * float((depth >= min_depth).mean()), 2),
    }


# Bins kept for each histogram returned by bam_qc
_QC_BINS = 100


def _stats_sections(report: str) -> dict:
    """Group ``samtools stats`` output lines by their section code (SN, IS, ...)."""
    sections = {}
    for line in report.splitlines():
        if line and not line.startswith("#"):
            fields = line.split("\t")
            sections.setdefault(fields[0], []).append(fields[1:])
    return sections


def bam_qc(path: str, region: str | None = None) -> dict:
    """Detailed alignment QC for a BAM file, from ``samtools stats`` and ``idxstats``.

    Covers what a quality report usually shows: pairing, error rate, base
    quality, mapping quality, insert size, read length, the distribution of
    coverage depth, and reads per contig. ``samtools stats`` reads every
    alignment, so for large files pass a ``region`` to limit the work.

    Parameters
    ----------
    path : str
        Path to a ``.bam`` file. A ``.bai`` or ``.csi`` index is needed for
        ``region`` and for the per-contig counts.
    region : str, optional
        Restrict the statistics to a region in samtools format, e.g.
        ``"chr22"`` or ``"chr17:43044295-43125483"``.

    Returns
    -------
    dict
        ``file``, ``region``
            The path and the region analysed (``None`` for the whole file).
        ``summary``
            Read and base counts: ``reads``, ``mapped``, ``mapped_pct``,
            ``paired``, ``properly_paired_pct``, ``duplicates``, ``mq0``
            (mapped with MAPQ 0), ``secondary``, ``supplementary``,
            ``bases_mapped``, ``mismatches``, ``error_rate`` (mismatches per
            mapped base), ``avg_length``, ``max_length``, ``avg_quality``,
            ``insert_size_mean``, ``insert_size_sd``, ``inward_pairs``,
            ``outward_pairs`` and ``pairs_on_different_contigs``.
        ``mapq``
            Reads per mapping-quality value, as ``mapq`` and ``reads`` lists.
        ``insert_size``
            Read pairs per insert size, as ``start``, ``end`` and ``count``
            lists (at most 100 bins). Empty for single-end data.
        ``read_length``
            Reads per read length, in the same form.
        ``coverage``
            Covered bases per depth, as ``depth`` and ``bases`` lists, and
            ``pct_bases_at_least`` giving the share of covered bases at or
            above 1, 5, 10, 20, 30, 50 and 100×. Only covered bases are
            counted, so these are percentages of covered bases, not of the
            whole reference. As in ``samtools stats``, duplicate reads count
            toward depth here (unlike ``bam_depth_profile``).
        ``contigs``
            Mapped and unmapped reads per contig from the index, largest
            first (at most 50), plus ``unplaced_unmapped`` for reads with no
            position. Empty if the BAM has no index.

    Raises
    ------
    ValueError
        If the region's contig is not in the BAM header.
    pysam.utils.SamtoolsError
        If the file is missing or invalid, or a region is given without an index.
    """
    args = [path]
    if region:
        with pysam.AlignmentFile(path) as bam:
            _require_contig(_region_contig(region, bam.references), bam.references)
        args.append(region)
    sections = _stats_sections(samtools.stats(*args))

    sn = {fields[0].rstrip(":"): fields[1] for fields in sections.get("SN", [])}

    def number(key, kind=float):
        return kind(float(sn.get(key, 0)))

    reads = number("sequences", int)
    mapped = number("reads mapped", int)
    paired = number("reads paired", int)
    summary = {
        "reads": reads,
        "mapped": mapped,
        "mapped_pct": round(100 * mapped / reads, 2) if reads else 0.0,
        "paired": paired,
        "properly_paired_pct": (
            round(100 * number("reads properly paired") / paired, 2) if paired else 0.0
        ),
        "duplicates": number("reads duplicated", int),
        "mq0": number("reads MQ0", int),
        "secondary": number("non-primary alignments", int),
        "supplementary": number("supplementary alignments", int),
        "bases_mapped": number("bases mapped (cigar)", int),
        "mismatches": number("mismatches", int),
        "error_rate": float(f"{number('error rate'):.6g}"),
        "avg_length": number("average length"),
        "max_length": number("maximum length", int),
        "avg_quality": number("average quality"),
        "insert_size_mean": number("insert size average"),
        "insert_size_sd": number("insert size standard deviation"),
        "inward_pairs": number("inward oriented pairs", int),
        "outward_pairs": number("outward oriented pairs", int),
        "pairs_on_different_contigs": number("pairs on different chromosomes", int),
    }

    mapq = sorted((int(q), int(n)) for q, n in sections.get("MAPQ", []))
    insert = [(int(f[0]), int(f[1])) for f in sections.get("IS", []) if int(f[1])]
    lengths = [(int(f[0]), int(f[1])) for f in sections.get("RL", [])]

    # COV lines look like "[10-10]  10  2345"; the last one is "[1000<]  1000  n"
    depth_bases = sorted((int(f[1]), int(f[2])) for f in sections.get("COV", []))
    covered = sum(b for _, b in depth_bases)
    at_least = {
        threshold: round(100 * sum(b for d, b in depth_bases if d >= threshold) / covered, 2)
        if covered else 0.0
        for threshold in (1, 5, 10, 20, 30, 50, 100)
    }
    shown = depth_bases if len(depth_bases) <= 500 else depth_bases[:500]

    contigs, unplaced = [], 0
    try:
        for line in samtools.idxstats(path).splitlines():
            name, length, n_mapped, n_unmapped = line.split("\t")
            if name == "*":
                unplaced = int(n_unmapped)
            else:
                contigs.append({"contig": name, "length": int(length),
                                "mapped": int(n_mapped), "unmapped": int(n_unmapped)})
    except pysam.utils.SamtoolsError:
        pass  # no index: per-contig counts are unavailable
    contigs.sort(key=lambda c: c["mapped"], reverse=True)

    return {
        "file": path,
        "region": region,
        "summary": summary,
        "mapq": {"mapq": [q for q, _ in mapq], "reads": [n for _, n in mapq]},
        "insert_size": bin_counts([s for s, _ in insert], [n for _, n in insert], _QC_BINS),
        "read_length": bin_counts([s for s, _ in lengths], [n for _, n in lengths], _QC_BINS),
        "coverage": {
            "depth": [d for d, _ in shown],
            "bases": [b for _, b in shown],
            "pct_bases_at_least": {f"{t}x": v for t, v in at_least.items()},
        },
        "contigs": contigs[:50],
        "unplaced_unmapped": unplaced,
    }
