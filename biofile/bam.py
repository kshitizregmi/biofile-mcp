"""Sequence alignment (BAM) summaries."""
import json
import re

import pysam
from pysam import samtools

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
        contig = _region_contig(region, contigs)
        if contig not in contigs:
            shown = list(contigs[:_CONTIGS_TO_SHOW])
            more = f" and {len(contigs) - len(shown)} more" if len(contigs) > len(shown) else ""
            raise ValueError(
                f"Contig '{contig}' is not in the BAM header. Contigs in this file: {shown}{more}"
            )

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
