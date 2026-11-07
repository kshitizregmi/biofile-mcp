"""Sequence alignment (BAM) summaries."""
import json

from pysam import samtools


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
    pysam.utils.SamtoolsError
        If the file is missing or invalid, the index is missing, or the region
        cannot be parsed or does not exist in the BAM header.
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
