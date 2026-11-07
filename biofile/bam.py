"""BAM tools."""
import json

from pysam import samtools


def bam_stats(path: str, region: str | None = None) -> dict:
    """Summarize a BAM file: total, mapped and duplicate read counts.
    If region is given (e.g. 'chr1:1001-1100'), also return read count,
    % of bases covered and mean depth for that region. Needs a .bai index for regions."""
    data = json.loads(samtools.flagstat("-O", "json", path))
    stats = data["QC-passed reads"]

    summary = {
        "file": path,
        "total_reads": stats["total"],
        "mapped_reads": stats["mapped"],
        "mapped_pct": stats["mapped %"],
        "duplicates": stats["duplicates"],
    }

    if region:
        text = samtools.coverage("-r", region, path)
        lines = text.strip().splitlines()
        header = lines[0].lstrip("#").split("\t")
        values = lines[1].split("\t")
        cov = dict(zip(header, values))

        summary["region"] = region
        summary["reads_in_region"] = int(cov["numreads"])
        summary["pct_bases_covered"] = float(cov["coverage"])
        summary["mean_depth"] = float(cov["meandepth"])

    return summary
