"""VCF tools."""
from pysam import bcftools


def vcf_stats(path: str) -> dict:
    """Summarize a VCF or VCF.gz file using bcftools stats: number of samples,
    records, SNPs, MNPs, indels, multiallelic sites, transitions (ts),
    transversions (tv) and the Ts/Tv ratio."""
    text = bcftools.stats(path)

    summary = {"file": path}

    for line in text.splitlines():
        parts = line.split("\t")

        if parts[0] == "SN":
            key = parts[2].replace("number of ", "").rstrip(":")
            summary[key] = int(parts[3])

        if parts[0] == "TSTV":
            summary["ts"] = int(parts[2])
            summary["tv"] = int(parts[3])
            summary["ts_tv"] = float(parts[4])

    return summary
