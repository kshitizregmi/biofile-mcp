"""Variant Call Format (VCF/BCF) summaries."""
from pysam import bcftools


def vcf_stats(path: str) -> dict:
    """Summarize a VCF or BCF file with ``bcftools stats``.

    Reads every record in the file, so run time grows with file size.

    Parameters
    ----------
    path : str
        Path to a ``.vcf``, ``.vcf.gz`` or ``.bcf`` file.

    Returns
    -------
    dict
        ``file``
            The path that was read.
        ``samples``, ``records``
            Number of samples and variant records.
        ``SNPs``, ``MNPs``, ``indels``, ``others``, ``no-ALTs``
            Record counts by variant type.
        ``multiallelic sites``, ``multiallelic SNP sites``
            Sites with more than one alternate allele.
        ``ts``, ``tv``, ``ts_tv``
            Transition and transversion counts and their ratio. Roughly 2.0-2.1
            is expected for human whole-genome data and 2.8-3.3 for exomes.

    Raises
    ------
    pysam.utils.SamtoolsError
        If the file is missing or is not a valid VCF/BCF.
    """
    report = bcftools.stats(path)

    summary = {"file": path}

    for line in report.splitlines():
        fields = line.split("\t")

        # "SN" lines hold the overall counts, e.g. "SN  0  number of SNPs:  3"
        if fields[0] == "SN":
            key = fields[2].replace("number of ", "").rstrip(":")
            summary[key] = int(fields[3])

        # "TSTV" line: transitions, transversions, ratio
        if fields[0] == "TSTV":
            summary["ts"] = int(fields[2])
            summary["tv"] = int(fields[3])
            summary["ts_tv"] = float(fields[4])

    return summary
