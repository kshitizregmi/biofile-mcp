"""Variant Call Format (VCF/BCF) summaries."""
import pysam
from pysam import bcftools

from biofile._histogram import bin_counts

# Most samples reported individually in per_sample
_MAX_SAMPLES = 100


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
        ``substitutions``
            Count of each of the 12 single-base changes, keyed ``"A>C"`` to
            ``"T>G"`` as written in the file (not collapsed by strand).
        ``qual``
            Sites per QUAL value as ``start``, ``end`` and ``count`` lists
            (at most 50 bins).
        ``indel_lengths``
            Indel sites per length as ``length`` and ``sites`` lists;
            deletions have negative lengths.
        ``genotype_depth``
            Genotypes per FORMAT/DP value as ``depth`` (labels, e.g. ``"30"``)
            and ``genotypes`` lists. A ``">500"`` label collects every depth
            above bcftools' upper limit.
        ``per_sample``
            One entry per sample (at most 100): ``sample``, ``hom_ref``,
            ``het``, ``hom_alt``, ``missing``, ``ts``, ``tv``, ``indels``,
            ``singletons`` and ``avg_depth``.

    Raises
    ------
    pysam.utils.SamtoolsError
        If the file is missing or is not a valid VCF/BCF.
    """
    # "-s -" adds per-sample statistics; bcftools rejects it for sites-only files
    # Silence htslib's "could not retrieve index" note; only the header is read
    verbosity = pysam.set_verbosity(0)
    try:
        with pysam.VariantFile(path) as vcf:
            has_samples = len(vcf.header.samples) > 0
    finally:
        pysam.set_verbosity(verbosity)
    report = bcftools.stats(*(["-s", "-"] if has_samples else []), path)

    summary = {"file": path}
    substitutions, qual, indels, depth, samples = {}, [], [], [], []

    for line in report.splitlines():
        fields = line.split("\t")
        kind = fields[0]

        # "SN" lines hold the overall counts, e.g. "SN  0  number of SNPs:  3"
        if kind == "SN":
            key = fields[2].replace("number of ", "").rstrip(":")
            summary[key] = int(fields[3])

        # "TSTV" line: transitions, transversions, ratio
        elif kind == "TSTV":
            summary["ts"] = int(fields[2])
            summary["tv"] = int(fields[3])
            summary["ts_tv"] = float(fields[4])

        # "ST" lines: one per substitution type, e.g. "ST  0  A>G  1"
        elif kind == "ST":
            substitutions[fields[2]] = int(fields[3])

        # "QUAL" lines: quality, SNPs, transitions, transversions, indels.
        # Records without a QUAL value (".") are left out.
        elif kind == "QUAL" and fields[2] != ".":
            qual.append((float(fields[2]), int(fields[3]) + int(fields[6])))

        # "IDD" lines: indel length (deletions negative), number of sites
        elif kind == "IDD":
            indels.append((int(fields[2]), int(fields[3])))

        # "DP" lines: depth bin (">500" for the last), number of genotypes
        elif kind == "DP":
            depth.append((fields[2], int(fields[3])))

        # "PSC" lines: per-sample counts
        elif kind == "PSC" and len(samples) < _MAX_SAMPLES:
            samples.append({
                "sample": fields[2],
                "hom_ref": int(fields[3]),
                "het": int(fields[5]),
                "hom_alt": int(fields[4]),
                "missing": int(fields[13]) if len(fields) > 13 else 0,
                "ts": int(fields[6]),
                "tv": int(fields[7]),
                "indels": int(fields[8]),
                "singletons": int(fields[10]),
                "avg_depth": float(fields[9]),
            })

    summary["substitutions"] = substitutions
    summary["qual"] = bin_counts([q for q, _ in qual], [n for _, n in qual], 50)
    summary["indel_lengths"] = {
        "length": [L for L, _ in sorted(indels)],
        "sites": [n for _, n in sorted(indels)],
    }
    summary["genotype_depth"] = {
        "depth": [d for d, _ in depth],
        "genotypes": [n for _, n in depth],
    }
    summary["per_sample"] = samples
    return summary
