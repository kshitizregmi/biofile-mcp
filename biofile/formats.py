"""Identify which genomics format a file is in."""
import gzip
from pathlib import Path

EXTENSIONS = {
    "VCF": (".vcf", ".vcf.gz", ".bcf"),
    "BAM": (".bam",),
    "FASTQ": (".fastq", ".fq", ".fastq.gz", ".fq.gz"),
    "h5ad": (".h5ad",),
}
INDEX_SUFFIXES = (".bai", ".csi", ".tbi", ".crai")

# First bytes of each format (after gzip/BGZF decompression where relevant)
_HDF5 = b"\x89HDF\r\n\x1a\n"
_GZIP = b"\x1f\x8b"


def is_index(path: str) -> bool:
    """True for index files (.bai, .csi, .tbi, .crai), which are never summarized."""
    return path.lower().endswith(INDEX_SUFFIXES)


def detect_format(path: str) -> str | None:
    """Return the format of a file: ``"VCF"``, ``"BAM"``, ``"FASTQ"`` or ``"h5ad"``.

    The file extension is checked first. If it is not recognized, the first
    bytes of the file are read instead, so a FASTQ saved as ``reads.txt`` or a
    VCF without its ``.vcf`` extension is still identified.

    Parameters
    ----------
    path : str
        Path to the file.

    Returns
    -------
    str or None
        The format name, or ``None`` for index files and unrecognized files.
    """
    name = Path(path).name.lower()
    if is_index(name):
        return None
    for fmt, suffixes in EXTENSIONS.items():
        if name.endswith(suffixes):
            return fmt
    return _sniff(path)


def _sniff(path: str) -> str | None:
    """Guess the format from the file's first bytes."""
    try:
        with open(path, "rb") as f:
            head = f.read(16)
        if head.startswith(_HDF5):
            return "h5ad"
        # BAM, BCF and compressed VCF/FASTQ are all gzip-compressed (BGZF)
        if head.startswith(_GZIP):
            with gzip.open(path, "rb") as f:
                head = f.read(16)
    except (OSError, EOFError):
        return None

    if head.startswith(b"BAM\x01"):
        return "BAM"
    if head.startswith((b"BCF", b"##fileformat=VCF")):
        return "VCF"
    if head.startswith(b"@"):
        return "FASTQ"
    return None
