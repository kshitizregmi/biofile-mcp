"""Tests for :func:`biofile.detect_format`."""
import shutil

import pytest

from biofile import detect_format


@pytest.mark.parametrize(
    "name, expected",
    [
        ("test.vcf", "VCF"),
        ("test.vcf.gz", "VCF"),
        ("test.bam", "BAM"),
        ("sample_R1.fastq", "FASTQ"),
        ("sample_R1.fastq.gz", "FASTQ"),
        ("small.h5ad", "h5ad"),
        ("test.bam.bai", None),
        ("test.vcf.gz.tbi", None),
    ],
)
def test_detects_by_extension(data_dir, name, expected):
    assert detect_format(str(data_dir / name)) == expected


@pytest.mark.parametrize(
    "source, expected",
    [
        ("test.vcf", "VCF"),
        ("test.vcf.gz", "VCF"),
        ("test.bam", "BAM"),
        ("sample_R1.fastq", "FASTQ"),
        ("sample_R1.fastq.gz", "FASTQ"),
        ("small.h5ad", "h5ad"),
    ],
)
def test_detects_by_content_when_extension_is_unknown(data_dir, tmp_path, source, expected):
    renamed = tmp_path / "renamed.dat"
    shutil.copy(data_dir / source, renamed)
    assert detect_format(str(renamed)) == expected


def test_unrecognized_and_empty_files(tmp_path):
    text = tmp_path / "notes.txt"
    text.write_text("just some notes\n")
    empty = tmp_path / "empty.dat"
    empty.write_bytes(b"")
    assert detect_format(str(text)) is None
    assert detect_format(str(empty)) is None
    assert detect_format(str(tmp_path / "missing.dat")) is None
