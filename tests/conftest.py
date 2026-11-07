"""Shared pytest fixtures."""
import gzip
import shutil

import pytest
import pysam

from make_test_data import make_test_data


@pytest.fixture(scope="session")
def data_dir(tmp_path_factory):
    """Build the known-answer test files once per test session.

    Adds a bgzipped, tabix-indexed ``test.vcf.gz`` and a gzipped
    ``sample_R1.fastq.gz`` to the files written by ``make_test_data``.
    """
    folder = make_test_data(tmp_path_factory.mktemp("data"))

    # tabix_index compresses the file in place, so work on a copy
    shutil.copy(folder / "test.vcf", folder / "copy.vcf")
    pysam.tabix_index(str(folder / "copy.vcf"), preset="vcf", force=True)
    (folder / "copy.vcf.gz").rename(folder / "test.vcf.gz")
    (folder / "copy.vcf.gz.tbi").rename(folder / "test.vcf.gz.tbi")

    with open(folder / "sample_R1.fastq", "rb") as src, gzip.open(folder / "sample_R1.fastq.gz", "wb") as dst:
        shutil.copyfileobj(src, dst)

    return folder
