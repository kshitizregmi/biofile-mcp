"""Shared fixtures: build the known-answer test files once per test run."""
import gzip
import shutil

import pytest
import pysam

from make_test_data import make_test_data


@pytest.fixture(scope="session")
def data_dir(tmp_path_factory):
    """Folder with test.vcf, test.vcf.gz, test.bam(+.bai), sample_R1.fastq(.gz), small.h5ad."""
    folder = make_test_data(tmp_path_factory.mktemp("data"))

    # compressed + indexed VCF, made from the plain one
    shutil.copy(folder / "test.vcf", folder / "copy.vcf")
    pysam.tabix_index(str(folder / "copy.vcf"), preset="vcf", force=True)
    (folder / "copy.vcf.gz").rename(folder / "test.vcf.gz")
    (folder / "copy.vcf.gz.tbi").rename(folder / "test.vcf.gz.tbi")

    # gzipped FASTQ
    with open(folder / "sample_R1.fastq", "rb") as src, gzip.open(folder / "sample_R1.fastq.gz", "wb") as dst:
        shutil.copyfileobj(src, dst)

    return folder
