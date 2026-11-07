"""Tests for :func:`biofile.bam_stats`."""
import pytest
from pysam.utils import SamtoolsError

from biofile import bam_stats
from biofile.bam import _region_contig


def test_read_counts(data_dir):
    r = bam_stats(str(data_dir / "test.bam"))
    assert r["total_reads"] == 25
    assert r["mapped_reads"] == 23
    assert r["mapped_pct"] == 92.0
    assert r["duplicates"] == 3
    assert "region" not in r


def test_region_coverage_skips_duplicates(data_dir):
    r = bam_stats(str(data_dir / "test.bam"), region="chr1:1001-1100")
    assert r["reads_in_region"] == 20  # duplicates are excluded
    assert r["pct_bases_covered"] == 50.0
    assert r["mean_depth"] == 10.0


def test_region_with_no_reads(data_dir):
    r = bam_stats(str(data_dir / "test.bam"), region="chr1:5001-6000")
    assert r["reads_in_region"] == 0
    assert r["mean_depth"] == 0.0


def test_whole_contig_region(data_dir):
    r = bam_stats(str(data_dir / "test.bam"), region="chr1")
    assert r["reads_in_region"] == 20


def test_wrong_contig_name_lists_available_contigs(data_dir):
    # The BAM names the contig "chr1", not "1"
    with pytest.raises(ValueError, match=r"Contig '1' is not in the BAM header.*\['chr1'\]"):
        bam_stats(str(data_dir / "test.bam"), region="1:1001-1100")


def test_malformed_region_raises(data_dir):
    with pytest.raises(SamtoolsError):
        bam_stats(str(data_dir / "test.bam"), region="chr1:200-100")


HEADER = ("chr1", "chr2", "HLA-A*01:01:01:01")


@pytest.mark.parametrize(
    "region, contig",
    [
        ("chr1:1001-1100", "chr1"),
        ("chr1:1,001-1,100", "chr1"),
        ("chr1:1001", "chr1"),
        ("chr1", "chr1"),
        ("HLA-A*01:01:01:01", "HLA-A*01:01:01:01"),
        ("HLA-A*01:01:01:01:1-100", "HLA-A*01:01:01:01"),
        ("1:1001-1100", "1"),
    ],
)
def test_region_contig_parsing(region, contig):
    assert _region_contig(region, HEADER) == contig
