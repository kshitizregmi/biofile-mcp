"""Tests for :func:`biofile.bam_stats`."""
import pytest
from pysam.utils import SamtoolsError

from biofile import bam_depth_profile, bam_qc, bam_stats
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


# bam_depth_profile

def test_depth_profile_matches_known_coverage(data_dir):
    r = bam_depth_profile(str(data_dir / "test.bam"), "chr1:1001-1100")
    assert r["region"] == "chr1:1001-1100"
    assert len(r["depth"]) == 100
    assert r["depth"][:50] == [20.0] * 50   # 20 reads over 1001-1050; duplicates excluded
    assert r["depth"][50:] == [0.0] * 50
    assert r["mean_depth"] == 10.0
    assert r["pct_bases_at_min_depth"] == 50.0


def test_depth_profile_bins_long_regions(data_dir):
    r = bam_depth_profile(str(data_dir / "test.bam"), "chr1:1-20000", max_points=500)
    assert len(r["depth"]) == 500
    assert r["position_start"][1] == 41 and r["position_end"][-1] == 20000
    assert r["mean_depth"] == round(20 * 50 / 20000, 2)


def test_depth_profile_open_ended_region_runs_to_the_contig_end(data_dir):
    # chr1:1001 means 1001 to the end of chr1 (248 Mb), which exceeds the limit
    with pytest.raises(ValueError, match="limit is 10,000,000 bp"):
        bam_depth_profile(str(data_dir / "test.bam"), "chr1:1001")


@pytest.mark.parametrize(
    "region, message",
    [
        ("1:1001-1100", r"Contig '1' is not in the BAM header"),
        ("chr1:200-100", "Invalid region"),
        ("chr1", "limit is 10,000,000 bp"),
    ],
)
def test_depth_profile_rejects_bad_regions(data_dir, region, message):
    with pytest.raises(ValueError, match=message):
        bam_depth_profile(str(data_dir / "test.bam"), region)


# bam_qc

def test_qc_summary_matches_known_reads(data_dir):
    r = bam_qc(str(data_dir / "test.bam"))
    s = r["summary"]
    assert (s["reads"], s["mapped"], s["duplicates"]) == (25, 23, 3)
    assert s["mapped_pct"] == 92.0
    assert s["paired"] == 0 and s["properly_paired_pct"] == 0.0
    assert (s["avg_length"], s["max_length"]) == (50.0, 50)
    assert s["error_rate"] == 0.0


def test_qc_histograms(data_dir):
    r = bam_qc(str(data_dir / "test.bam"))
    assert r["mapq"] == {"mapq": [60], "reads": [20]}  # duplicates not counted
    assert r["insert_size"]["count"] == []              # single-end reads
    assert r["read_length"]["start"] == [50]
    # samtools stats counts duplicates in coverage: 23 reads over 50 bases
    assert r["coverage"]["depth"] == [23] and r["coverage"]["bases"] == [50]
    assert r["coverage"]["pct_bases_at_least"]["20x"] == 100.0
    assert r["coverage"]["pct_bases_at_least"]["30x"] == 0.0


def test_qc_contigs_from_index(data_dir):
    r = bam_qc(str(data_dir / "test.bam"))
    assert r["contigs"] == [{"contig": "chr1", "length": 248956422, "mapped": 23, "unmapped": 0}]
    assert r["unplaced_unmapped"] == 2


def test_qc_region_and_bad_contig(data_dir):
    assert bam_qc(str(data_dir / "test.bam"), region="chr1:1001-1100")["region"] == "chr1:1001-1100"
    with pytest.raises(ValueError, match="Contig '1' is not in the BAM header"):
        bam_qc(str(data_dir / "test.bam"), region="1:1001-1100")
