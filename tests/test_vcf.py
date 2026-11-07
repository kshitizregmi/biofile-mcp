"""Tests for :func:`biofile.vcf_stats`."""
import pytest

from biofile import vcf_stats


def test_counts_match_known_answers(data_dir):
    r = vcf_stats(str(data_dir / "test.vcf"))
    assert r["samples"] == 2
    assert r["records"] == 5
    assert r["SNPs"] == 3
    assert r["indels"] == 2
    assert r["multiallelic sites"] == 0


def test_transitions_transversions(data_dir):
    r = vcf_stats(str(data_dir / "test.vcf"))
    assert (r["ts"], r["tv"]) == (2, 1)
    assert r["ts_tv"] == 2.0


def test_compressed_vcf_gives_same_result(data_dir):
    plain = vcf_stats(str(data_dir / "test.vcf"))
    gz = vcf_stats(str(data_dir / "test.vcf.gz"))
    plain.pop("file"), gz.pop("file")
    assert plain == gz


def test_missing_file_raises(data_dir):
    with pytest.raises(Exception):
        vcf_stats(str(data_dir / "does_not_exist.vcf"))
