"""Tests for :func:`biofile.fastq_stats`."""
import pytest

from biofile import fastq_stats


def test_read_and_base_counts(data_dir):
    r = fastq_stats(str(data_dir / "sample_R1.fastq"))
    assert r["reads"] == 1000
    assert r["sampled"] is False
    assert (r["min_len"], r["max_len"]) == (240, 251)
    assert r["n_bases"] == 20          # one N in every 50th read
    assert r["lowercase_bases"] == 40  # "acgt" in every 100th read
    assert r["other_characters"] == 0
    assert sum(r["base_counts"].values()) == r["total_bases"]


def test_quality_and_gc_are_percentages(data_dir):
    r = fastq_stats(str(data_dir / "sample_R1.fastq"))
    assert 0 <= r["q30_pct"] <= r["q20_pct"] <= 100
    assert 0 <= r["gc_pct"] <= 100


def test_gzipped_fastq_gives_same_result(data_dir):
    plain = fastq_stats(str(data_dir / "sample_R1.fastq"))
    gz = fastq_stats(str(data_dir / "sample_R1.fastq.gz"))
    plain.pop("file"), gz.pop("file")
    assert plain == gz


def test_max_reads_samples_the_file(data_dir):
    r = fastq_stats(str(data_dir / "sample_R1.fastq"), max_reads=100)
    assert r["reads"] == 100
    assert r["sampled"] is True


def test_empty_file(tmp_path):
    empty = tmp_path / "empty.fastq"
    empty.write_text("")
    r = fastq_stats(str(empty))
    assert r["reads"] == 0
    assert r["warning"] == "file is empty"


def test_non_fastq_file_raises(tmp_path):
    bad = tmp_path / "notes.fastq"
    bad.write_text("this is not a fastq file\n")
    with pytest.raises(RuntimeError):
        fastq_stats(str(bad))
