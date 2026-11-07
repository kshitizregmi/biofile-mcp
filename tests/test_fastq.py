"""Tests for :func:`biofile.fastq_stats`."""
import pytest

from biofile import fastq_quality_profile, fastq_stats


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


# fastq_quality_profile

def test_quality_profile_is_per_base_for_short_reads(data_dir):
    r = fastq_quality_profile(str(data_dir / "sample_R1.fastq"))
    assert r["reads"] == 1000
    assert r["position_start"][:3] == [1, 2, 3]
    assert r["position_start"] == r["position_end"]  # one position per bin
    assert len(r["median"]) == 251                    # longest read
    keys = ["mean", "median", "q1", "q3", "p10", "p90", "reads_covering"]
    assert all(len(r[k]) == 251 for k in keys)


def test_quality_profile_finds_the_quality_drop(data_dir):
    # Test reads are Q37 for the first 70% of their length, then Q10-30
    r = fastq_quality_profile(str(data_dir / "sample_R1.fastq"))
    assert r["median"][0] == 37
    assert 168 <= r["median_below_q30_from"] <= 177
    assert r["median_below_q20_from"] > r["median_below_q30_from"]
    assert all(lo <= mid <= hi for lo, mid, hi in zip(r["q1"], r["median"], r["q3"]))


def test_quality_profile_counts_fall_where_reads_end(data_dir):
    r = fastq_quality_profile(str(data_dir / "sample_R1.fastq"))
    assert r["reads_covering"][0] == 1000
    assert r["reads_covering"][-1] < 1000  # not every read reaches 251 bp


def test_read_length_histogram(data_dir):
    r = fastq_quality_profile(str(data_dir / "sample_R1.fastq"))
    assert r["length_start"][0] == 240 and r["length_end"][-1] == 251
    assert sum(r["length_reads"]) == 1000


def test_quality_profile_bins_long_reads(tmp_path):
    long_fastq = tmp_path / "long.fastq"
    long_fastq.write_text("@r1\n" + "A" * 1000 + "\n+\n" + "I" * 1000 + "\n")
    r = fastq_quality_profile(str(long_fastq), max_points=100)
    assert len(r["median"]) == 100
    assert (r["position_start"][1], r["position_end"][0]) == (11, 10)
    assert set(r["median"]) == {40}


def test_quality_profile_sampling_and_empty_file(data_dir, tmp_path):
    assert fastq_quality_profile(str(data_dir / "sample_R1.fastq"), max_reads=10)["sampled"] is True
    empty = tmp_path / "empty.fastq"
    empty.write_text("")
    assert fastq_quality_profile(str(empty))["warning"] == "file is empty"


def test_base_composition_per_position(data_dir):
    r = fastq_quality_profile(str(data_dir / "sample_R1.fastq"))
    pct = r["base_pct"]
    assert set(pct) == {"A", "C", "G", "T", "N"}
    for i in range(len(r["position_start"])):
        assert abs(sum(pct[b][i] for b in pct) - 100) < 0.05
    assert pct["N"][10] > 0  # every 50th read has an N at position 11


def test_per_read_quality_and_gc(data_dir):
    r = fastq_quality_profile(str(data_dir / "sample_R1.fastq"))
    assert sum(r["read_quality"]["reads"]) == 1000
    assert sum(r["read_gc"]["reads"]) == 1000
    assert all(0 <= g <= 100 for g in r["read_gc"]["gc_pct"])
