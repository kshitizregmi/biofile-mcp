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


def test_substitution_counts(data_dir):
    subs = vcf_stats(str(data_dir / "test.vcf"))["substitutions"]
    assert len(subs) == 12
    assert (subs["A>G"], subs["C>T"], subs["G>C"]) == (1, 1, 1)
    assert sum(subs.values()) == 3


def test_compressed_vcf_gives_same_result(data_dir):
    plain = vcf_stats(str(data_dir / "test.vcf"))
    gz = vcf_stats(str(data_dir / "test.vcf.gz"))
    plain.pop("file"), gz.pop("file")
    assert plain == gz


def test_missing_file_raises(data_dir):
    with pytest.raises(Exception):
        vcf_stats(str(data_dir / "does_not_exist.vcf"))


def test_qual_and_indel_lengths(data_dir):
    r = vcf_stats(str(data_dir / "test.vcf"))
    assert sum(r["qual"]["count"]) == 5
    assert (r["qual"]["start"][0], r["qual"]["end"][-1]) == (20, 70)
    assert r["indel_lengths"] == {"length": [-1, 2], "sites": [1, 1]}


def test_per_sample_genotypes(data_dir):
    by_sample = {s["sample"]: s for s in vcf_stats(str(data_dir / "test.vcf"))["per_sample"]}
    assert set(by_sample) == {"TUMOR", "NORMAL"}
    assert (by_sample["TUMOR"]["het"], by_sample["TUMOR"]["hom_alt"]) == (2, 1)
    assert (by_sample["NORMAL"]["hom_ref"], by_sample["NORMAL"]["het"]) == (3, 1)
    assert by_sample["NORMAL"]["missing"] == 1


def test_sites_only_vcf(tmp_path):
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\t.\tA\tG\t.\tPASS\t.\n"
        "chr1\t20\t.\tAT\tA\t30\tPASS\t.\n"
    )
    r = vcf_stats(str(vcf))
    assert (r["samples"], r["records"]) == (0, 2)
    assert r["per_sample"] == []
    assert sum(r["qual"]["count"]) == 1  # the record without QUAL is left out
