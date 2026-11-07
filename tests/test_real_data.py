"""Tests on small real public datasets.

The files are downloaded once into ``real_data/`` (ignored by git) and reused
on later runs. These tests are skipped by default; run them with::

    python -m pytest -m real_data

If a file cannot be downloaded (for example, with no internet connection),
the tests that need it are skipped rather than failed. Files already in
``real_data/`` are reused, so they can also be downloaded beforehand with the
``curl`` command in the README.

Sources:

- 16S amplicon reads: nf-core/ampliseq test data
- Human reads on part of chr22: nf-core/test-datasets (modules branch)
- NA12878 chr22 variants (Genome in a Bottle): nf-core/test-datasets
- PBMC 3k single-cell data: CZ CELLxGENE example dataset
"""
import asyncio
import json
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

import pytest
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

from biofile import (
    bam_depth_profile,
    bam_stats,
    detect_format,
    fastq_quality_profile,
    fastq_stats,
    h5ad_summary,
    vcf_stats,
)

pytestmark = pytest.mark.real_data

PROJECT_DIR = Path(__file__).resolve().parents[1]
REAL_DIR = Path(os.environ.get("BIOFILE_REAL_DATA", PROJECT_DIR / "real_data"))

NF_CORE = "https://raw.githubusercontent.com/nf-core/test-datasets"
URLS = {
    "1_S103_L001_R1_001.fastq.gz": f"{NF_CORE}/ampliseq/testdata/1_S103_L001_R1_001.fastq.gz",
    "test.paired_end.sorted.bam":
        f"{NF_CORE}/modules/data/genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam",
    "test.paired_end.sorted.bam.bai":
        f"{NF_CORE}/modules/data/genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam.bai",
    "NA12878_GIAB.chr22.vcf.gz":
        f"{NF_CORE}/modules/data/genomics/homo_sapiens/illumina/vcf/NA12878_GIAB.chr22.vcf.gz",
    "pbmc3k.h5ad":
        "https://raw.githubusercontent.com/chanzuckerberg/cellxgene/main/example-dataset/pbmc3k.h5ad",
}


def download(url: str, dest: Path) -> None:
    """Download url to dest, preferring curl.

    curl falls back quickly between IPv6 and IPv4 and uses the system's proxy
    and certificate settings, where Python's urllib can stall for minutes on
    some networks. urllib is only used when curl is not installed.
    """
    if shutil.which("curl"):
        subprocess.run(
            ["curl", "-fsSL", "--retry", "2", "--connect-timeout", "20", "-o", str(dest), url],
            check=True, timeout=600, capture_output=True,
        )
    else:
        with urllib.request.urlopen(url, timeout=30) as response:
            dest.write_bytes(response.read())


def fetch(name: str) -> str:
    """Return the local path of a real data file, downloading it if needed."""
    target = REAL_DIR / name
    if not target.exists():
        REAL_DIR.mkdir(parents=True, exist_ok=True)
        partial = target.with_name(target.name + ".part")
        try:
            download(URLS[name], partial)
        except (OSError, subprocess.SubprocessError) as e:
            partial.unlink(missing_ok=True)
            pytest.skip(f"could not download {name}: {e}")
        partial.rename(target)  # only complete downloads get the real name
    return str(target)


@pytest.fixture(scope="module")
def fastq():
    return fetch("1_S103_L001_R1_001.fastq.gz")


@pytest.fixture(scope="module")
def bam():
    fetch("test.paired_end.sorted.bam.bai")
    return fetch("test.paired_end.sorted.bam")


@pytest.fixture(scope="module")
def vcf():
    return fetch("NA12878_GIAB.chr22.vcf.gz")


@pytest.fixture(scope="module")
def h5ad():
    return fetch("pbmc3k.h5ad")


# Format detection

def test_formats_are_detected(fastq, bam, vcf, h5ad):
    assert [detect_format(p) for p in (fastq, bam, vcf, h5ad)] == ["FASTQ", "BAM", "VCF", "h5ad"]


# FASTQ: 16S amplicon reads

def test_16s_fastq_stats(fastq):
    r = fastq_stats(fastq)
    assert r["reads"] == 2500
    assert (r["min_len"], r["max_len"]) == (35, 250)
    assert r["total_bases"] == 487364
    assert r["q30_pct"] == 83.62
    assert r["gc_pct"] == 51.42
    assert r["n_bases"] == 21805
    assert r["lowercase_bases"] == 0 and r["other_characters"] == 0


def test_16s_quality_profile(fastq):
    r = fastq_quality_profile(fastq)
    assert r["reads"] == 2500
    assert len(r["median"]) == 250  # per-base for reads up to 250 bp
    assert r["median"][0] == 32
    assert r["median_below_q30_from"] == 250
    assert sum(r["length_reads"]) == 2500


# BAM: human reads on part of chr22

def test_bam_read_counts(bam):
    r = bam_stats(bam)
    assert r["total_reads"] == 5644
    assert r["mapped_reads"] == 5642
    assert r["mapped_pct"] == 99.96
    assert r["duplicates"] == 0


def test_bam_region_coverage(bam):
    r = bam_stats(bam, region="chr22:2000-6000")
    assert r["reads_in_region"] == 5640
    assert r["mean_depth"] == pytest.approx(165.265)


def test_bam_depth_profile_across_the_contig(bam):
    r = bam_depth_profile(bam, "chr22", min_depth=10)
    assert r["region"] == "chr22:1-40001"
    assert len(r["depth"]) == 494  # 40,001 bp in 81 bp bins
    assert r["mean_depth"] == 16.77
    assert r["max_depth"] == 2532


def test_bam_region_with_wrong_contig_name(bam):
    with pytest.raises(ValueError, match=r"Contigs in this file: \['chr22'\]"):
        bam_stats(bam, region="22:2000-6000")


# VCF: NA12878 chr22 variants

def test_giab_vcf_stats(vcf):
    r = vcf_stats(vcf)
    assert r["samples"] == 1
    assert r["records"] == 81
    assert (r["SNPs"], r["indels"]) == (70, 11)
    assert (r["ts"], r["tv"], r["ts_tv"]) == (52, 18, 2.89)
    assert r["multiallelic sites"] == 3
    assert sum(r["substitutions"].values()) == 70


# h5ad: PBMC 3k

def test_pbmc3k_summary(h5ad):
    r = h5ad_summary(h5ad, column="louvain")
    assert (r["n_cells"], r["n_genes"]) == (2638, 1838)
    assert r["X_looks_like"] == "normalized or scaled"
    counts = r["cells_per_louvain"]
    assert len(counts) == 8
    assert counts["CD4 T cells"] == 1144
    assert counts["Megakaryocytes"] == 15
    assert sum(counts.values()) == 2638


# The MCP server on real data

def test_mcp_server_on_real_data(fastq, bam, vcf, h5ad):
    params = StdioServerParameters(
        command=sys.executable,
        args=[str(PROJECT_DIR / "server.py")],
        cwd=str(PROJECT_DIR),
        env={**os.environ, "BIOFILE_ALLOWED_DIR": str(REAL_DIR.resolve())},
    )
    calls = [
        ("fastq_stats", {"path": Path(fastq).name}),
        ("bam_stats", {"path": Path(bam).name, "region": "chr22:2000-6000"}),
        ("vcf_stats", {"path": Path(vcf).name}),
        ("h5ad_summary", {"path": Path(h5ad).name, "column": "louvain"}),
    ]

    async def run():
        with open(os.devnull, "w") as quiet:
            async with stdio_client(params, errlog=quiet) as (read, write):
                async with ClientSession(read, write) as session:
                    await session.initialize()
                    return [await session.call_tool(name, args) for name, args in calls]

    results = asyncio.run(run())
    assert not any(r.is_error for r in results)
    fq, bm, vc, ad = (json.loads(r.content[0].text) for r in results)
    assert fq["reads"] == 2500
    assert bm["reads_in_region"] == 5640
    assert vc["ts_tv"] == 2.89
    assert ad["cells_per_louvain"]["B cells"] == 342


# The web interface on real data

def test_web_interface_on_real_data(h5ad):
    from streamlit.testing.v1 import AppTest

    at = AppTest.from_file(str(PROJECT_DIR / "app.py"), default_timeout=120).run()
    at.sidebar.radio[0].set_value("Path on this computer").run()
    at.sidebar.text_input[0].input(h5ad).run()
    at.selectbox[0].select("louvain").run()
    at.button[0].click().run()
    assert not at.exception and not at.error
    shown = {m.label: m.value for m in at.metric}
    assert (shown["Cells"], shown["Genes"]) == ("2,638", "1,838")
    assert shown["Matrix"] == "Normalized"


# Detailed statistics on real data

def test_bam_qc_on_real_reads(bam):
    from biofile import bam_qc
    r = bam_qc(bam)
    s = r["summary"]
    assert s["properly_paired_pct"] == 99.93
    assert (s["insert_size_mean"], s["insert_size_sd"]) == (124.8, 31.2)
    assert s["error_rate"] == 0.00201493
    assert r["coverage"]["pct_bases_at_least"]["20x"] == 79.34
    assert r["mapq"]["reads"][-1] == 5635  # almost every read at MAPQ 60


def test_vcf_detail_on_real_variants(vcf):
    r = vcf_stats(vcf)
    [sample] = r["per_sample"]
    assert (sample["sample"], sample["het"], sample["hom_alt"]) == ("NA12878", 32, 38)
    assert sum(r["indel_lengths"]["sites"]) == 13


def test_fastq_detail_on_16s_reads(fastq):
    r = fastq_quality_profile(fastq)
    by_quality = dict(zip(r["read_quality"]["quality"], r["read_quality"]["reads"]))
    assert by_quality[2] == 623  # a quarter of the reads are all-N, quality 2
    assert sum(r["read_gc"]["reads"]) == 2500


def test_h5ad_detail_on_pbmc3k(h5ad):
    from biofile import h5ad_column_stats, h5ad_embedding
    genes = h5ad_column_stats(h5ad, "n_genes")
    assert (genes["median"], genes["min"], genes["max"]) == (820.0, 212.0, 2455.0)
    emb = h5ad_embedding(h5ad, color="louvain")
    assert emb["basis"] == "X_umap" and len(emb["x"]) == 2638
