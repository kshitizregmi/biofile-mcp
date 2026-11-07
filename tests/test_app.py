"""Tests for the Streamlit interface, run headless with Streamlit's AppTest."""
from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest  # noqa: E402

APP = str(Path(__file__).resolve().parents[1] / "app.py")


def summarize(example: str, configure=None) -> AppTest:
    at = AppTest.from_file(APP, default_timeout=60).run()
    at.sidebar.selectbox[0].select(example).run()
    if configure:
        configure(at)
    at.button[0].click().run()
    assert not at.exception
    return at


def metrics(at: AppTest) -> dict:
    return {m.label: m.value for m in at.metric}


def test_vcf_view():
    at = summarize("test.vcf")
    assert metrics(at) == {"Records": "5", "SNPs": "3", "Indels": "2", "Ts/Tv": "2.0",
                           "Samples": "2"}
    assert not at.error


def test_bam_view_with_region():
    at = summarize("test.bam", lambda a: a.text_input[0].input("chr1:1001-1100").run())
    assert metrics(at)["Mapped"] == "92.0%"
    assert not at.error


def test_bam_view_reports_wrong_contig():
    at = summarize("test.bam", lambda a: a.text_input[0].input("1:1-100").run())
    assert "Contig '1' is not in the BAM header" in at.error[0].value


def test_fastq_view():
    at = summarize("sample_R1.fastq")
    assert metrics(at)["Reads"] == "1,000"
    assert not at.error


def test_h5ad_view_with_column():
    at = summarize("small.h5ad", lambda a: a.selectbox[0].select("cell_type").run())
    assert metrics(at) == {"Cells": "200", "Genes": "100", "Matrix": "Raw counts",
                           "Embeddings": "1"}


def test_missing_path_is_reported():
    at = AppTest.from_file(APP, default_timeout=60).run()
    at.sidebar.radio[0].set_value("Path on this computer").run()
    at.sidebar.text_input[0].input("/no/such/file.bam").run()
    assert at.sidebar.error[0].value == "No file at that path."


def test_upload_only_mode_hides_paths(monkeypatch):
    monkeypatch.setenv("BIOFILE_UPLOAD_ONLY", "1")
    at = AppTest.from_file(APP, default_timeout=60).run()
    assert "Path on this computer" not in at.sidebar.radio[0].options
