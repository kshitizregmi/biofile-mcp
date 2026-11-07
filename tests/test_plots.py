"""Tests for the Plotly figures in :mod:`plots`."""
import pytest

import plots
from biofile import (
    bam_depth_profile,
    bam_stats,
    fastq_quality_profile,
    h5ad_summary,
    vcf_stats,
)


@pytest.fixture(params=["light", "dark"])
def theme(request):
    return request.param


def title_of(fig) -> str:
    return fig.layout.title.text


def test_quality_profile_figure(data_dir, theme):
    profile = fastq_quality_profile(str(data_dir / "sample_R1.fastq"))
    fig = plots.quality_profile_figure(profile, theme)
    assert "Base quality by read position" in title_of(fig)
    assert [t.name for t in fig.data if t.showlegend is not False][-1] == "Median"
    # Q30 drop is annotated at the position the profile reports
    notes = [a.text for a in fig.layout.annotations]
    assert any(f"position {profile['median_below_q30_from']}" in n for n in notes)
    assert {"Q30", "Q20"} <= set(notes)
    fig.to_json()


def test_read_length_figure(data_dir, theme):
    profile = fastq_quality_profile(str(data_dir / "sample_R1.fastq"))
    fig = plots.read_length_figure(profile, theme)
    assert sum(fig.data[0].y) == 1000
    fig.to_json()


def test_depth_figure_marks_low_coverage(data_dir, theme):
    profile = bam_depth_profile(str(data_dir / "test.bam"), "chr1:901-1200")
    fig = plots.depth_figure(profile, theme)
    assert "chr1:901–1,200" in title_of(fig)
    assert len(fig.layout.shapes) >= 2  # threshold line + shaded low-coverage runs
    assert any(t.name == "Below 20×" for t in fig.data)
    fig.to_json()


def test_read_breakdown_segments_add_up(data_dir, theme):
    stats = bam_stats(str(data_dir / "test.bam"))
    fig = plots.read_breakdown_figure(stats, theme)
    assert [t.name for t in fig.data] == ["Mapped", "Duplicates", "Unmapped"]
    assert sum(t.x[0] for t in fig.data) == stats["total_reads"]


def test_cell_counts_figure_sorts_and_folds(theme):
    counts = {f"type {i}": i + 1 for i in range(30)}
    fig = plots.cell_counts_figure(counts, "cell_type", theme, max_categories=10)
    names = list(fig.data[0].y)
    assert len(names) == 10
    assert names[0].startswith("Other (21 categories)")  # drawn at the bottom
    assert names[-1] == "type 29"                        # largest at the top
    assert sum(fig.data[0].x) == sum(counts.values())


def test_cell_counts_figure_from_h5ad(data_dir, theme):
    result = h5ad_summary(str(data_dir / "small.h5ad"), column="cell_type")
    fig = plots.cell_counts_figure(result["cells_per_cell_type"], "cell_type", theme)
    assert list(fig.data[0].x) == [40, 60, 100]


def test_substitution_figure_collapses_strands(data_dir, theme):
    stats = vcf_stats(str(data_dir / "test.vcf"))
    fig = plots.substitution_figure(stats, theme)
    by_class = {x: y for trace in fig.data for x, y in zip(trace.x, trace.y)}
    # A>G is counted as T>C, G>C as C>G
    assert by_class == {"C>T": 1, "T>C": 1, "C>A": 0, "C>G": 1, "T>A": 0, "T>G": 0}
    assert [t.name for t in fig.data] == ["Transition", "Transversion"]


def test_variant_types_figure(data_dir, theme):
    stats = vcf_stats(str(data_dir / "test.vcf"))
    fig = plots.variant_types_figure(stats, theme)
    assert "Variants by type" in title_of(fig)
    assert dict(zip(fig.data[0].y, fig.data[0].x)) == {"SNPs": 3, "Indels": 2, "MNPs": 0, "Other": 0}


# Detailed views

def test_bam_qc_figures(data_dir, theme):
    from biofile import bam_qc
    qc = bam_qc(str(data_dir / "test.bam"))
    assert "100.0% of mapped reads at MAPQ 30" in title_of(plots.mapq_figure(qc, theme))
    assert plots.insert_size_figure(qc, theme) is None  # single-end test reads
    cov = plots.coverage_distribution_figure(qc, theme)
    assert cov.layout.xaxis.type == "log"
    # On a log axis annotation x is log10(depth): 20x sits at 1.3, not at 20
    xs = [a.x for a in cov.layout.annotations]
    assert all(x < 4 for x in xs)
    assert plots.contigs_figure(qc, theme).data[0].x[0] == 23
    assert plots.bam_read_length_figure(qc, theme) is not None


def test_fastq_detail_figures(data_dir, theme):
    profile = fastq_quality_profile(str(data_dir / "sample_R1.fastq"))
    comp = plots.base_composition_figure(profile, theme)
    assert [t.name for t in comp.data] == ["A", "C", "G", "T", "N"]
    quality = plots.read_quality_figure(profile, theme)
    assert sum(quality.data[0].y) == 1000
    assert quality.data[0].width < 1  # one bar per quality value, not merged
    gc = plots.read_gc_figure(profile, theme)
    assert sum(gc.data[0].y) == 1000


def test_vcf_detail_figures(data_dir, theme):
    stats = vcf_stats(str(data_dir / "test.vcf"))
    assert plots.qual_figure(stats, theme) is not None
    indels = plots.indel_length_figure(stats, theme)
    assert [t.name for t in indels.data] == ["Deletion", "Insertion"]
    samples = plots.per_sample_figure(stats, theme)
    assert {t.name for t in samples.data} == {"Hom-ref", "Het", "Hom-alt", "Missing"}
    assert plots.genotype_depth_figure(stats, theme) is None  # test VCF has no DP


def test_h5ad_detail_figures(data_dir, theme):
    from biofile import h5ad_embedding
    emb = h5ad_embedding(str(data_dir / "small.h5ad"), color="cell_type")
    fig = plots.embedding_figure(emb, theme, highlight="B cell")
    names = [t.name for t in fig.data]
    assert names == ["Other cells", "B cell"]
    assert len(fig.data[1].x) == 60
    labels = {a.text for a in fig.layout.annotations}
    assert "<b>B cell</b>" in labels and "T cell" in labels
    plain = plots.embedding_figure(h5ad_embedding(str(data_dir / "small.h5ad")), theme)
    assert len(plain.data) == 1
