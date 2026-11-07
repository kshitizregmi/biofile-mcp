"""Web interface to the biofile summaries.

Usage::

    streamlit run app.py

Pick a file (upload it, give a path on this computer, or use the built-in
example data), set the options for that file type and press Summarize.

Set ``BIOFILE_UPLOAD_ONLY=1`` when the app is hosted for other people, so
visitors cannot read files from the server by path.
"""
import json
import os
import tempfile
from pathlib import Path

import pandas as pd
import pysam
import streamlit as st

import plots
from biofile import (
    bam_depth_profile,
    bam_qc,
    bam_stats,
    detect_format,
    fastq_quality_profile,
    fastq_stats,
    h5ad_column_stats,
    h5ad_embedding,
    h5ad_summary,
    vcf_stats,
)
from biofile.formats import is_index
from make_test_data import make_test_data

UPLOAD_ONLY = os.environ.get("BIOFILE_UPLOAD_ONLY") == "1"

# Reads used for the quality plot when "all reads" is chosen; enough for a
# stable profile without reading a whole sequencing run
PROFILE_READS = 100_000


@st.cache_resource
def example_folder() -> Path:
    """Build the example files once per app session."""
    return make_test_data(Path(tempfile.mkdtemp()) / "example_data")


def save_uploads(files) -> Path:
    """Copy uploaded files into this browser session's own temporary folder.

    Each file is written once; later reruns reuse the copy. Original names are
    kept because samtools finds a BAM's index by looking for ``<name>.bam.bai``.
    Index files are given a newer modification time than their data files, so
    htslib does not warn that the index is older than the file it indexes.
    """
    if "upload_dir" not in st.session_state:
        st.session_state.upload_dir = tempfile.mkdtemp(prefix="biofile_")
    folder = Path(st.session_state.upload_dir)
    written = False
    for f in files:
        target = folder / f.name
        if not target.exists() or target.stat().st_size != f.size:
            target.write_bytes(f.getbuffer())
            written = True

    # Files are written in upload order, so an index can end up older than
    # its BAM or VCF. Touch the indexes last to put them back in order.
    if written:
        for index in folder.iterdir():
            if is_index(index.name):
                index.touch()
    return folder


def has_bam_index(path: str) -> bool:
    return any(Path(path + suffix).exists() for suffix in (".bai", ".csi")) or \
        Path(path).with_suffix(".bai").exists()


FUNCTIONS = {
    "vcf_stats": vcf_stats,
    "bam_stats": bam_stats,
    "bam_qc": bam_qc,
    "bam_depth_profile": bam_depth_profile,
    "fastq_stats": fastq_stats,
    "fastq_quality_profile": fastq_quality_profile,
    "h5ad_summary": h5ad_summary,
    "h5ad_column_stats": h5ad_column_stats,
    "h5ad_embedding": h5ad_embedding,
}


@st.cache_data(show_spinner=False, max_entries=64)
def _cached(name: str, path: str, modified: float, *args, **kwargs):
    return FUNCTIONS[name](path, *args, **kwargs)


def cached(name: str, path: str, *args, **kwargs):
    """Call a biofile function, reusing the result until the file changes."""
    return _cached(name, path, os.path.getmtime(path), *args, **kwargs)


def chart(fig) -> None:
    st.plotly_chart(fig, theme=None, config=plots.CONFIG, width="stretch")


# Page header
st.set_page_config(page_title="biofile", page_icon="🧬", layout="wide")
theme = "dark" if st.context.theme.type == "dark" else "light"

st.title("🧬 biofile")
st.caption("Summarize VCF, BAM, FASTQ and h5ad files. Files are only read, never changed.")

# Step 1: choose a file
sources = ["Example data", "Upload files"] + ([] if UPLOAD_ONLY else ["Path on this computer"])
source = st.sidebar.radio("File source", sources)

path = None
if source == "Example data":
    folder = example_folder()
    names = sorted(p.name for p in folder.iterdir() if not is_index(p.name))
    path = str(folder / st.sidebar.selectbox("Example file", names))

elif source == "Upload files":
    uploads = st.sidebar.file_uploader(
        "One or more files. For a BAM, add its .bai index if you have it.",
        accept_multiple_files=True,
    )
    if uploads:
        folder = save_uploads(uploads)
        names = [f.name for f in uploads if not is_index(f.name)]
        formats = {name: detect_format(str(folder / name)) for name in names}

        if not names:
            # Only index files so far (e.g. a .bai on its own, or one that
            # finished uploading before its BAM)
            st.sidebar.info("Only index files are uploaded. Add the BAM or VCF they belong to.")
        else:
            # When files are added, switch to the newest one
            previous = st.session_state.get("uploaded_names", [])
            added = [n for n in names if n not in previous]
            if added or st.session_state.get("picked") not in names:
                st.session_state.picked = (added or names)[-1]

            if len(names) > 1:
                st.radio(
                    "Uploaded files: choose one to summarize",
                    names,
                    key="picked",
                    horizontal=True,
                    format_func=lambda n: f"{n}  ({formats[n] or 'not recognized'})",
                )
            path = str(folder / st.session_state.picked)

            # A BAM uploaded without its index gets one, built on the uploaded
            # copy only; region and depth queries need it
            if formats[st.session_state.picked] == "BAM" and not has_bam_index(path):
                pysam.index(path)
                st.caption("No .bai index was uploaded for this BAM, so one was created.")
        st.session_state.uploaded_names = names

else:
    typed = st.sidebar.text_input("Full path", placeholder="/Users/me/data/sample.bam").strip()
    if typed:
        if Path(typed).expanduser().is_file():
            path = str(Path(typed).expanduser())
        else:
            st.sidebar.error("No file at that path.")

if not path:
    st.info("Choose a file in the sidebar to begin.")
    st.stop()

fmt = detect_format(path)
if fmt is None:
    st.error(f"{Path(path).name} is not a recognized VCF, BAM, FASTQ or h5ad file.")
    st.stop()

st.subheader(f"{Path(path).name}  ·  {fmt}")

# Step 2: options for this file type
try:
    if fmt == "BAM":
        left, right = st.columns([3, 1])
        region = left.text_input("Region (optional)", placeholder="chr1:1001-1100").strip() or None
        min_depth = int(right.number_input("Target depth", min_value=1, value=20, step=5))
        options = (region, min_depth)

    elif fmt == "FASTQ":
        limit = int(st.number_input("Reads to check (0 = all)", min_value=0, value=0, step=100_000))
        options = (limit,)

    elif fmt == "h5ad":
        columns = cached("h5ad_summary", path)["cell_columns"]
        column = st.selectbox("Group cells by", ["(none)"] + columns)
        column = None if column == "(none)" else column
        options = (column,)

    else:
        options = ()

    # Results stay on screen until the file or an option changes, so the
    # controls inside the charts can be used without pressing Summarize again
    run_key = (path, fmt, options)
    if st.button("Summarize", type="primary"):
        st.session_state.run_key = run_key
    if st.session_state.get("run_key") != run_key:
        st.stop()

    if fmt == "BAM" and region and not has_bam_index(path):
        st.error(f"A region needs the BAM's index. Create it with: samtools index {Path(path).name}")
        st.stop()

    # Step 3: read the file (cached, so reruns are instant)
    with st.spinner("Reading file..."):
        if fmt == "VCF":
            result = cached("vcf_stats", path)
        elif fmt == "BAM":
            result = cached("bam_stats", path, region=region)
            qc = cached("bam_qc", path, region=region)
            depth = cached("bam_depth_profile", path, region, min_depth=min_depth) if region else None
        elif fmt == "FASTQ":
            result = cached("fastq_stats", path, max_reads=limit or None)
            profile = (cached("fastq_quality_profile", path, max_reads=limit or PROFILE_READS)
                       if result["reads"] else None)
        else:
            result = cached("h5ad_summary", path, column=column)
            numeric = []
            for c in result["cell_columns"][:12]:
                stats = cached("h5ad_column_stats", path, c)
                if stats["kind"] == "numeric" and "histogram" in stats:
                    numeric.append(stats)
            embedding = (cached("h5ad_embedding", path, color=column)
                         if result["embeddings"] else None)

except Exception as e:
    st.error(f"{type(e).__name__}: {e}")
    st.stop()


def pair(*figures) -> None:
    """Show up to two charts side by side, skipping any that are None."""
    figures = [f for f in figures if f is not None]
    for column, fig in zip(st.columns(len(figures)) if len(figures) > 1 else [st.container()],
                           figures):
        with column:
            chart(fig)


# Step 4: show the results
if "warning" in result:
    st.warning(result["warning"])

if fmt == "VCF":
    cols = st.columns(5)
    cols[0].metric("Records", f"{result['records']:,}")
    cols[1].metric("SNPs", f"{result['SNPs']:,}")
    cols[2].metric("Indels", f"{result['indels']:,}")
    cols[3].metric("Ts/Tv", result.get("ts_tv", "n/a"))
    cols[4].metric("Samples", f"{result['samples']:,}")

    overview, quality, indels, samples = st.tabs(["Overview", "Quality & depth", "Indels", "Samples"])
    with overview:
        pair(plots.substitution_figure(result, theme), plots.variant_types_figure(result, theme))
    with quality:
        pair(plots.qual_figure(result, theme), plots.genotype_depth_figure(result, theme))
    with indels:
        fig = plots.indel_length_figure(result, theme)
        if fig:
            chart(fig)
        else:
            st.info("No indels in this file.")
    with samples:
        fig = plots.per_sample_figure(result, theme)
        if fig:
            chart(fig)
            st.dataframe(pd.DataFrame(result["per_sample"]).set_index("sample"), width="stretch")
        else:
            st.info("This VCF has no sample columns (sites only).")

elif fmt == "BAM":
    summary = qc["summary"]
    cols = st.columns(5)
    cols[0].metric("Total reads", f"{result['total_reads']:,}")
    cols[1].metric("Mapped", f"{result['mapped_pct']}%")
    cols[2].metric("Properly paired", f"{summary['properly_paired_pct']}%" if summary["paired"] else "single-end")
    cols[3].metric("Error rate", f"{100 * summary['error_rate']:.2f}%")
    cols[4].metric("Mean insert", f"{summary['insert_size_mean']:,.0f} bp" if summary["paired"] else "n/a")
    if region:
        st.caption(f"Detailed statistics below are for {region}.")

    overview, coverage, mapq, inserts, lengths = st.tabs(
        ["Overview", "Coverage", "Mapping quality", "Insert size", "Read length"])
    with overview:
        chart(plots.read_breakdown_figure(result, theme))
        fig = plots.contigs_figure(qc, theme)
        if fig:
            chart(fig)
        details = {
            "Primary reads": summary["reads"], "Secondary": summary["secondary"],
            "Supplementary": summary["supplementary"], "Duplicates": summary["duplicates"],
            "MAPQ 0": summary["mq0"], "Bases mapped": summary["bases_mapped"],
            "Mismatches": summary["mismatches"], "Average base quality": summary["avg_quality"],
            "Pairs on different contigs": summary["pairs_on_different_contigs"],
        }
        st.dataframe(pd.Series(details, name="value").to_frame(), width="content")
    with coverage:
        if depth:
            chart(plots.depth_figure(depth, theme))
        else:
            st.info("Enter a region above (for example a contig name) to see read depth along it.")
        fig = plots.coverage_distribution_figure(qc, theme)
        if fig:
            chart(fig)
    with mapq:
        chart(plots.mapq_figure(qc, theme))
    with inserts:
        fig = plots.insert_size_figure(qc, theme)
        if fig:
            chart(fig)
        else:
            st.info("No read pairs: this looks like single-end data.")
    with lengths:
        fig = plots.bam_read_length_figure(qc, theme)
        if fig:
            chart(fig)

elif fmt == "FASTQ" and result["reads"]:
    cols = st.columns(4)
    cols[0].metric("Reads", f"{result['reads']:,}" + (" (sample)" if result["sampled"] else ""))
    cols[1].metric("Avg length", f"{result['avg_len']} bp")
    cols[2].metric("Q30", f"{result['q30_pct']}%")
    cols[3].metric("GC", f"{result['gc_pct']}%")

    quality, composition, lengths = st.tabs(["Quality", "Composition", "Read length"])
    with quality:
        chart(plots.quality_profile_figure(profile, theme))
        chart(plots.read_quality_figure(profile, theme))
    with composition:
        chart(plots.base_composition_figure(profile, theme))
        left, right = st.columns([3, 2])
        with left:
            chart(plots.read_gc_figure(profile, theme))
        with right:
            st.markdown("**Base counts**")
            # A table rather than a chart: N and lowercase counts are often thousands
            # of times smaller than A/C/G/T and would not show up as bars
            bases = pd.DataFrame(
                {"count": result["base_counts"].values()}, index=result["base_counts"].keys()
            ).sort_values("count", ascending=False)
            bases["% of bases"] = (100 * bases["count"] / result["total_bases"]).round(3)
            st.dataframe(bases, width="stretch")
    with lengths:
        if len(profile["length_reads"]) > 1:
            chart(plots.read_length_figure(profile, theme))
        else:
            st.info(f"Every read is {profile['length_start'][0]} bp long.")

elif fmt == "h5ad":
    cols = st.columns(4)
    cols[0].metric("Cells", f"{result['n_cells']:,}")
    cols[1].metric("Genes", f"{result['n_genes']:,}")
    cols[2].metric("Matrix", "Raw counts" if result["X_looks_like"] == "raw counts" else "Normalized",
                   help="Raw counts are needed by methods such as DESeq2 or scVI.")
    cols[3].metric("Embeddings", len(result["embeddings"]),
                   help=", ".join(k.removeprefix("X_") for k in result["embeddings"]) or "none")

    groups, qc_tab, embed_tab = st.tabs(["Groups", "QC metrics", "Embedding"])
    with groups:
        if column:
            chart(plots.cell_counts_figure(result[f"cells_per_{column}"], column, theme))
        else:
            st.info("Choose a column in “Group cells by” above to count cells per group.")
    with qc_tab:
        if numeric:
            for i in range(0, len(numeric), 2):
                pair(*(plots.obs_histogram_figure(c, theme) for c in numeric[i:i + 2]))
        else:
            st.info("No numeric cell metadata (such as n_counts or percent_mito) in this file.")
    with embed_tab:
        if embedding is None:
            st.info("This file has no embedding (such as X_umap) to draw.")
        else:
            highlight = None
            if column and embedding.get("values") and isinstance(embedding["values"][0], str):
                group_names = result[f"cells_per_{column}"]
                highlight = st.selectbox("Highlight", list(group_names))
            chart(plots.embedding_figure(embedding, theme, highlight))

with st.expander("Full result (JSON)"):
    st.json(result)

st.download_button(
    "Download JSON",
    data=json.dumps(result, indent=2),
    file_name=f"{Path(path).name}.summary.json",
    mime="application/json",
)
