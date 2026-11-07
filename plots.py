"""Plotly figures for the biofile web interface.

Every function takes the dict returned by a ``biofile`` summary function and a
theme name (``"light"`` or ``"dark"``) and returns a ``plotly.graph_objects.Figure``.
Colors come from one validated palette with separate steps for light and dark
backgrounds; marks carry the color, and all text uses neutral ink colors.
"""
import math

import plotly.graph_objects as go

FONT = 'system-ui, -apple-system, "Segoe UI", sans-serif'

THEMES = {
    "light": {
        "ink": "#0b0b0b",
        "ink_2": "#52514e",
        "muted": "#898781",
        "grid": "#e1e0d9",
        "axis": "#c3c2b7",
        "surface": "#fcfcfb",
        "border": "rgba(11,11,11,0.10)",
        "series": ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"],
        "critical": "#d03b3b",
        "context": "#c3c2b7",
        "sequential": ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b"],
    },
    "dark": {
        "ink": "#f0efec",
        "ink_2": "#c3c2b7",
        "muted": "#898781",
        "grid": "#2c2c2a",
        "axis": "#383835",
        "surface": "#1a1a19",
        "border": "rgba(255,255,255,0.10)",
        "series": ["#3987e5", "#d95926", "#199e70", "#c98500"],
        "critical": "#d03b3b",
        "context": "#52514e",
        "sequential": ["#184f95", "#256abf", "#3987e5", "#6da7ec", "#9ec5f4", "#cde2fb"],
    },
}

# Config for st.plotly_chart: keep the PNG download, drop selection tools
CONFIG = {
    "displaylogo": False,
    "modeBarButtonsToRemove": ["select2d", "lasso2d", "autoScale2d", "toggleSpikelines"],
    "toImageButtonOptions": {"format": "png", "scale": 2},
}

BAR_PX = 22  # target bar thickness; bars never fill their whole slot


def _rgba(hex_color: str, alpha: float) -> str:
    r, g, b = (int(hex_color[i:i + 2], 16) for i in (1, 3, 5))
    return f"rgba({r},{g},{b},{alpha})"


def _base(fig: go.Figure, t: dict, title: str, subtitle: str, height: int) -> go.Figure:
    """Apply the shared look: title block, recessive axes, quiet hover labels."""
    fig.update_layout(
        height=height,
        margin=dict(l=8, r=48, t=80, b=8),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family=FONT, size=13, color=t["ink_2"]),
        title=dict(
            text=f"<b>{title}</b><br><span style='font-size:13px;color:{t['ink_2']}'>{subtitle}</span>",
            x=0, xref="container", xanchor="left", y=1, yref="container", yanchor="top",
            pad=dict(t=16, l=0), font=dict(size=16, color=t["ink"]),
        ),
        hoverlabel=dict(bgcolor=t["surface"], bordercolor=t["border"],
                        font=dict(family=FONT, size=13, color=t["ink"])),
        showlegend=False,
        legend=dict(orientation="h", x=0, xanchor="left", y=1.0, yanchor="bottom",
                    traceorder="normal", font=dict(color=t["ink_2"]), bgcolor="rgba(0,0,0,0)"),
    )
    axis = dict(gridcolor=t["grid"], gridwidth=1, linecolor=t["axis"], zeroline=False,
                tickfont=dict(color=t["muted"]), title_font=dict(color=t["ink_2"]),
                ticks="", automargin=True)
    fig.update_xaxes(**axis)
    fig.update_yaxes(**axis)
    return fig


def _reference_line(fig: go.Figure, t: dict, y: float, label: str) -> None:
    """Horizontal hairline with a small label at its right end."""
    fig.add_hline(y=y, line=dict(color=t["axis"], width=1), layer="below")
    fig.add_annotation(x=1, xref="paper", y=y, text=label, showarrow=False, xshift=6,
                       xanchor="left", yanchor="middle", font=dict(size=12, color=t["muted"]))


def _count_axis(fig: go.Figure, top: float) -> None:
    """Y axis for whole-number counts: no fractional ticks, a little headroom."""
    fig.update_yaxes(range=[0, top * 1.15 if top else 1], tickformat=",d",
                     dtick=1 if top <= 5 else None)


def _bar_gap(n: int, slot_px: float) -> float:
    """Gap fraction that keeps bars near BAR_PX thick inside slots of slot_px."""
    return max(0.15, 1 - BAR_PX / max(slot_px, BAR_PX))


# FASTQ

def quality_profile_figure(profile: dict, theme: str = "light") -> go.Figure:
    """Per-position quality: median line over IQR and 10th-90th percentile bands."""
    t = THEMES[theme]
    blue = t["series"][0]
    binned = profile["position_start"] != profile["position_end"]
    x = [(s + e) / 2 for s, e in zip(profile["position_start"], profile["position_end"])]

    fig = go.Figure()
    # Washes need more opacity to read on a dark background
    outer, inner = (0.10, 0.22) if theme == "light" else (0.16, 0.32)
    for low, high, alpha, name, rank in [("p10", "p90", outer, "10th–90th percentile", 3),
                                         ("q1", "q3", inner, "Middle 50%", 2)]:
        fig.add_trace(go.Scatter(x=x, y=profile[low], mode="lines", line=dict(width=0),
                                 hoverinfo="skip", showlegend=False))
        fig.add_trace(go.Scatter(x=x, y=profile[high], mode="lines", line=dict(width=0),
                                 fill="tonexty", fillcolor=_rgba(blue, alpha),
                                 name=name, legendrank=rank, hoverinfo="skip"))

    position = (
        [f"{s:,}–{e:,}" for s, e in zip(profile["position_start"], profile["position_end"])]
        if binned else [f"{s:,}" for s in profile["position_start"]]
    )
    custom = list(zip(position, profile["mean"], profile["q1"], profile["q3"],
                      profile["p10"], profile["p90"], profile["reads_covering"]))
    fig.add_trace(go.Scatter(
        x=x, y=profile["median"], mode="lines", name="Median", legendrank=1,
        line=dict(color=blue, width=2),
        customdata=custom,
        hovertemplate=(
            "<b>Position %{customdata[0]}</b><br>"
            "Median Q%{y} · mean Q%{customdata[1]:.1f}<br>"
            "Middle 50%: Q%{customdata[2]}–Q%{customdata[3]}<br>"
            "10th–90th: Q%{customdata[4]}–Q%{customdata[5]}<br>"
            "%{customdata[6]:,} bases<extra></extra>"
        ),
    ))

    _reference_line(fig, t, 30, "Q30")
    _reference_line(fig, t, 20, "Q20")

    drop = profile["median_below_q30_from"]
    if drop is not None:
        i = profile["position_start"].index(drop)
        fig.add_trace(go.Scatter(
            x=[x[i]], y=[profile["median"][i]], mode="markers", hoverinfo="skip",
            showlegend=False,
            marker=dict(size=10, color=blue, line=dict(color=t["surface"], width=2)),
        ))
        fig.add_annotation(
            x=x[i], y=profile["median"][i], ax=-70, ay=-58,
            text=f"Median drops below Q30<br>at position {drop:,}",
            showarrow=True, arrowhead=0, arrowwidth=1, arrowcolor=t["muted"],
            align="right", font=dict(size=12, color=t["ink"]),
            bgcolor=t["surface"], borderpad=3,
        )

    reads = f"{profile['reads']:,} reads" + (" (sample)" if profile["sampled"] else "")
    _base(fig, t, "Base quality by read position",
          f"Phred score of every base, summarised per position · {reads}", height=420)
    top = max(41, max(profile["p90"]) + 2)
    fig.update_yaxes(title_text="Phred quality", range=[0, top], dtick=10)
    fig.update_xaxes(title_text="Position in read (bp)", showgrid=False, showline=True,
                     showspikes=True, spikemode="across", spikecolor=t["axis"],
                     spikethickness=1, spikedash="solid")
    fig.update_layout(showlegend=True, hovermode="x", margin=dict(t=104))
    return fig


def read_length_figure(profile: dict, theme: str = "light", plot_px: int = 560) -> go.Figure:
    """Histogram of read lengths with the most common length labelled."""
    t = THEMES[theme]
    starts, ends, reads = profile["length_start"], profile["length_end"], profile["length_reads"]
    centers = [(s + e) / 2 for s, e in zip(starts, ends)]
    bin_width = ends[0] - starts[0] + 1
    labels = [f"{s:,}" if s == e else f"{s:,}–{e:,}" for s, e in zip(starts, ends)]
    total = sum(reads)
    span = ends[-1] - starts[0] + 1
    slot_px = plot_px * bin_width / span

    fig = go.Figure(go.Bar(
        x=centers, y=reads, width=bin_width * min(0.7, BAR_PX / slot_px),
        marker=dict(color=t["series"][0], cornerradius=4),
        customdata=list(zip(labels, [100 * n / total for n in reads])),
        hovertemplate="<b>%{customdata[0]} bp</b><br>%{y:,} reads (%{customdata[1]:.1f}%)<extra></extra>",
    ))
    peak = max(range(len(reads)), key=reads.__getitem__)
    fig.add_annotation(x=centers[peak], y=reads[peak], yshift=2, showarrow=False,
                       yanchor="bottom", text=f"{reads[peak]:,}",
                       font=dict(size=12, color=t["ink_2"]))

    _base(fig, t, "Read length distribution",
          f"Most common: {labels[peak]} bp · range {starts[0]:,}–{ends[-1]:,} bp", height=320)
    fig.update_xaxes(title_text="Read length (bp)", showgrid=False, showline=True,
                     range=[starts[0] - bin_width, ends[-1] + bin_width], tickformat=",d")
    fig.update_yaxes(title_text="Reads")
    _count_axis(fig, max(reads))
    return fig


# BAM

def depth_figure(profile: dict, theme: str = "light") -> go.Figure:
    """Coverage track with the depth threshold and low-coverage stretches marked."""
    t = THEMES[theme]
    blue = t["series"][0]
    threshold = profile["min_depth_threshold"]
    starts, ends, depth = profile["position_start"], profile["position_end"], profile["depth"]
    per_base = starts == ends

    # Shade stretches below the threshold (merged runs of consecutive bins)
    runs, run_start = [], None
    for s, e, d in zip(starts, ends, depth):
        if d < threshold and run_start is None:
            run_start = s
        if d >= threshold and run_start is not None:
            runs.append((run_start, s - 1))
            run_start = None
    if run_start is not None:
        runs.append((run_start, ends[-1]))

    fig = go.Figure()
    for lo, hi in runs[:50]:
        fig.add_vrect(x0=lo - 0.5, x1=hi + 0.5, fillcolor=_rgba(t["critical"], 0.08),
                      line_width=0, layer="below")
    if runs:
        fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers", name=f"Below {threshold}×",
                                 marker=dict(symbol="square", size=12,
                                             color=_rgba(t["critical"], 0.25))))

    fig.add_trace(go.Scatter(
        x=starts, y=depth, mode="lines", name="Depth",
        line=dict(color=blue, width=2, shape="hv" if per_base else "linear"),
        fill="tozeroy", fillcolor=_rgba(blue, 0.10),
        customdata=list(zip(starts, ends)),
        hovertemplate=(
            ("<b>%{customdata[0]:,}</b>" if per_base else "<b>%{customdata[0]:,}–%{customdata[1]:,}</b>")
            + "<br>Depth %{y:.1f}×<extra></extra>"
        ),
    ))
    _reference_line(fig, t, threshold, f"{threshold}× target")

    contig, span = profile["region"].split(":")
    lo, hi = (int(v) for v in span.split("-"))
    _base(fig, t, f"Read depth across {contig}:{lo:,}–{hi:,}",
          f"Mean {profile['mean_depth']:.1f}× · {profile['pct_bases_at_min_depth']:.1f}% "
          f"of bases at or above {threshold}× · duplicates excluded", height=380)
    fig.update_yaxes(title_text="Depth (×)", range=[0, max(max(depth), threshold) * 1.2])
    fig.update_xaxes(title_text=f"Position on {contig}", showgrid=False, showline=True,
                     tickformat=",d", showspikes=True, spikemode="across",
                     spikecolor=t["axis"], spikethickness=1, spikedash="solid")
    fig.update_layout(showlegend=bool(runs), hovermode="x", margin=dict(t=108 if runs else 80, r=88))
    return fig


def read_breakdown_figure(stats: dict, theme: str = "light") -> go.Figure:
    """One stacked bar splitting all reads into mapped, duplicate and unmapped."""
    t = THEMES[theme]
    total = stats["total_reads"]
    duplicates = stats["duplicates"]
    mapped_unique = max(stats["mapped_reads"] - duplicates, 0)
    unmapped = total - stats["mapped_reads"]
    parts = [("Mapped", mapped_unique), ("Duplicates", duplicates), ("Unmapped", unmapped)]

    fig = go.Figure()
    for (name, n), color in zip(parts, t["series"]):
        share = 100 * n / total if total else 0
        fig.add_trace(go.Bar(
            y=["reads"], x=[n], name=name, orientation="h",
            marker=dict(color=color, line=dict(color=t["surface"], width=2)),
            # label inside only when the segment is wide enough to hold it
            text=[f"{share:.0f}%" if share >= 10 else ""], textposition="inside",
            insidetextanchor="middle",
            textfont=dict(color="#ffffff", size=13), cliponaxis=False,
            hovertemplate=f"<b>{name}</b><br>%{{x:,}} reads ({share:.1f}%)<extra></extra>",
        ))
    _base(fig, t, "Where the reads went",
          f"{total:,} QC-passed reads · {stats['mapped_pct']}% mapped", height=178)
    fig.update_layout(barmode="stack", bargap=0.45, showlegend=True, margin=dict(t=108, b=8))
    fig.update_xaxes(visible=False, range=[0, total])
    fig.update_yaxes(visible=False)
    return fig


# h5ad

def ranked_bars_figure(counts: dict, title: str, subtitle: str, unit: str,
                       theme: str = "light", max_categories: int = 25) -> go.Figure:
    """Horizontal bars, largest first, each labelled with its count and share."""
    t = THEMES[theme]
    items = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)
    if len(items) > max_categories:
        kept, rest = items[:max_categories - 1], items[max_categories - 1:]
        items = kept + [(f"Other ({len(rest)} categories)", sum(n for _, n in rest))]
    names = [k for k, _ in items][::-1]  # Plotly draws the first category at the bottom
    values = [v for _, v in items][::-1]
    total = sum(values) or 1

    fig = go.Figure(go.Bar(
        y=names, x=values, orientation="h",
        marker=dict(color=t["series"][0], cornerradius=4),
        text=[f"{v:,}  ·  {100 * v / total:.1f}%" if v else "0" for v in values],
        textposition="outside", textfont=dict(color=t["ink_2"], size=12), cliponaxis=False,
        hovertemplate="<b>%{y}</b><br>%{x:,} " + unit + "<extra></extra>",
    ))
    slot = 34
    _base(fig, t, title, subtitle, height=110 + slot * len(names))
    fig.update_layout(bargap=_bar_gap(len(names), slot), margin=dict(r=110))
    fig.update_xaxes(visible=False, range=[0, max(values + [1]) * 1.02])
    fig.update_yaxes(showgrid=False, tickfont=dict(color=t["ink_2"], size=13))
    return fig


def cell_counts_figure(counts: dict, column: str, theme: str = "light",
                       max_categories: int = 25) -> go.Figure:
    """Horizontal bars of cells per category, largest first, each labelled."""
    return ranked_bars_figure(counts, f"Cells per {column}",
                              f"{len(counts):,} categories · {sum(counts.values()):,} cells",
                              "cells", theme, max_categories)


def histogram_figure(start: list, end: list, count: list, title: str, subtitle: str,
                     x_title: str, unit: str, theme: str = "light",
                     marker_x: float | None = None, marker_label: str = "",
                     log_y: bool = False, plot_px: int = 560,
                     bin_width: float | None = None) -> go.Figure:
    """Columns for a binned distribution, with an optional labelled vertical marker.

    ``bin_width`` sets the width of one bin; by default it is taken from the
    first bin. Pass 1 for whole-number values listed without their empty bins.
    """
    t = THEMES[theme]
    centers = [(a + b) / 2 for a, b in zip(start, end)]
    if bin_width is None:
        bin_width = (start[1] - start[0]) if len(start) > 1 else ((end[0] - start[0]) or 1)
    width = bin_width
    span = (end[-1] - start[0]) + width
    slot_px = plot_px * width / span
    total = sum(count) or 1
    labels = [f"{a:,.4g}" if a == b else f"{a:,.4g}–{b:,.4g}" for a, b in zip(start, end)]

    fig = go.Figure(go.Bar(
        x=centers, y=count, width=width * min(0.85, max(0.4, BAR_PX / max(slot_px, 1))),
        marker=dict(color=t["series"][0], cornerradius=3 if slot_px > 8 else 0),
        customdata=list(zip(labels, [100 * n / total for n in count])),
        hovertemplate="<b>%{customdata[0]}</b><br>%{y:,} " + unit
                      + " (%{customdata[1]:.1f}%)<extra></extra>",
    ))
    if marker_x is not None:
        fig.add_vline(x=marker_x, line=dict(color=t["ink_2"], width=1))
        fig.add_annotation(x=marker_x, y=1, yref="paper", text=marker_label, showarrow=False,
                           xanchor="left", yanchor="top", xshift=4,
                           font=dict(size=12, color=t["ink"]))
    _base(fig, t, title, subtitle, height=330)
    fig.update_xaxes(title_text=x_title, showgrid=False, showline=True,
                     range=[start[0] - width, end[-1] + width])
    if log_y:
        fig.update_yaxes(title_text=f"{unit.capitalize()} (log scale)", type="log")
    else:
        fig.update_yaxes(title_text=unit.capitalize(), rangemode="tozero")
    return fig


def _weighted_median(values: list, weights: list) -> float:
    half, running = sum(weights) / 2, 0
    for v, w in sorted(zip(values, weights)):
        running += w
        if running >= half:
            return v
    return values[-1] if values else 0


# VCF

# Collapse the 12 substitutions to the 6 pyrimidine-reference classes:
# a change on one strand is the same event as its complement on the other
_SPECTRUM = {
    "C>A": ("C>A", "G>T"), "C>G": ("C>G", "G>C"), "C>T": ("C>T", "G>A"),
    "T>A": ("T>A", "A>T"), "T>C": ("T>C", "A>G"), "T>G": ("T>G", "A>C"),
}
_TRANSITIONS = {"C>T", "T>C"}


def substitution_figure(stats: dict, theme: str = "light") -> go.Figure:
    """Six-class substitution spectrum, transitions and transversions colored apart."""
    t = THEMES[theme]
    subs = stats["substitutions"]
    classes = list(_SPECTRUM)
    counts = [sum(subs.get(s, 0) for s in _SPECTRUM[c]) for c in classes]
    total = sum(counts) or 1

    fig = go.Figure()
    for kind, color in [("Transition", t["series"][0]), ("Transversion", t["series"][1])]:
        picked = [(c, n) for c, n in zip(classes, counts) if (c in _TRANSITIONS) == (kind == "Transition")]
        fig.add_trace(go.Bar(
            x=[c for c, _ in picked], y=[n for _, n in picked], name=kind,
            marker=dict(color=color, cornerradius=4),
            customdata=[100 * n / total for _, n in picked],
            hovertemplate="<b>%{x}</b><br>%{y:,} SNVs (%{customdata:.1f}%)<extra>" + kind + "</extra>",
        ))
    ts_tv = stats.get("ts_tv")
    _base(fig, t, "Substitution spectrum",
          f"{sum(counts):,} single-base changes, strands combined"
          + (f" · Ts/Tv {ts_tv:.2f}" if ts_tv is not None else ""), height=360)
    fig.update_layout(showlegend=True, bargap=_bar_gap(6, 420 / 6), margin=dict(t=108))
    fig.update_xaxes(categoryorder="array", categoryarray=classes, showgrid=False,
                     showline=True, tickfont=dict(color=t["ink_2"], size=13))
    fig.update_yaxes(title_text="SNVs")
    _count_axis(fig, max(counts))
    return fig


def variant_types_figure(stats: dict, theme: str = "light") -> go.Figure:
    """Records by variant type as labelled horizontal bars."""
    labels = {"SNPs": "SNPs", "indels": "Indels", "MNPs": "MNPs", "others": "Other"}
    counts = {name: stats.get(key, 0) for key, name in labels.items()}
    fig = cell_counts_figure(counts, "variant type", theme)
    fig.update_layout(title_text=(
        f"<b>Variants by type</b><br><span style='font-size:13px;color:{THEMES[theme]['ink_2']}'>"
        f"{stats['records']:,} records · {stats['samples']} samples</span>"
    ))
    fig.update_traces(hovertemplate="<b>%{y}</b><br>%{x:,} records<extra></extra>")
    return fig


# BAM quality report (bam_qc)

def mapq_figure(qc: dict, theme: str = "light") -> go.Figure:
    """Reads per mapping quality, log scale, with MAPQ 30 marked."""
    mapq, reads = qc["mapq"]["mapq"], qc["mapq"]["reads"]
    total = sum(reads) or 1
    high = 100 * sum(n for q, n in zip(mapq, reads) if q >= 30) / total
    return histogram_figure(
        mapq, mapq, reads, "Mapping quality",
        f"{high:.1f}% of mapped reads at MAPQ 30 or higher · MAPQ 0 = could map equally well elsewhere",
        "MAPQ", "reads", theme, marker_x=30, marker_label="MAPQ 30", log_y=True, bin_width=1,
    )


def insert_size_figure(qc: dict, theme: str = "light") -> go.Figure | None:
    """Distribution of insert sizes for read pairs, with the median marked."""
    hist = qc["insert_size"]
    if not hist["count"]:
        return None
    t = THEMES[theme]
    blue = t["series"][0]
    x = [(a + b) / 2 for a, b in zip(hist["start"], hist["end"])]
    median = _weighted_median(x, hist["count"])
    summary = qc["summary"]

    fig = go.Figure(go.Scatter(
        x=x, y=hist["count"], mode="lines", line=dict(color=blue, width=2),
        fill="tozeroy", fillcolor=_rgba(blue, 0.10),
        customdata=list(zip(hist["start"], hist["end"])),
        hovertemplate="<b>%{customdata[0]:,}–%{customdata[1]:,} bp</b><br>%{y:,} pairs<extra></extra>",
    ))
    fig.add_vline(x=median, line=dict(color=t["ink_2"], width=1))
    fig.add_annotation(x=median, y=1, yref="paper", text=f"Median {median:,.0f} bp",
                       showarrow=False, xanchor="left", yanchor="top", xshift=4,
                       font=dict(size=12, color=t["ink"]))
    _base(fig, t, "Insert size",
          f"Mean {summary['insert_size_mean']:,.1f} bp ± {summary['insert_size_sd']:,.1f} · "
          f"{summary['inward_pairs']:,} inward and {summary['outward_pairs']:,} outward pairs",
          height=330)
    fig.update_xaxes(title_text="Insert size (bp)", showgrid=False, showline=True, tickformat=",d")
    fig.update_yaxes(title_text="Read pairs", rangemode="tozero")
    fig.update_layout(hovermode="x")
    return fig


def coverage_distribution_figure(qc: dict, theme: str = "light") -> go.Figure | None:
    """Share of covered bases reaching each depth, with common thresholds labelled."""
    depth, bases = qc["coverage"]["depth"], qc["coverage"]["bases"]
    if not depth:
        return None
    t = THEMES[theme]
    blue = t["series"][0]
    total = sum(bases)
    # Share of covered bases at or above each depth (reverse cumulative sum)
    at_least, running = [], total
    for b in bases:
        at_least.append(100 * running / total)
        running -= b

    fig = go.Figure(go.Scatter(
        x=depth, y=at_least, mode="lines", line=dict(color=blue, width=2, shape="hv"),
        fill="tozeroy", fillcolor=_rgba(blue, 0.10),
        hovertemplate="<b>%{y:.1f}%</b> of covered bases at ≥ %{x:,}×<extra></extra>",
    ))
    for threshold in (10, 20, 30, 100):
        share = qc["coverage"]["pct_bases_at_least"].get(f"{threshold}x")
        if share is None or threshold > max(depth):
            continue
        fig.add_trace(go.Scatter(
            x=[threshold], y=[share], mode="markers", hoverinfo="skip",
            marker=dict(size=9, color=blue, line=dict(color=t["surface"], width=2)),
        ))
        # On a log axis, annotation positions are given as log10 of the value
        fig.add_annotation(x=math.log10(threshold), y=share, text=f"{share:.0f}% ≥ {threshold}×",
                           showarrow=False, xanchor="left", yanchor="bottom", xshift=6, yshift=4,
                           font=dict(size=12, color=t["ink"]))
    _base(fig, t, "Coverage depth",
          "Share of covered bases at or above each depth · bases with no reads are not counted "
          "· duplicates included",
          height=350)
    fig.update_xaxes(title_text="Depth (×, log scale)", type="log", showgrid=False, showline=True)
    fig.update_yaxes(title_text="Covered bases (%)", range=[0, 105], ticksuffix="%")
    fig.update_layout(hovermode="x")
    return fig


def contigs_figure(qc: dict, theme: str = "light") -> go.Figure | None:
    """Mapped reads per contig, largest first."""
    contigs = [c for c in qc["contigs"] if c["mapped"]]
    if not contigs:
        return None
    counts = {c["contig"]: c["mapped"] for c in contigs}
    subtitle = f"{len(contigs):,} contigs with reads"
    if qc.get("unplaced_unmapped"):
        subtitle += f" · {qc['unplaced_unmapped']:,} unmapped reads without a position"
    return ranked_bars_figure(counts, "Mapped reads per contig", subtitle, "reads", theme)


def bam_read_length_figure(qc: dict, theme: str = "light") -> go.Figure | None:
    hist = qc["read_length"]
    if not hist["count"]:
        return None
    return histogram_figure(hist["start"], hist["end"], hist["count"], "Read length",
                            f"Average {qc['summary']['avg_length']:,.0f} bp · "
                            f"longest {qc['summary']['max_length']:,} bp",
                            "Read length (bp)", "reads", theme)


# FASTQ composition and per-read distributions (fastq_quality_profile)

def base_composition_figure(profile: dict, theme: str = "light") -> go.Figure:
    """Percentage of A, C, G, T (and N) at each read position, lines labelled at the end."""
    t = THEMES[theme]
    x = [(s + e) / 2 for s, e in zip(profile["position_start"], profile["position_end"])]
    fig = go.Figure()
    colors = dict(zip("ACGT", t["series"]))
    lines = [(b, colors[b]) for b in "ACGT"]
    if max(profile["base_pct"]["N"], default=0) > 0:
        lines.append(("N", t["muted"]))
    for base, color in lines:
        y = profile["base_pct"][base]
        fig.add_trace(go.Scatter(
            x=x, y=y, mode="lines", name=base, line=dict(color=color, width=2),
            hovertemplate=f"<b>{base}</b> %{{y:.1f}}% at %{{x:,.0f}}<extra></extra>",
        ))
    # Name lines at their right end where the labels have room; the legend
    # covers the rest, so identity never rests on color alone
    placed = []
    for base, _ in sorted(lines, key=lambda bl: profile["base_pct"][bl[0]][-1]):
        end_y = profile["base_pct"][base][-1]
        if all(abs(end_y - other) >= 6 for other in placed):
            placed.append(end_y)
            fig.add_annotation(x=x[-1], y=end_y, text=base, showarrow=False, xanchor="left",
                               xshift=6, font=dict(size=12, color=t["ink"]))
    _base(fig, t, "Base composition by read position",
          "Whole-genome reads give flat lines near 25%. Amplicons such as 16S swing "
          "position by position, since every read covers the same gene.",
          height=360)
    fig.update_xaxes(title_text="Position in read (bp)", showgrid=False, showline=True)
    fig.update_yaxes(title_text="Bases (%)", range=[0, 100], ticksuffix="%")
    fig.update_layout(showlegend=True, hovermode="x unified", margin=dict(t=104))
    return fig


def read_quality_figure(profile: dict, theme: str = "light") -> go.Figure:
    """Reads by their mean quality, with Q20 marked and low-quality reads counted."""
    q, reads = profile["read_quality"]["quality"], profile["read_quality"]["reads"]
    total = sum(reads) or 1
    low = sum(n for v, n in zip(q, reads) if v < 20)
    return histogram_figure(
        q, q, reads, "Mean quality per read",
        f"{low:,} reads ({100 * low / total:.1f}%) average below Q20",
        "Mean Phred quality of the read", "reads", theme, marker_x=20, marker_label="Q20",
        bin_width=1,
    )


def read_gc_figure(profile: dict, theme: str = "light") -> go.Figure:
    """Reads by their GC content, with the median marked."""
    gc, reads = profile["read_gc"]["gc_pct"], profile["read_gc"]["reads"]
    median = _weighted_median(gc, reads)
    return histogram_figure(
        gc, gc, reads, "GC content per read",
        "A single smooth peak is expected; a second peak can mean contamination",
        "GC (%)", "reads", theme, marker_x=median, marker_label=f"Median {median:.0f}%",
        bin_width=1,
    )


# VCF detail (vcf_stats)

def qual_figure(stats: dict, theme: str = "light") -> go.Figure | None:
    hist = stats["qual"]
    if not hist["count"]:
        return None
    return histogram_figure(hist["start"], hist["end"], hist["count"], "Variant quality (QUAL)",
                            f"{sum(hist['count']):,} sites with a QUAL value",
                            "QUAL", "sites", theme)


def indel_length_figure(stats: dict, theme: str = "light") -> go.Figure | None:
    """Indel sites by length: deletions left of zero, insertions right."""
    lengths, sites = stats["indel_lengths"]["length"], stats["indel_lengths"]["sites"]
    if not lengths:
        return None
    t = THEMES[theme]
    fig = go.Figure()
    for name, color, keep in [("Deletion", t["series"][1], lambda L: L < 0),
                              ("Insertion", t["series"][0], lambda L: L > 0)]:
        picked = [(L, n) for L, n in zip(lengths, sites) if keep(L)]
        if picked:
            fig.add_trace(go.Bar(
                x=[L for L, _ in picked], y=[n for _, n in picked], name=name,
                marker=dict(color=color, cornerradius=3), width=0.7,
                hovertemplate="<b>%{x:+d} bp</b><br>%{y:,} sites<extra>" + name + "</extra>",
            ))
    deletions = sum(n for L, n in zip(lengths, sites) if L < 0)
    insertions = sum(n for L, n in zip(lengths, sites) if L > 0)
    _base(fig, t, "Indel lengths",
          f"{deletions:,} deletions · {insertions:,} insertions · "
          "in coding regions, lengths that are multiples of 3 keep the reading frame",
          height=340)
    fig.update_xaxes(title_text="Length (bp; deletions negative)", showgrid=False,
                     showline=True, zeroline=True, zerolinecolor=t["axis"])
    fig.update_yaxes(title_text="Sites")
    _count_axis(fig, max(sites))
    fig.update_layout(showlegend=True, margin=dict(t=104))
    return fig


def genotype_depth_figure(stats: dict, theme: str = "light") -> go.Figure | None:
    labels, genotypes = stats["genotype_depth"]["depth"], stats["genotype_depth"]["genotypes"]
    if not labels:
        return None
    t = THEMES[theme]
    fig = go.Figure(go.Bar(
        x=labels, y=genotypes, marker=dict(color=t["series"][0], cornerradius=3),
        hovertemplate="<b>Depth %{x}</b><br>%{y:,} genotypes<extra></extra>",
    ))
    _base(fig, t, "Genotype depth (FORMAT/DP)",
          f"{sum(genotypes):,} genotypes · bcftools groups depths above 500 together", height=330)
    fig.update_layout(bargap=_bar_gap(len(labels), 560 / max(len(labels), 1)))
    fig.update_xaxes(title_text="Read depth", type="category", showgrid=False, showline=True)
    fig.update_yaxes(title_text="Genotypes")
    _count_axis(fig, max(genotypes))
    return fig


def per_sample_figure(stats: dict, theme: str = "light", max_samples: int = 30) -> go.Figure | None:
    """Genotype calls per sample as stacked bars: hom-ref, het, hom-alt, missing."""
    samples = stats["per_sample"][:max_samples]
    if not samples:
        return None
    t = THEMES[theme]
    names = [s["sample"] for s in samples][::-1]
    parts = [("Hom-ref", "hom_ref", t["series"][0]), ("Het", "het", t["series"][1]),
             ("Hom-alt", "hom_alt", t["series"][2]), ("Missing", "missing", t["context"])]
    fig = go.Figure()
    for label, key, color in parts:
        values = [s[key] for s in samples][::-1]
        if any(values):
            fig.add_trace(go.Bar(
                y=names, x=values, name=label, orientation="h",
                marker=dict(color=color, line=dict(color=t["surface"], width=2)),
                hovertemplate="<b>%{y}</b><br>" + label + ": %{x:,}<extra></extra>",
            ))
    ratios = [f"{s['sample']} het/hom-alt {s['het'] / s['hom_alt']:.2f}"
              for s in samples[:3] if s["hom_alt"]]
    _base(fig, t, "Genotypes per sample", " · ".join(ratios) or f"{len(samples)} samples",
          height=150 + 34 * len(names))
    fig.update_layout(barmode="stack", bargap=_bar_gap(len(names), 34), showlegend=True,
                      margin=dict(t=104))
    fig.update_xaxes(title_text="Variant sites", showline=True)
    fig.update_yaxes(showgrid=False, tickfont=dict(color=t["ink_2"]))
    return fig


# h5ad QC metrics and embeddings

def obs_histogram_figure(col: dict, theme: str = "light") -> go.Figure | None:
    """Histogram of a numeric obs column with its median marked."""
    if col.get("kind") != "numeric" or "histogram" not in col:
        return None
    h = col["histogram"]
    return histogram_figure(
        h["start"], h["end"], h["count"], col["column"],
        f"Median {col['median']:,.4g} · middle 50% {col['q1']:,.4g}–{col['q3']:,.4g} · "
        f"range {col['min']:,.4g}–{col['max']:,.4g}",
        col["column"], "cells", theme, marker_x=col["median"], marker_label="Median",
    )


def embedding_figure(emb: dict, theme: str = "light", highlight: str | None = None) -> go.Figure:
    """2-D embedding. Categories: one highlighted, the rest gray, every group labelled.
    Numbers: a single-hue color scale. No color column: plain points."""
    t = THEMES[theme]
    x, y = emb["x"], emb["y"]
    values = emb.get("values")
    point = dict(size=5, line=dict(width=0))
    fig = go.Figure()
    name = emb["basis"].removeprefix("X_").upper()

    if values is None:
        fig.add_trace(go.Scattergl(x=x, y=y, mode="markers", hoverinfo="skip",
                                   marker=dict(color=t["series"][0], opacity=0.6, **point)))
        subtitle = f"{emb['n_cells']:,} cells"
    elif isinstance(values[0], (int, float)):
        fig.add_trace(go.Scattergl(
            x=x, y=y, mode="markers",
            marker=dict(color=values, colorscale=t["sequential"], opacity=0.85,
                        colorbar=dict(title=dict(text=emb["color"], side="right"),
                                      thickness=10, outlinewidth=0,
                                      tickfont=dict(color=t["muted"])), **point),
            hovertemplate=emb["color"] + " %{marker.color:,.4g}<extra></extra>",
        ))
        subtitle = f"Colored by {emb['color']} · {emb['n_cells']:,} cells"
    else:
        groups = {}
        for xi, yi, v in zip(x, y, values):
            groups.setdefault(v, []).append((xi, yi))
        if highlight not in groups:
            highlight = max(groups, key=lambda g: len(groups[g]))
        rest = [(p, g) for g, pts in groups.items() if g != highlight for p in pts]
        mine = groups[highlight]
        fig.add_trace(go.Scattergl(
            x=[p[0] for p, _ in rest], y=[p[1] for p, _ in rest], mode="markers",
            name="Other cells", text=[g for _, g in rest], hovertemplate="%{text}<extra></extra>",
            marker=dict(color=t["context"], opacity=0.5, **point)))
        fig.add_trace(go.Scattergl(
            x=[p[0] for p in mine], y=[p[1] for p in mine], mode="markers",
            name=str(highlight), hovertemplate=f"<b>{highlight}</b><extra></extra>",
            marker=dict(color=t["series"][0], opacity=0.9, **point)))
        # Label groups at their median position; very small groups are left to
        # the hover text (unless highlighted) so labels don't pile up
        for g, pts in groups.items():
            if g != highlight and len(pts) < 0.02 * len(values):
                continue
            gx = sorted(p[0] for p in pts)[len(pts) // 2]
            gy = sorted(p[1] for p in pts)[len(pts) // 2]
            fig.add_annotation(x=gx, y=gy, text=f"<b>{g}</b>" if g == highlight else str(g),
                               showarrow=False, bgcolor=_rgba(t["surface"], 0.75), borderpad=2,
                               font=dict(size=12, color=t["ink"] if g == highlight else t["ink_2"]))
        share = 100 * len(mine) / len(values)
        subtitle = (f"{highlight}: {len(mine):,} cells ({share:.1f}%) highlighted · "
                    f"{len(groups)} groups in {emb['color']}")
        fig.update_layout(showlegend=True)

    if emb.get("sampled"):
        subtitle += f" · random sample of {len(x):,} shown"
    _base(fig, t, name, subtitle, height=560)
    fig.update_xaxes(showticklabels=False, showgrid=False, zeroline=False, title_text=f"{name} 1")
    fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False, title_text=f"{name} 2",
                     scaleanchor="x", scaleratio=1)
    fig.update_layout(margin=dict(t=104 if values is not None else 80))
    return fig
