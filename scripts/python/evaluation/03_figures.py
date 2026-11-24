"""Evaluation step 3 - the figures, drawn from the tables of step 2 (data/evaluation/).

Run (repo root, bioenv active):   python scripts/python/evaluation/03_figures.py
Reads:   data/evaluation/  summary.tsv, scores.tsv, calls.tsv, profile_bins.tsv.gz, profile_snps.tsv.gz
         and, for labels only, the metadata, arrays and SNP bands named in inputs.py
Writes:  results/figures/evaluation/*.png  (FIGURES below lists them)

"""
import importlib
import re
import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))
import inputs as I  # noqa: E402
ev = importlib.import_module("02_evaluate")      # its parent-karyotype reader, so bands are placed the same way

# what to draw
OUT = REPO / "results/figures/evaluation"
MAIN = "Pipeline end to end"                     # the setup shown alone in the single-setup figures
EXAMPLES = {"SNP031CXJ2": [9, 12],               # embryo -> chromosomes for the close-up (empty: none)
            "SNP011ZJ3": [15],
            "SNP013TYY5": []}
DPI = 200
SIZE = (13.33, 7.5)                              # inches: 16:9, the shape of a slide

FIGURES = """
  01_profile_<embryo>.png         genome-wide LRR (both inputs) and the calls of every setup, for each EXAMPLES embryo
  02_close_up_<embryo>.png        the EXAMPLES chromosomes: SNPs, segment levels, authors' region, parent's band
  03_cnv_map_<setup>.png          every call of one setup against the authors' regions, one figure per setup
  04_benchmark.png                MAIN against the authors: embryo and event scores with 95% intervals
  05_setup_comparison.png         the same scores for all setups, plus the parent checks
  06_parents.png                  MAIN's calls on the parents' translocation chromosomes, with the breakpoint bands
  07_errors_vs_noise.png          missed and extra calls per embryo against that embryo's noise, one panel per setup
  08_baf_by_copy_number.png       BAF of known 1-, 2- and 3-copy regions, for each input
  09_x_dosage.png                 X against autosome LRR in XX and XY embryos: the one-copy step of each input
"""

# look
C = dict(loss="#1A4F94", mosaic_loss="#86B6EE", gain="#D03B3B", mosaic_gain="#F2A582", gain2="#6E1010",
         chrom="#EFEEEA", snp="#B9B7B0", bin="#2B2B2B", band="#F3D36B", text="#222222", muted="#666666")
CALLER = {"ASPCF + rules (06c)": "#4A3AA8", "ASCAT": "#D9822B", "ASPCF + GBoost": "#1A9E77"}
INPUT_MARKER = {"pipeline": "o", "CRMA v2": "s"}
plt.rcParams.update({"font.size": 11, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": "#999999", "axes.labelcolor": C["muted"],
                     "xtick.color": C["muted"], "ytick.color": C["muted"], "legend.frameon": False})


def at(p) -> Path:
    return REPO / p


def slug(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def figure(title: str, subtitle: str, note: str = "", **grid):
    """A 16:9 figure with the title, subtitle and footnote in the same places every time."""
    fig = plt.figure(figsize=SIZE)
    fig.text(0.015, 0.965, title, fontsize=20, weight="bold", color=C["text"], va="top")
    fig.text(0.015, 0.905, textwrap.fill(subtitle, 135), fontsize=12.5, color=C["muted"], va="top")
    if note:
        fig.text(0.015, 0.02, textwrap.fill(note, 170), fontsize=9.5, color=C["muted"], va="bottom")
    grid = dict(left=0.17, right=0.9, top=0.8, bottom=0.12) | grid
    return fig, grid


def save(fig, name: str) -> None:
    fig.savefig(OUT / name, dpi=DPI)
    plt.close(fig)
    print("  ", name)


def pct(v) -> str:
    return f"{100 * v:.0f}%"


# the tables of step 2
T = at(I.TABLE_FOLDER)
summary = pd.read_csv(T / "summary.tsv", sep="\t")
scores = pd.read_csv(T / "scores.tsv", sep="\t")
calls = pd.read_csv(T / "calls.tsv", sep="\t")
bins = pd.read_csv(T / "profile_bins.tsv.gz", sep="\t")
meta = pd.read_csv(at(I.METADATA), sep="\t", dtype={"couple": str}).set_index("sample")
n_x = pd.read_csv(at(I.ARRAYS), sep="\t").set_index("sample")["n_X"]

SETUPS = [s for s in summary["setup"] if s != "Authors"]
INPUT = dict(zip(summary["setup"], summary["input"]))
CALLER_OF = dict(zip(summary["setup"], summary["caller"]))
authors = calls[calls["setup"] == "Authors"]

# chromosomes laid end to end, as genome_mb in the tables (1 first, X last)
ends = bins[bins["input"] == bins["input"].iloc[0]].groupby("chr").agg(first=("start", "min"), last=("end", "max"))
ends["offset"] = (ends["last"] / 1e6).cumsum().shift(fill_value=0.0)
ends["g0"] = ends["offset"] + ends["first"] / 1e6
ends["g1"] = ends["offset"] + ends["last"] / 1e6
GENOME_MB = ends["g1"].iloc[-1]


def chr_name(c) -> str:
    return "X" if int(c) == 23 else str(int(c))


def call_colour(direction: int, label) -> str:
    label = str(label)
    copies = [int(n) for n in re.findall(r"(\d+) copies", label)]
    if direction > 0:
        if copies and max(copies) >= 4:
            return C["gain2"]
        return C["mosaic_gain"] if "mosaic" in label else C["gain"]
    return C["mosaic_loss"] if "mosaic" in label else C["loss"]


def found(a: pd.DataFrame, b: pd.DataFrame) -> np.ndarray:
    """For every row of a: does a row of b overlap it (same embryo, chromosome and direction)?"""
    m = a.reset_index().merge(b[["sample", "chr", "direction", "start", "end"]],
                              on=["sample", "chr", "direction"], suffixes=("", "_b"))
    hit = m.loc[(m["start"] < m["end_b"]) & (m["start_b"] < m["end"]), "index"].unique()
    return a.index.isin(hit)


def errors(setup: str) -> tuple:
    """(authors' regions with found / missed, the setup's calls with extra = no authors' region)."""
    c = calls[calls["setup"] == setup]
    return authors.assign(found=found(authors, c)), c.assign(extra=c["overlap"].fillna(0) == 0)


def chromosome_axis(ax, two_rows: bool = True, top: bool = False) -> None:
    mids = (ends["g0"] + ends["g1"]) / 2
    ax.set_xlim(0, GENOME_MB)
    ax.set_xticks(mids)
    second = (lambda n: n + "\n") if top else (lambda n: "\n" + n)
    labels = [chr_name(c) if (not two_rows or i % 2 == 0) else second(chr_name(c)) for i, c in enumerate(ends.index)]
    ax.set_xticklabels(labels, fontsize=9)
    ax.tick_params(axis="x", length=0, labeltop=top, labelbottom=not top)


def shade_chromosomes(ax) -> None:
    for i, r in enumerate(ends.itertuples()):
        if i % 2:
            ax.axvspan(r.g0, r.g1, color="#F6F5F2", zorder=0, lw=0)


def call_legend(fig, y: float = 0.055, extra: tuple = ()) -> None:
    handles = [Patch(color=C["loss"], label="Loss"), Patch(color=C["mosaic_loss"], label="Mosaic loss"),
               Patch(color=C["mosaic_gain"], label="Mosaic gain"), Patch(color=C["gain"], label="Gain"),
               Patch(color=C["gain2"], label="Gain of 2+ copies")] + list(extra)
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), bbox_to_anchor=(0.5, y), fontsize=11)


# 01 genome-wide profile of one embryo

def profile(sample: str, snps: pd.DataFrame) -> None:
    inputs = [i for i in INPUT_MARKER if i in set(snps["input"])]
    sources = ["Authors"] + SETUPS
    m = calls[calls["setup"] == MAIN]
    called = ", ".join(f"{'+' if r.direction > 0 else '−'}{chr_name(r.chr)} ({r.size_mb:.0f} Mb)"
                       for r in m[m["sample"] == sample].itertuples()) or "no change"
    fig, g = figure(f"Genome-wide LRR and calls: {sample}",
                    f"Parent {meta.loc[sample, 'donor_karyotype']}. Authors: {meta.loc[sample, 'reported_result']}. "
                    f"{MAIN}: {called}.",
                    "Grey = every 10th SNP, black = median of 50 SNPs. Bars = calls of each setup; the first row is the authors'.",
                    left=0.19, bottom=0.11)
    gs = fig.add_gridspec(len(inputs) + 1, 1, height_ratios=[3] * len(inputs) + [0.42 * len(sources)], hspace=0.12, **g)
    for k, inp in enumerate(inputs):
        ax = fig.add_subplot(gs[k])
        d = snps[(snps["input"] == inp) & (snps["sample"] == sample)]
        b = bins[(bins["input"] == inp) & (bins["sample"] == sample)]
        shade_chromosomes(ax)
        ax.scatter(d["genome_mb"], d["lrr"], s=0.6, color=C["snp"], alpha=0.5, lw=0, rasterized=True)
        ax.scatter(b["genome_mb"], b["lrr"], s=2.5, color=C["bin"], lw=0, rasterized=True)
        ax.axhline(0, color=C["muted"], lw=0.6)
        ax.set_ylim(-1.2, 1.2)
        ax.set_ylabel(f"LRR\n{inp}")
        chromosome_axis(ax)
        ax.tick_params(labelbottom=False)
    ax = fig.add_subplot(gs[-1])
    shade_chromosomes(ax)
    for i, src in enumerate(sources):
        for r in calls[(calls["setup"] == src) & (calls["sample"] == sample)].itertuples():
            w = max(r.genome_end_mb - r.genome_start_mb, 4)
            ax.add_patch(Rectangle((r.genome_start_mb, i - 0.35), w, 0.7, color=call_colour(r.direction, r.label)))
    ax.set_ylim(len(sources) - 0.5, -0.5)
    ax.set_yticks(range(len(sources)))
    ax.set_yticklabels(sources, fontsize=9.5)
    ax.tick_params(axis="y", length=0)
    ax.spines[["left", "bottom"]].set_visible(False)
    chromosome_axis(ax)
    save(fig, f"01_profile_{sample}.png")


# 02 close-up of chromosomes of one embryo

def close_up(sample: str, chromosomes: list, snps: pd.DataFrame, parents: pd.DataFrame) -> None:
    inp = INPUT[MAIN]
    m = calls[(calls["setup"] == MAIN) & (calls["sample"] == sample)]
    a = authors[authors["sample"] == sample]
    parts = []
    for c in chromosomes:
        mc, ac = m[m["chr"] == c], a[a["chr"] == c]
        p = parents[(parents["sample"] == sample) & (parents["chr"] == c)]
        bit = f"chr {chr_name(c)}: authors " + (", ".join(f"{r.start / 1e6:.1f}–{r.end / 1e6:.1f} Mb" for r in ac.itertuples()) or "none")
        bit += ", called " + (", ".join(f"{r.start / 1e6:.1f}–{r.end / 1e6:.1f} Mb" for r in mc.itertuples()) or "none")
        if len(p):
            bit += f", parent's band {chr_name(c)}{p['band'].iloc[0]} at {p['band_start'].iloc[0] / 1e6:.1f}–{p['band_end'].iloc[0] / 1e6:.1f} Mb"
        parts.append(bit)
    fig, g = figure(f"Segmentation close-up: {sample}", "; ".join(parts) + ".",
                    f"{inp} LRR. Grey = every 10th SNP, black = median of 50 SNPs, coloured line = median LRR of the call "
                    f"({MAIN}). Yellow = the parent's breakpoint band. Outline above = the authors' region.",
                    left=0.07, right=0.98, bottom=0.15, wspace=0.12)
    gs = fig.add_gridspec(1, len(chromosomes), **g)
    for k, c in enumerate(chromosomes):
        ax = fig.add_subplot(gs[k])
        d = snps[(snps["input"] == inp) & (snps["sample"] == sample) & (snps["chr"] == c)]
        b = bins[(bins["input"] == inp) & (bins["sample"] == sample) & (bins["chr"] == c)]
        x = (b["start"] + b["end"]) / 2e6
        p = parents[(parents["sample"] == sample) & (parents["chr"] == c)]
        for r in p.itertuples():
            lo, hi = r.band_start / 1e6, max(r.band_end / 1e6, r.band_start / 1e6 + 0.8)
            ax.axvspan(lo, hi, color=C["band"], alpha=0.6, lw=0, zorder=0)
        ax.scatter(d["pos"] / 1e6, d["lrr"], s=2, color=C["snp"], alpha=0.6, lw=0, rasterized=True)
        ax.scatter(x, b["lrr"], s=7, color=C["bin"], lw=0)
        inside = np.zeros(len(b), bool)
        for r in m[m["chr"] == c].itertuples():
            here = ((b["start"] >= r.start) & (b["end"] <= r.end)).to_numpy()
            inside |= here
            ax.plot([r.start / 1e6, r.end / 1e6], [b["lrr"][here].median()] * 2, color=call_colour(r.direction, r.label), lw=4)
        if (~inside).any():
            ax.plot([b["start"].min() / 1e6, b["end"].max() / 1e6], [b["lrr"][~inside].median()] * 2,
                    color=C["muted"], lw=1.5, ls="--")
        for r in a[a["chr"] == c].itertuples():
            ax.add_patch(Rectangle((r.start / 1e6, 0.92), (r.end - r.start) / 1e6, 0.12, fill=False, lw=1.3,
                                   edgecolor=C["text"]))
            ax.text(r.start / 1e6, 1.08, f"authors: {r.label}", fontsize=9.5, color=C["text"], va="bottom")
        ax.axhline(0, color=C["muted"], lw=0.6)
        ax.set_ylim(-1.2, 1.25)
        ax.set_xlabel(f"Chromosome {chr_name(c)} (Mb, hg19)")
        if k == 0:
            ax.set_ylabel(f"LRR ({inp})")
    save(fig, f"02_close_up_{sample}.png")


# 03 every call of one setup in the cohort

def cnv_map(setup: str, parents: pd.DataFrame) -> None:
    a, c = errors(setup)
    rows = [s for s in meta.index if s in set(a["sample"])]
    parent = parents.groupby("sample")["parent"].first()
    rows = sorted(rows, key=lambda s: (parent[s].startswith("rob"), -len(a[a["sample"] == s]), s))
    only_setup = sorted(set(c["sample"]) - set(rows))
    rows += only_setup
    n_found, n_extra = int(a["found"].sum()), int(c["extra"].sum())
    fig, g = figure(f"Copy-number calls in the cohort: {setup}",
                    f"{n_found} of the {len(a)} reported changes are found. "
                    + (f"{n_extra} of {len(c)} calls are not reported by the authors"
                       + (f", {len(only_setup)} of the embryos they call normal get a call." if only_setup else ".")
                       if n_extra else "There is no call that the authors do not report."),
                    f"Fill = call, outline = the authors' region, × = missed. Calls under 3 Mb are drawn 3 Mb wide. "
                    f"Not shown: {len(meta) - len(rows)} embryos normal in both.",
                    left=0.2, right=0.9, top=0.78, bottom=0.14)
    ax = fig.add_axes([g["left"], g["bottom"], g["right"] - g["left"], g["top"] - g["bottom"]])
    for i, s in enumerate(rows):
        for r in ends.itertuples():
            ax.add_patch(Rectangle((r.g0 + 1.5, i - 0.32), r.g1 - r.g0 - 3, 0.64, color=C["chrom"], lw=0))
        for r in c[c["sample"] == s].itertuples():
            w = max(r.genome_end_mb - r.genome_start_mb, 3)
            ax.add_patch(Rectangle((r.genome_start_mb, i - 0.32), w, 0.64, color=call_colour(r.direction, r.label), lw=0))
        for r in a[a["sample"] == s].itertuples():
            w = max(r.genome_end_mb - r.genome_start_mb, 3)
            ax.add_patch(Rectangle((r.genome_start_mb, i - 0.42), w, 0.84, fill=False, lw=1.2, edgecolor=C["text"]))
            if not r.found:
                ax.plot(r.genome_start_mb + w / 2, i, marker="x", color=C["text"], ms=6, mew=1.3)
    if only_setup:
        ax.axhline(len(rows) - len(only_setup) - 0.5, color=C["muted"], lw=0.8, ls=":")
        ax.text(GENOME_MB + 8, len(rows) - len(only_setup) / 2 - 0.5, "Normal\nin GEO", fontsize=10, weight="bold", va="center")
    ax.text(GENOME_MB + 8, (len(rows) - len(only_setup)) / 2 - 0.5, "Abnormal\nin GEO", fontsize=10, weight="bold", va="center")
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{s}  {parent[s]}" for s in rows], fontsize=8.5 if len(rows) > 20 else 10)
    ax.tick_params(length=0)
    ax.spines[["left", "bottom"]].set_visible(False)
    chromosome_axis(ax, top=True)
    call_legend(fig, extra=(Patch(facecolor="none", edgecolor=C["text"], label="Reported by authors"),
                            Line2D([], [], marker="x", color=C["text"], ls="", label="Missed")))
    save(fig, f"03_cnv_map_{slug(setup)}.png")


# 04 / 05 scores against the authors

METRICS = [("embryo", "sensitivity", "Embryo sensitivity"), ("embryo", "specificity", "Embryo specificity"),
           ("event, any overlap", "recall", "Event recall"), ("event, any overlap", "precision", "Event precision"),
           (f"event, overlap >= {I.MIN_RECIPROCAL_OVERLAP:.0%}", "recall", "Strict recall"),
           (f"event, overlap >= {I.MIN_RECIPROCAL_OVERLAP:.0%}", "precision", "Strict precision")]
INTERVAL_NOTE = (f"Point and 95% Wilson interval. Event = a reported region, found with any overlap or, strict, with a "
                 f"reciprocal overlap of at least {I.MIN_RECIPROCAL_OVERLAP:.0%}. The reference is the authors' results in GEO, "
                 f"another algorithm on the same arrays, so this measures agreement.")


def score(setup: str, level: str, metric: str) -> pd.Series:
    return scores[(scores["setup"] == setup) & (scores["level"] == level) & (scores["metric"] == metric)].iloc[0]


def benchmark() -> None:
    s = {m: score(MAIN, lv, m_) for lv, m_, m in METRICS}
    fig, g = figure(f"Benchmark against the authors' results: {MAIN}",
                    f"Embryo sensitivity {s['Embryo sensitivity'].x:.0f}/{s['Embryo sensitivity'].n:.0f} and specificity "
                    f"{s['Embryo specificity'].x:.0f}/{s['Embryo specificity'].n:.0f}. Event recall {pct(s['Event recall'].value)} "
                    f"and precision {pct(s['Event precision'].value)}.",
                    INTERVAL_NOTE, left=0.17, right=0.82)
    ax = fig.add_axes([g["left"], g["bottom"], g["right"] - g["left"], g["top"] - g["bottom"]])
    y = np.array([0, 1, 2.6, 3.6, 5.2, 6.2])
    for yi, (name, r) in zip(y, s.items()):
        ax.plot([r.low, r.high], [yi, yi], color="#555555", lw=2.5, solid_capstyle="butt")
        ax.plot(r.value, yi, "o", color=C["text"], ms=9, clip_on=False)
        ax.text(1.04, yi, f"{r.x:.0f}/{r.n:.0f}  ({pct(r.value)})", va="center", fontsize=13, transform=ax.get_yaxis_transform())
    for yi, head in zip([-0.7, 1.9, 4.5], ["Per embryo", "Per event (any overlap)",
                                           f"Per event (reciprocal overlap ≥ {I.MIN_RECIPROCAL_OVERLAP:.0%})"]):
        ax.text(0, yi, head, weight="bold", fontsize=13, va="center")
    ax.set_yticks(y)
    ax.set_yticklabels([n.replace("Embryo ", "").replace("Event ", "").replace("Strict ", "").capitalize() for n in s],
                       fontsize=13)
    ax.set_ylim(6.7, -1.1)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=12)
    ax.grid(axis="x", color="#E3E2DC", lw=1.5)
    ax.tick_params(length=0)
    ax.spines[["left", "bottom"]].set_visible(False)
    save(fig, "04_benchmark.png")


def setup_marker(ax, setup: str, x, y, **kw) -> None:
    inp = INPUT[setup]
    ax.plot(x, y, INPUT_MARKER[inp], color=CALLER[CALLER_OF[setup]], ms=8,
            mfc=CALLER[CALLER_OF[setup]] if inp == "pipeline" else "white", mew=1.8, clip_on=False, **kw)


def setup_comparison() -> None:
    sm = summary.set_index("setup")
    best = max(SETUPS, key=lambda s: score(s, "event, any overlap", "precision").value
               + score(s, "embryo", "specificity").value)
    fig, g = figure("Six setups against the authors and the parents",
                    f"Only the rules (06c) keep embryo specificity and event precision at 100%, on either input. "
                    f"ASCAT and GBoost find more of the reported regions but add calls the authors do not make.",
                    INTERVAL_NOTE + " Filled = pipeline LRR/BAF, open = CRMA v2. Parents: calls on the two chromosomes of the "
                    "parent's translocation, against chance (karyotypes shuffled among couples).",
                    left=0.2, right=0.97, top=0.8, bottom=0.22, wspace=0.14, hspace=0.5)
    gs = fig.add_gridspec(2, 4, **g)
    y = np.arange(len(SETUPS))
    panels = METRICS + [("parents", None, "On a parent's chromosome"), ("breakpoint", None, "Edge to breakpoint")]
    for k, (level, metric, name) in enumerate(panels):
        ax = fig.add_subplot(gs[k // 4, k % 4])
        for yi, s in zip(y, SETUPS):
            if level == "parents":
                ax.plot([sm.loc[s, "chance_calls"]] * 2, [yi - 0.3, yi + 0.3], color=C["muted"], lw=2)
                setup_marker(ax, s, sm.loc[s, "calls_on_parent"], yi)
            elif level == "breakpoint":
                setup_marker(ax, s, sm.loc[s, "median_breakpoint_mb"], yi)
                ax.text(sm.loc[s, "median_breakpoint_mb"] + 0.8, yi, f"n={sm.loc[s, 'edges_measured']}", va="center", fontsize=8.5,
                        color=C["muted"])
            else:
                r = score(s, level, metric)
                ax.plot([r.low, r.high], [yi, yi], color=CALLER[CALLER_OF[s]], lw=2, alpha=0.6)
                setup_marker(ax, s, r.value, yi)
        ax.set_title(name, fontsize=11.5, weight="bold", loc="left")
        ax.set_ylim(len(SETUPS) - 0.5, -0.5)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{s}  ({sm.loc[s, 'calls']} calls)" for s in SETUPS] if k % 4 == 0 else [], fontsize=10)
        ax.tick_params(axis="y", length=0)
        ax.grid(axis="x", color="#E3E2DC")
        if level == "parents":
            ax.set_xlabel("calls (| = chance)")
            ax.set_xlim(0, sm["calls_on_parent"].max() * 1.1)
        elif level == "breakpoint":
            ax.set_xlabel("median Mb")
            ax.set_xlim(0, sm["median_breakpoint_mb"].max() * 1.35)
        else:
            ax.set_xlim(-0.02, 1.02)
            ax.set_xticks([0, 0.5, 1])
            ax.set_xticklabels(["0", "50%", "100%"])
    handles = [Line2D([], [], marker="o", ls="", color=col, ms=8, label=caller) for caller, col in CALLER.items()]
    handles += [Line2D([], [], marker=mk, ls="", color=C["muted"], mfc=C["muted"] if i == "pipeline" else "white", mew=1.8,
                       ms=8, label=f"{i} LRR/BAF") for i, mk in INPUT_MARKER.items()]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), bbox_to_anchor=(0.55, 0.095), fontsize=11)
    save(fig, "05_setup_comparison.png")


# 06 calls against the parents' karyotypes

def parents_figure(parents: pd.DataFrame) -> None:
    c = calls[calls["setup"] == MAIN]
    sm = summary.set_index("setup")
    key = lambda t: set(zip(t["sample"], t["chr"]))  # noqa: E731
    keep = key(c) | key(authors)
    p = parents[[k in keep for k in zip(parents["sample"], parents["chr"])]].copy()
    p = p.sort_values(["robertsonian", "sample", "chr"])
    rows, y, yi = [], [], 0.0
    for rob, part in p.groupby("robertsonian", sort=True):
        rows.append(("head", rob))
        y.append(yi)
        yi += 1
        for r in part.itertuples():
            rows.append(("row", r))
            y.append(yi)
            yi += 1
        yi += 0.6
    a_, d_ = sm.loc[MAIN], sm.loc["Authors"]
    fig, g = figure(f"Calls against the parents' karyotypes: {MAIN}",
                    f"{a_.calls_on_parent} of {a_.calls} calls fall on a chromosome of the parent's translocation, against "
                    f"{a_.chance_calls:.1f} by chance (p = {a_.p_calls:.4f}). Measurable edges are a median "
                    f"{a_.median_breakpoint_mb:.1f} Mb from the parent's band (authors: {d_.median_breakpoint_mb:.1f} Mb).",
                    "Triangle and bar = the parent's breakpoint band. Outline = authors. Right = distance of the called edge "
                    "from the band; edges at a chromosome end are not measured.",
                    left=0.15, right=0.82, top=0.8, bottom=0.2)
    ax = fig.add_axes([g["left"], g["bottom"], g["right"] - g["left"], g["top"] - g["bottom"]])
    for (what, r), yi in zip(rows, y):
        if what == "head":
            ax.text(0, yi, "Robertsonian (whole chromosome expected)" if r else "Reciprocal (one edge expected at the breakpoint)",
                    weight="bold", fontsize=12.5, va="center")
            continue
        e = ends.loc[r.chr]
        ax.add_patch(Rectangle((e["first"] / 1e6, yi - 0.22), (e["last"] - e["first"]) / 1e6, 0.44, color=C["chrom"], lw=0))
        mine = c[(c["sample"] == r.sample) & (c["chr"] == r.chr)]
        for m in mine.itertuples():
            ax.add_patch(Rectangle((m.start / 1e6, yi - 0.22), (m.end - m.start) / 1e6, 0.44, color=call_colour(m.direction, m.label), lw=0))
        theirs = authors[(authors["sample"] == r.sample) & (authors["chr"] == r.chr)]
        for m in theirs.itertuples():
            ax.add_patch(Rectangle((m.start / 1e6, yi - 0.32), (m.end - m.start) / 1e6, 0.64, fill=False, lw=1.3, edgecolor=C["text"]))
        mid = (r.band_start + r.band_end) / 2e6
        ax.plot([r.band_start / 1e6, r.band_end / 1e6], [yi - 0.4] * 2, color=C["text"], lw=3, solid_capstyle="butt")
        ax.plot(mid, yi - 0.5, "v", color=C["text"], ms=7)
        d = mine["breakpoint_mb"].dropna()
        text = ("missed" if len(theirs) else "not called") if mine.empty else \
            ("whole chromosome" if d.empty else ("in band" if d.min() == 0 else f"{d.min():.1f} Mb off"))
        ax.text(1.02, yi, text, transform=ax.get_yaxis_transform(), va="center", fontsize=11.5)
    ax.set_yticks([yi for (w, _), yi in zip(rows, y) if w == "row"])
    ax.set_yticklabels([f"{r.sample}  {chr_name(r.chr)}{r.band}" for w, r in rows if w == "row"], fontsize=10)
    ax.set_ylim(max(y) + 0.7, -0.7)
    ax.set_xlim(0, ends["last"].max() / 1e6)
    ax.set_xlabel("Position on chromosome (Mb, hg19)")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    call_legend(fig, y=0.04)
    save(fig, "06_parents.png")


# 07 errors against noise

def bin_noise() -> pd.DataFrame:
    """Per input and embryo: median |difference| of neighbouring 50-SNP medians on autosomes, the noise at the
    scale calls are made (a real step adds one big difference and hardly moves the median)."""
    b = bins[bins["chr"] < 23]
    d = b.groupby(["input", "sample", "chr"])["lrr"].diff().abs()
    return d.groupby([b["input"], b["sample"]]).median().rename("noise").reset_index()


def errors_vs_noise() -> None:
    noise = bin_noise()
    fig, g = figure("Where do the errors come from?",
                    "Missed regions and extra calls per embryo, against that embryo's noise. With the rules the few "
                    "errors are misses; ASCAT and GBoost add calls across the noise range.",
                    "Noise = median |difference| of neighbouring 50-SNP LRR medians on the autosomes, per input. Up = calls the "
                    "authors do not report, down = reported regions not found. Labels: the two embryos with most errors (3 or more).",
                    left=0.06, right=0.98, top=0.8, bottom=0.14, wspace=0.15, hspace=0.5)
    gs = fig.add_gridspec(2, 3, **g)
    top = 0
    per = {}
    for s in SETUPS:
        a, c = errors(s)
        e = pd.DataFrame({"missed": a[~a["found"]].groupby("sample").size(), "extra": c[c["extra"]].groupby("sample").size()})
        e = noise[noise["input"] == INPUT[s]].set_index("sample")[["noise"]].join(e).fillna(0)
        per[s] = e
        top = max(top, e["extra"].max(), e["missed"].max())
    for k, s in enumerate(SETUPS):
        ax = fig.add_subplot(gs[k // 3, k % 3])
        e = per[s]
        ax.vlines(e["noise"], 0, e["extra"], color="#D9822B", lw=2)
        ax.vlines(e["noise"], 0, -e["missed"], color=C["loss"], lw=2)
        ax.plot(e["noise"], np.zeros(len(e)), "|", color=C["muted"], ms=5)
        worst = e.assign(n=e["missed"] + e["extra"]).query("n >= 3").nlargest(2, "n")
        for name, r in worst.iterrows():
            ax.text(r["noise"], (r["extra"] if r["extra"] >= r["missed"] else -r["missed"]) * 1.05 + (0.3 if r["extra"] >= r["missed"] else -0.3),
                    name.replace("SNP", ""), fontsize=8, ha="center", va="bottom" if r["extra"] >= r["missed"] else "top")
        ax.axhline(0, color=C["muted"], lw=0.6)
        ax.set_ylim(-top * 1.3 - 0.5, top * 1.3 + 0.5)
        ax.set_title(f"{s}: {int(e['extra'].sum())} extra, {int(e['missed'].sum())} missed", fontsize=11, weight="bold", loc="left")
        ax.set_xlabel(f"noise ({INPUT[s]} LRR)", fontsize=9.5)
        if k % 3 == 0:
            ax.set_ylabel("missed  |  extra")
    save(fig, "07_errors_vs_noise.png")


# 08 BAF by copy number

def baf_by_copy_number(snps: pd.DataFrame) -> None:
    """SNPs of regions whose copy number is known without our calls: the authors' losses and gains, male X
    outside the pseudo-autosomal regions, and the autosomes of embryos the authors call normal."""
    lo, hi = I.HET_BAF
    normal = set(meta.index) - set(authors["sample"])
    male = set(n_x[n_x == 1].index)
    groups = {}
    for inp in [i for i in INPUT_MARKER if i in set(snps["input"])]:
        s = snps[snps["input"] == inp]
        def in_regions(direction):
            m = np.zeros(len(s), bool)
            for r in authors[authors["direction"] == direction].itertuples():
                m |= ((s["sample"] == r.sample) & (s["chr"] == r.chr) & (s["pos"] >= r.start) & (s["pos"] <= r.end)).to_numpy()
            return m
        x_male = (s["sample"].isin(male) & (s["chr"] == 23) & (s["pos"] > ev.PAR1_END) & (s["pos"] < ev.PAR2_START)).to_numpy()
        auto_normal = (s["sample"].isin(normal) & (s["chr"] < 23)).to_numpy()
        groups[inp] = {"1 copy: authors' losses": s["baf"][in_regions(-1)], "1 copy: male X": s["baf"][x_male],
                       "2 copies: normal embryos": s["baf"][auto_normal], "3 copies: authors' gains": s["baf"][in_regions(1)]}
    colours = {"1 copy: authors' losses": C["loss"], "1 copy: male X": C["mosaic_loss"],
               "2 copies: normal embryos": C["muted"], "3 copies: authors' gains": C["gain"]}
    # one column per copy number, with the BAF values a SNP can take there if the measurement were perfect
    columns = [("1 copy", ["1 copy: authors' losses", "1 copy: male X"], [0, 1], "0 and 1"),
               ("2 copies", ["2 copies: normal embryos"], [0, 0.5, 1], "0, ½ and 1"),
               ("3 copies", ["3 copies: authors' gains"], [0, 1 / 3, 2 / 3, 1], "0, ⅓, ⅔ and 1")]
    edges = np.linspace(0, 1, 51)
    hist = {(i, k): np.histogram(v, bins=edges, density=True)[0] for i, g_ in groups.items() for k, v in g_.items()}
    share = {i: {k: ((v > lo) & (v < hi)).mean() for k, v in g_.items()} for i, g_ in groups.items()}
    p = share.get("pipeline", next(iter(share.values())))
    losses, male_x, two, three = (p[k] for k in colours)
    fig, g = figure("Does BAF follow copy number?",
                    "If it did, one copy would leave SNPs only at 0 and 1, two copies would add a peak at ½, and three "
                    "copies peaks at ⅓ and ⅔ (dashed lines). The pipeline's BAF has the peaks at 0 and 1 but none in between, "
                    f"and its middle band holds {pct(losses)} of SNPs in losses, {pct(two)} with two copies and {pct(three)} "
                    "in gains: it follows noise, not copy number. CRMA v2's BAF is one broad hump for every copy number.",
                    f"Every 10th SNP, same axes in every panel. Shaded = BAF {lo}–{hi} (the share in it is in the legend). "
                    "Copy number from the authors' results and the XY embryos' X, not from the calls of any setup.",
                    left=0.08, right=0.98, top=0.76, bottom=0.15, wspace=0.06, hspace=0.12)
    gs = fig.add_gridspec(len(groups), len(columns), **g)
    first = None
    for r, inp in enumerate(groups):
        for c, (title, names, peaks, where) in enumerate(columns):
            ax = fig.add_subplot(gs[r, c], sharex=first, sharey=first)
            first = first or ax
            ax.axvspan(lo, hi, color="#F1F0EC", lw=0)
            for x in peaks:
                ax.axvline(x, color=C["text"], ls="--", lw=1)
            for name in names:
                ax.stairs(hist[(inp, name)], edges, color=colours[name], lw=2,
                          label=f"{name.split(': ')[1]}: {pct(share[inp][name])} in band")
            ax.set_yscale("log")
            ax.legend(fontsize=9.5, loc="upper center", frameon=True, facecolor="white", edgecolor="none", framealpha=1)
            if r == 0:
                ax.set_title(f"{title}: expected at {where}", fontsize=12, weight="bold", loc="left")
            if c == 0:
                ax.set_ylabel(f"{inp} BAF\ndensity (log)")
            else:
                ax.tick_params(labelleft=False)
            if r == len(groups) - 1:
                ax.set_xlabel("BAF")
            else:
                ax.tick_params(labelbottom=False)
    low = min(h[h > 0].min() for h in hist.values())
    first.set_ylim(low * 0.7, max(h.max() for h in hist.values()) * 25)
    first.set_xlim(-0.03, 1.03)
    save(fig, "08_baf_by_copy_number.png")


# 09 X dosage: the one-copy step of each input

def x_dosage() -> None:
    b = bins
    x = b[(b["chr"] == 23) & (b["start"] > ev.PAR1_END) & (b["end"] < ev.PAR2_START)].groupby(["input", "sample"])["lrr"].median()
    auto = b[b["chr"] < 23].groupby(["input", "sample"])["lrr"].median()
    d = (x - auto).rename("x").reset_index()
    d["n_X"] = d["sample"].map(n_x)
    inputs = [i for i in INPUT_MARKER if i in set(d["input"])]
    k = {i: d[(d["input"] == i) & (d["n_X"] == 1)]["x"].median() - d[(d["input"] == i) & (d["n_X"] == 2)]["x"].median()
         for i in inputs}
    fig, g = figure("One copy less: the X of male embryos",
                    "One lost copy should move LRR by −1. It moves it by "
                    + " and ".join(f"{v:+.2f} on the {i} LRR".replace("-", "−") for i, v in k.items())
                    + ": the arrays compress the scale, and each input needs its own calibration.",
                    "Median LRR of X (outside the pseudo-autosomal regions) minus the median of the autosomes, per embryo. "
                    "One or two X from step 04 (arrays.tsv). Line = median of the group. Both panels share the y axis.",
                    left=0.08, right=0.96, wspace=0.2)
    gs = fig.add_gridspec(1, len(inputs), **g)
    rng = np.random.default_rng(1)
    first = None
    for j, inp in enumerate(inputs):
        ax = fig.add_subplot(gs[j], sharey=first)
        first = first or ax
        for pos, nx, col in ((0, 2, C["gain"]), (1, 1, C["loss"])):
            v = d[(d["input"] == inp) & (d["n_X"] == nx)]["x"]
            ax.scatter(pos + rng.uniform(-0.15, 0.15, len(v)), v, s=30, color=col, alpha=0.8, lw=0)
            ax.plot([pos - 0.28, pos + 0.28], [v.median()] * 2, color=C["text"], lw=2)
            ax.text(pos + 0.32, v.median(), f"{v.median() + 0:+.2f}".replace("-0.00", "0.00"), va="center", fontsize=11)
        ax.axhline(0, color=C["muted"], lw=0.6)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"XX (n={int((d[(d['input'] == inp)]['n_X'] == 2).sum())})",
                            f"XY (n={int((d[(d['input'] == inp)]['n_X'] == 1).sum())})"], fontsize=12)
        ax.set_xlim(-0.6, 1.8)
        ax.set_title(f"{inp} LRR: one-copy step {k[inp]:+.2f}", fontsize=12.5, weight="bold", loc="left")
        if j == 0:
            ax.set_ylabel("X − autosomes (median LRR)")
        else:
            ax.tick_params(labelleft=True)
    save(fig, "09_x_dosage.png")


# run

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"figures -> {OUT.relative_to(REPO)}")
    snps = pd.read_csv(T / "profile_snps.tsv.gz", sep="\t",
                       dtype={"input": "category", "sample": "category", "chr": np.int8, "pos": np.int64,
                              "genome_mb": np.float32, "lrr": np.float32, "baf": np.float32})
    parents = ev.load_parents(meta.reset_index(), ev.load_snps())
    for sample, chromosomes in EXAMPLES.items():
        profile(sample, snps)
        if chromosomes:
            close_up(sample, chromosomes, snps, parents)
    for s in SETUPS:
        cnv_map(s, parents)
    benchmark()
    setup_comparison()
    parents_figure(parents)
    errors_vs_noise()
    baf_by_copy_number(snps)
    x_dosage()


if __name__ == "__main__":
    main()
