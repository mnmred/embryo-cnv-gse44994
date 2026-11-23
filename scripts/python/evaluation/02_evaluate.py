"""Evaluation step 2 - score every setup against the authors' results and the parents' karyotypes.

Run (repo root, ml_env or bioenv active):   python scripts/python/evaluation/02_evaluate.py   (~25 s)
Reads:   the call files of the setups in inputs.py (a setup whose file is missing is listed as
         "not run yet"), and the references in inputs.py
Writes:  results/evaluation.md              the tables, for reading
         data/evaluation/ (for figures; the call tables start with setup, input, caller)
           summary.tsv         one row per setup: every headline number with its 95% interval (_low, _high)
           scores.tsv          the same numbers in long form: level, metric, x, n, value, low, high
           calls.tsv           every call of every setup and of the authors: position, direction, overlap
                               with the authors, on a parent chromosome, Mb to the breakpoint
           per_chromosome.tsv  every source x embryo x chromosome: state -1 / 0 / +1 (2 = gain and loss)
           per_embryo.tsv      one row per embryo: gained and lost chromosomes, one column per source
           profile_bins.tsv.gz every input x embryo x 50-SNP bin: median LRR, share of heterozygous-looking
                               SNPs, median |BAF - 0.5| of those
           profile_snps.tsv.gz every input x embryo x 10th SNP: LRR and BAF, for scatter plots
         Positions: `genome_mb` lays the chromosomes end to end (1 first, X last), the same in every table.

The same rules for every setup:
  - calls shorter than MIN_CALL_MB are dropped, and touching calls of the same direction are
    joined (this matters only for ASCAT, which reports every copy-number step as its own segment)
  - a call and an authors' region match when they are in the same embryo, on the same chromosome,
    in the same direction (gain / loss), and overlap. "Strict": reciprocal overlap >= 50%.
  - mosaic and full changes count the same: only the direction matters.

Three levels against the authors:
  embryo       47 embryos: abnormal (any change) or not.  What a clinic would decide on.
  chromosome   47 x 23 chromosomes: changed or not.       Gives a specificity with many negatives.
  event        every region: found / confirmed.           Where the changes are.
Against the parents (no truth, only predictions): how many calls fall on the two chromosomes of
the parent's translocation, against chance, and how far their edges are from its breakpoints.
"""
import re
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))
import inputs as I  # noqa: E402

PAR1_END, PAR2_START = 2_699_520, 154_931_044      # hg19: X outside these two positions is not on Y
CHROMOSOMES = list(range(1, 24))                    # 23 = X
CALL_COLUMNS = ["sample", "chr", "start", "end", "direction", "label"]


def chr_name(c) -> str:
    return "X" if int(c) == 23 else str(int(c))


def at(p) -> Path:
    return REPO / p


# references

def load_snps() -> pd.DataFrame:
    snps = pd.read_parquet(at(I.SNPS), columns=["chr", "pos", "cytoband"]).rename(columns={"cytoband": "band"})
    return snps.astype({"pos": np.int64, "band": str})


def band_range(snps: pd.DataFrame, band: str) -> tuple:
    """First and last SNP position of a cytoband on one chromosome ("q32" includes q32.1 ...).
    The same rules as band_range() in scripts/python/ml/00_paper_results_to_table.py."""
    if band == "pter":
        return snps["pos"].min(), snps["pos"].min()
    if band == "qter":
        return snps["pos"].max(), snps["pos"].max()
    on_p = snps["band"].str.startswith("p")
    if band in ("p10", "q10"):                       # the centromere
        return snps.loc[~on_p, "pos"].min(), snps.loc[~on_p, "pos"].min()
    if band.startswith("p") and not on_p.any():      # no SNPs on the short arm (13, 14, 15, 22)
        return snps["pos"].min(), snps["pos"].min()
    wanted = band
    while True:
        if "." in wanted or len(wanted) <= 2:
            hit = snps["band"].str.startswith(wanted)
        else:
            hit = (snps["band"] == wanted) | snps["band"].str.startswith(wanted + ".")
        if hit.any():
            return snps.loc[hit, "pos"].min(), snps.loc[hit, "pos"].max()
        if len(wanted) <= 1:
            raise ValueError(f"cannot place band {band}")
        wanted = wanted.split(".")[0] if "." in wanted else wanted[:-1]


def load_authors() -> pd.DataFrame:
    a = pd.read_csv(at(I.AUTHORS), sep="\t")
    a["direction"] = np.where(a["direction"] == "gain", 1, -1)
    a["label"] = a["reported"]
    return a[CALL_COLUMNS].reset_index(drop=True)


PARENT = re.compile(r"(t|der|rob)\((\d+);(\d+)\)\(([pq][\d.]+);([pq][\d.]+)")


def load_parents(meta: pd.DataFrame, snps: pd.DataFrame) -> pd.DataFrame:
    """One row per embryo and chromosome of the parent's translocation, with the breakpoint band in bp.
    Robertsonian = both breakpoints q10 (GEO labels some of them 'reciprocal')."""
    rows = []
    for s, karyotype in zip(meta["sample"], meta["donor_karyotype"]):
        m = PARENT.search(karyotype)
        if m is None:
            raise ValueError(f"cannot parse the parent's karyotype of {s}: {karyotype}")
        chrs, bands = (int(m[2]), int(m[3])), (m[4], m[5])
        robertsonian = bands == ("q10", "q10")
        for c, b in zip(chrs, bands):
            lo, hi = band_range(snps[snps["chr"] == c], b)
            rows.append(dict(sample=s, chr=c, band=b, band_start=lo, band_end=hi, robertsonian=robertsonian,
                             parent=f"{'rob' if robertsonian else 't'}({chrs[0]};{chrs[1]})"))
    return pd.DataFrame(rows)


# calls of a setup

def read_rules(file) -> pd.DataFrame:
    c = pd.read_csv(file, sep="\t")
    return c.assign(direction=np.sign(c["change"]).astype(int), label=c["kind"])


def read_ml(file) -> pd.DataFrame:
    c = pd.read_csv(file, sep="\t")
    return c.assign(direction=np.where(c["direction"] == "gain", 1, -1), label=c["state"])


def read_ascat(file, n_x: dict) -> pd.DataFrame:
    """Every ASCAT segment -> gains and losses against the normal copy number there
    (1 on the X of a one-X embryo outside the regions also on Y, else 2)."""
    seg = pd.read_csv(file, sep="\t").sort_values(["sample", "chr", "start"]).reset_index(drop=True)
    x_only = (seg["chr"] == 23) & ~((seg["end"] < PAR1_END) | (seg["start"] > PAR2_START))
    normal = np.where(x_only & seg["sample"].map(n_x).eq(1), 1, 2)
    seg["direction"] = np.sign(seg["copies"] - normal).fillna(0).astype(int)
    seg["label"] = seg["copies"].astype("Int64").astype(str) + " copies"
    return seg[seg["direction"] != 0]


def join_touching(c: pd.DataFrame, snps: pd.DataFrame) -> pd.DataFrame:
    """Join calls of one embryo, chromosome and direction that have no SNP between them."""
    c = c.sort_values(["sample", "chr", "start"]).reset_index(drop=True)
    if c.empty:
        return c
    snp_key = (snps["chr"].to_numpy(np.int64) * 10**10 + snps["pos"].to_numpy(np.int64))
    snp_key.sort()
    key = lambda chrom, pos: np.asarray(chrom, np.int64) * 10**10 + np.asarray(pos, np.int64)  # noqa: E731
    # SNPs strictly between the end of the previous call and the start of this one
    between = (np.searchsorted(snp_key, key(c["chr"], c["start"]), "left")
               - np.searchsorted(snp_key, key(c["chr"].shift(fill_value=0), c["end"].shift(fill_value=0)), "right"))
    same = ((c["sample"] == c["sample"].shift()) & (c["chr"] == c["chr"].shift())
            & (c["direction"] == c["direction"].shift()) & (between <= 0))
    run = (~same).cumsum()
    return c.groupby(run).agg(sample=("sample", "first"), chr=("chr", "first"), start=("start", "min"),
                              end=("end", "max"), direction=("direction", "first"),
                              label=("label", lambda v: "/".join(dict.fromkeys(map(str, v))))).reset_index(drop=True)


def load_setup(setup: dict, n_x: dict, snps: pd.DataFrame) -> pd.DataFrame:
    file = at(setup["calls"])
    c = {"rules": read_rules, "ml": read_ml}[setup["reader"]](file) if setup["reader"] != "ascat" \
        else read_ascat(file, n_x)
    c = join_touching(c[CALL_COLUMNS], snps)
    c["size_mb"] = (c["end"] - c["start"]) / 1e6
    return c[c["size_mb"] >= I.MIN_CALL_MB].reset_index(drop=True)


#  scores

def wilson(x: int, n: int, z: float = 1.959964) -> tuple:
    if n == 0:
        return np.nan, np.nan
    p = x / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, centre - half), min(1.0, centre + half)


def score(level: str, metric: str, x: int, n: int) -> dict:
    lo, hi = wilson(x, n)
    return dict(level=level, metric=metric, x=int(x), n=int(n), value=x / n if n else np.nan, low=lo, high=hi)


def two_by_two(level: str, truth: np.ndarray, called: np.ndarray) -> list:
    tp, fp = int((truth & called).sum()), int((~truth & called).sum())
    fn, tn = int((truth & ~called).sum()), int((~truth & ~called).sum())
    return [score(level, "sensitivity", tp, tp + fn), score(level, "specificity", tn, tn + fp),
            score(level, "PPV", tp, tp + fp), score(level, "NPV", tn, tn + fn),
            score(level, "accuracy", tp + tn, tp + fp + fn + tn),
            dict(level=level, metric="F1", value=2 * tp / (2 * tp + fp + fn) if tp + fp + fn else np.nan),
            *[dict(level=level, metric=k, value=v) for k, v in dict(TP=tp, FP=fp, FN=fn, TN=tn).items()]]


def best_overlap(x: pd.DataFrame, y: pd.DataFrame) -> np.ndarray:
    """For every row of x: the best reciprocal overlap with a row of y of the same embryo, chromosome
    and direction (0 = none). Reciprocal overlap = shared length / length of the longer of the two."""
    m = x.reset_index().merge(y, on=["sample", "chr", "direction"], suffixes=("", "_y"))
    shared = (np.minimum(m["end"], m["end_y"]) - np.maximum(m["start"], m["start_y"]) + 1).clip(lower=0)
    longest = np.maximum(m["end"] - m["start"], m["end_y"] - m["start_y"]) + 1
    return (shared / longest).groupby(m["index"]).max().reindex(range(len(x)), fill_value=0).to_numpy()


def summary(calls: pd.DataFrame) -> dict:
    """embryo -> its changes per chromosome, e.g. '+9 -12' ('' = normal)."""
    out = {}
    for s, c in calls.groupby("sample"):
        parts = sorted({(int(ch), int(d)) for ch, d in zip(c["chr"], c["direction"])})
        out[s] = " ".join(f"{'+' if d > 0 else '−'}{chr_name(ch)}" for ch, d in parts)
    return out


def against_authors(calls: pd.DataFrame, authors: pd.DataFrame, samples: list) -> tuple:
    rows = []
    # embryo
    truth = np.isin(samples, authors["sample"])
    called = np.isin(samples, calls["sample"])
    rows += two_by_two("embryo", truth, called)
    a, c = summary(authors), summary(calls)
    rows.append(score("embryo", "same result as the authors", sum(a.get(s, "") == c.get(s, "") for s in samples),
                      len(samples)))
    # chromosome
    units = pd.MultiIndex.from_product([samples, CHROMOSOMES], names=["sample", "chr"])
    dirs = lambda t: t.groupby(["sample", "chr"])["direction"].agg(set).reindex(units)  # noqa: E731
    t_dirs, c_dirs = dirs(authors), dirs(calls)
    t_pos, c_pos = t_dirs.notna().to_numpy(), c_dirs.notna().to_numpy()
    rows += two_by_two("chromosome", t_pos, c_pos)
    wrong = sum(not (t & c) for t, c in zip(t_dirs[t_pos & c_pos], c_dirs[t_pos & c_pos]))
    rows.append(dict(level="chromosome", metric="changed on both sides, opposite direction", value=wrong))
    # event
    calls = calls.assign(overlap=best_overlap(calls, authors))
    authors = authors.assign(overlap=best_overlap(authors, calls))
    for name, threshold in (("any overlap", 0), (f"overlap >= {I.MIN_RECIPROCAL_OVERLAP:.0%}", I.MIN_RECIPROCAL_OVERLAP)):
        hit_t = authors["overlap"] > threshold if threshold == 0 else authors["overlap"] >= threshold
        hit_c = calls["overlap"] > threshold if threshold == 0 else calls["overlap"] >= threshold
        level = f"event, {name}"
        recall, precision = score(level, "recall", hit_t.sum(), len(authors)), score(level, "precision", hit_c.sum(), len(calls))
        r, p = recall["value"], precision["value"]
        rows += [recall, precision, dict(level=level, metric="F1", value=2 * r * p / (r + p) if r + p > 0 else np.nan)]
    return rows, calls, authors


def on_parent(calls: pd.DataFrame, parents: pd.DataFrame) -> np.ndarray:
    return pd.MultiIndex.from_frame(calls[["sample", "chr"]]).isin(pd.MultiIndex.from_frame(parents[["sample", "chr"]]))


def breakpoint_distance(calls: pd.DataFrame, parents: pd.DataFrame, snps: pd.DataFrame) -> np.ndarray:
    """Mb from the nearest edge of each call to the parent's breakpoint band on that chromosome (0 = inside).
    Only edges inside the chromosome count (not at its first or last SNP). NaN when there is nothing to
    measure: the call is not on a parent chromosome, it covers the whole chromosome, or the breakpoint
    itself is at the end of the array's data (Robertsonian q10 and 22p11.2: every call there touches it)."""
    ends = snps.groupby("chr")["pos"].agg(first="min", last="max")
    tol = I.END_TOLERANCE_MB * 1e6
    bp = parents.set_index(["sample", "chr"])
    out = np.full(len(calls), np.nan)
    for i, c in enumerate(calls.itertuples()):
        if (c.sample, c.chr) not in bp.index:
            continue
        b = bp.loc[(c.sample, c.chr)]
        first, last = ends.loc[c.chr, "first"], ends.loc[c.chr, "last"]
        if min(b.band_start - first, last - b.band_end) <= tol:
            continue
        edges = [e for e in (c.start, c.end) if min(e - first, last - e) > tol]
        if edges:
            out[i] = min(max(0, b.band_start - e, e - b.band_end) for e in edges) / 1e6
    return out


def against_parents(calls: pd.DataFrame, parents: pd.DataFrame, meta: pd.DataFrame, snps: pd.DataFrame) -> tuple:
    calls = calls.assign(on_parent=on_parent(calls, parents),
                         breakpoint_mb=breakpoint_distance(calls, parents, snps))
    # chance: give every couple the karyotype of a random couple (shuffled among couples)
    couple_of = dict(zip(meta["sample"], meta["couple"]))
    pair = parents.groupby("sample")["chr"].agg(list)
    couples = sorted(set(couple_of.values()))
    pairs = np.array([pair[next(s for s in meta["sample"] if couple_of[s] == k)] for k in couples])   # couples x 2
    rng = np.random.default_rng(I.SEED)
    perm = np.argsort(rng.random((I.PERMUTATIONS, len(couples))), axis=1)
    idx = np.array([couples.index(couple_of[s]) for s in calls["sample"]], dtype=int)
    chrs = calls["chr"].to_numpy()
    hit = (pairs[perm[:, idx], 0] == chrs) | (pairs[perm[:, idx], 1] == chrs)     # permutations x calls
    sample_codes, sample_idx = np.unique(calls["sample"], return_inverse=True)
    hit_embryo = np.column_stack([hit[:, sample_idx == j].any(axis=1) for j in range(len(sample_codes))]) \
        if len(calls) else np.zeros((I.PERMUTATIONS, 0), bool)
    observed_calls = int(calls["on_parent"].sum())
    observed_embryos = calls.loc[calls["on_parent"], "sample"].nunique()
    p = lambda chance, obs: (1 + (chance >= obs).sum()) / (1 + len(chance))  # noqa: E731
    d = calls["breakpoint_mb"].dropna()
    row = dict(calls=len(calls), calls_on_parent=observed_calls, chance_calls=hit.sum(1).mean(),
               p_calls=p(hit.sum(1), observed_calls),
               embryos_called=calls["sample"].nunique(), embryos_on_parent=observed_embryos,
               chance_embryos=hit_embryo.sum(1).mean(), p_embryos=p(hit_embryo.sum(1), observed_embryos),
               edges_measured=len(d), median_breakpoint_mb=d.median() if len(d) else np.nan,
               within_5mb=int((d <= 5).sum()))
    return row, calls


# md file

def fmt(v, digits=2):
    return "–" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{v:.{digits}f}"


def ci(r, digits=2) -> str:
    return f"{fmt(r['value'], digits)} ({fmt(r['low'], digits)}–{fmt(r['high'], digits)})"


def md_table(header: list, rows: list, left: int = 1) -> str:
    """`left`: how many columns, from the first, are text (left-aligned); the others are numbers."""
    header = list(map(str, header))
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" if i < left else "--:" for i in range(len(header))) + "|"]
    return "\n".join(lines + ["| " + " | ".join(map(str, r)) + " |" for r in rows])


def wide_summary(done: list, scores: pd.DataFrame, parents_rows: dict) -> pd.DataFrame:
    """One row per setup (and one for the authors, parent columns only): every headline number, with
    95% intervals as <name>_low / <name>_high. The table behind the report, made for plotting."""
    strict = f"event, overlap >= {I.MIN_RECIPROCAL_OVERLAP:.0%}"
    wanted = {"embryo_sensitivity": ("embryo", "sensitivity"), "embryo_specificity": ("embryo", "specificity"),
              "embryo_ppv": ("embryo", "PPV"), "embryo_npv": ("embryo", "NPV"),
              "same_result": ("embryo", "same result as the authors"),
              "chromosome_sensitivity": ("chromosome", "sensitivity"), "chromosome_specificity": ("chromosome", "specificity"),
              "chromosome_ppv": ("chromosome", "PPV"), "chromosome_npv": ("chromosome", "NPV"),
              "event_recall": ("event, any overlap", "recall"), "event_precision": ("event, any overlap", "precision"),
              "strict_recall": (strict, "recall"), "strict_precision": (strict, "precision")}
    counts = {f"{lv}_{k}": (lv, k) for lv in ("embryo", "chromosome") for k in ("TP", "FP", "FN", "TN")}
    rows = []
    for s in done:
        here = scores[scores.setup == s["name"]].set_index(["level", "metric"])
        row = dict(setup=s["name"], input=s["input"], caller=s["caller"])
        for name, key in wanted.items():
            row |= {name: here.loc[key, "value"], f"{name}_low": here.loc[key, "low"], f"{name}_high": here.loc[key, "high"]}
        row |= {name: int(here.loc[key, "value"]) for name, key in counts.items()}
        row |= parents_rows[s["name"]]
        rows.append(row)
    rows.append(dict(setup="Authors", input="", caller="CNAT 4.0.1") | parents_rows["Authors"])
    return pd.DataFrame(rows)


def write_report(missing: list, summary_table: pd.DataFrame, per_embryo: pd.DataFrame,
                 authors: pd.DataFrame, n_embryos: int) -> str:
    t = summary_table[summary_table.setup != "Authors"]
    names = list(t.setup)
    out = [f"# Evaluation\n",
           f"`scripts/python/evaluation/02_evaluate.py`, {date.today()}. Authors: {len(authors)} regions in "
           f"{authors['sample'].nunique()}/{n_embryos} embryos. Calls < {I.MIN_CALL_MB} Mb dropped; gain/loss only "
           f"(mosaic = full); ASCAT `{I.ASCAT_SETTING}`. Intervals and more columns: `{I.TABLE_FOLDER}/summary.tsv`.\n"]
    if missing:
        out.append("Not run yet: " + ", ".join(f"{s['name']} ({s['run']})" for s in missing) + ".\n")
    out += ["## Against the authors\n",
            md_table(["Setup", "Calls", "Embryo sens.", "Embryo spec.", "Same result", "Chr. sens.", "Chr. spec.",
                      "Recall", "Precision", "Strict recall", "Strict precision"],
                     [[r.setup, r.calls, fmt(r.embryo_sensitivity), fmt(r.embryo_specificity),
                       f"{round(r.same_result * n_embryos)}/{n_embryos}", fmt(r.chromosome_sensitivity),
                       fmt(r.chromosome_specificity, 3), fmt(r.event_recall), fmt(r.event_precision),
                       fmt(r.strict_recall), fmt(r.strict_precision)] for r in t.itertuples()]),
            "",
            "## Counts (TP / FP / FN / TN)\n",
            md_table(["Setup", "Embryos", "Chromosomes"],
                     [[r.setup, " / ".join(str(int(getattr(r, f"embryo_{k}"))) for k in ("TP", "FP", "FN", "TN")),
                       " / ".join(str(int(getattr(r, f"chromosome_{k}"))) for k in ("TP", "FP", "FN", "TN"))]
                      for r in t.itertuples()]),
            "",
            "## Against the parents\n",
            md_table(["Setup", "Calls on parent chr.", "Chance", "p", "Median Mb to breakpoint", "Edges"],
                     [[r.setup, f"{r.calls_on_parent}/{r.calls}", fmt(r.chance_calls, 1), fmt(r.p_calls, 4),
                       fmt(r.median_breakpoint_mb, 1), r.edges_measured] for r in summary_table.itertuples()]),
            "",
            "## Per embryo (`=` same as the authors)\n",
            md_table(["Embryo", "Parent", "Authors"] + names,
                     [[r["sample"], r["parent"], r["authors"] or "normal",
                       *("=" if r[n] == r["authors"] else (r[n] or "normal") for n in names)]
                      for _, r in per_embryo.iterrows()], left=3 + len(names)),
            ""]
    return "\n".join(out)


def per_chromosome(samples: list, tables: dict) -> pd.DataFrame:
    """One row per source, embryo and chromosome: state -1 (loss), 0 (normal), +1 (gain), 2 (both).
    For heatmaps of who calls what."""
    units = pd.MultiIndex.from_product([samples, CHROMOSOMES], names=["sample", "chr"])
    out = []
    for name, calls in tables.items():
        d = calls.groupby(["sample", "chr"])["direction"].agg(lambda v: int(v.iloc[0]) if v.nunique() == 1 else 2)
        out.append(d.reindex(units, fill_value=0).rename("state").reset_index().assign(setup=name))
    return pd.concat(out)[["setup", "sample", "chr", "state"]]


#LRR and BAF for figures

def genome_offsets(snps: pd.DataFrame) -> pd.Series:
    """Mb to add to a position on each chromosome to lay them end to end (1 first, X last)."""
    return (snps.groupby("chr")["pos"].max() / 1e6).cumsum().shift(fill_value=0.0)


def bin_ids(chrom: np.ndarray, size: int) -> np.ndarray:
    """`size` neighbouring SNPs of one chromosome per bin; the last bin of a chromosome takes the
    leftover SNPs (as in snparray/calling.py)."""
    out, n = np.empty(len(chrom), np.int64), 0
    edges = np.r_[0, np.flatnonzero(np.diff(chrom)) + 1, len(chrom)]
    for a, e in zip(edges[:-1], edges[1:]):
        k = max(1, (e - a) // size)
        out[a:e] = n + np.minimum(np.arange(e - a) // size, k - 1)
        n += k
    return out


def read_matrix(path, samples: list) -> pd.DataFrame:
    p = at(path)
    t = pd.read_parquet(p, columns=["chr", "pos"] + samples) if p.suffix == ".parquet" \
        else pd.read_csv(p, sep="\t", usecols=["chr", "pos"] + samples)
    return t[["chr", "pos"] + samples]


def write_profiles(samples: list, offsets: pd.Series, folder: Path) -> None:
    bins, points = [], []
    for name, files in I.PROFILES.items():
        if not (at(files["lrr"]).exists() and at(files["baf"]).exists()):
            print(f"profiles: no LRR/BAF for {name}, skipped")
            continue
        lrr, baf = read_matrix(files["lrr"], samples), read_matrix(files["baf"], samples)
        if not (lrr["chr"].equals(baf["chr"]) and lrr["pos"].equals(baf["pos"])):
            raise ValueError(f"LRR and BAF of {name} are not in the same SNP order")
        chrom, pos = lrr["chr"].to_numpy(), lrr["pos"].to_numpy(np.int64)
        genome = np.round(pos / 1e6 + offsets.reindex(chrom).to_numpy(), 4)
        L = lrr[samples].astype(np.float32)
        B = baf[samples].clip(0, 1).astype(np.float32)
        het = (B > I.HET_BAF[0]) & (B < I.HET_BAF[1])
        b = bin_ids(chrom, I.PROFILE_BIN_SNPS)
        where = pd.DataFrame({"chr": chrom, "pos": pos, "genome": genome}).groupby(b).agg(
            chr=("chr", "first"), start=("pos", "min"), end=("pos", "max"), n_snps=("pos", "size"),
            genome_mb=("genome", "mean"))
        values = {"lrr": L.groupby(b).median(), "het_share": het.groupby(b).mean(),
                  "baf_offset": (B - 0.5).abs().where(het).groupby(b).median()}
        nb = len(where)
        t = pd.DataFrame({"input": name, "sample": np.repeat(samples, nb)})
        for col in where.columns:
            t[col] = np.tile(where[col].to_numpy(), len(samples))
        for col, m in values.items():
            t[col] = m[samples].to_numpy().T.ravel()           # sample by sample, bins in genome order
        bins.append(t)
        keep = np.arange(0, len(chrom), I.PROFILE_EVERY_SNP)
        points.append(pd.DataFrame({"input": name, "sample": np.repeat(samples, len(keep)),
                                    "chr": np.tile(chrom[keep], len(samples)), "pos": np.tile(pos[keep], len(samples)),
                                    "genome_mb": np.tile(genome[keep], len(samples)),
                                    "lrr": L.to_numpy()[keep].T.ravel(), "baf": B.to_numpy()[keep].T.ravel()}))
    if bins:
        pd.concat(bins).round({"genome_mb": 3, "lrr": 4, "het_share": 3, "baf_offset": 4}).to_csv(
            folder / "profile_bins.tsv.gz", sep="\t", index=False)
        pd.concat(points).round({"lrr": 4, "baf": 4}).to_csv(folder / "profile_snps.tsv.gz", sep="\t", index=False)


def main():
    meta = pd.read_csv(at(I.METADATA), sep="\t", dtype={"couple": str})
    samples = list(meta["sample"])
    n_x = pd.read_csv(at(I.ARRAYS), sep="\t").set_index("sample")["n_X"].to_dict()
    snps = load_snps()
    authors = load_authors()
    parents = load_parents(meta, snps)

    done = [s for s in I.SETUPS if at(s["calls"]).exists()]
    missing = [s for s in I.SETUPS if s not in done]
    for s in missing:
        print(f"not run yet: {s['name']} ({s['calls']} missing) -> {s['run']}")

    scores, all_calls, parents_rows = [], [], {}
    parents_rows["Authors"], _ = against_parents(authors, parents, meta, snps)
    per_embryo = pd.DataFrame({"sample": samples})
    per_embryo["parent"] = per_embryo["sample"].map(parents.groupby("sample")["parent"].first())
    per_embryo["authors"] = per_embryo["sample"].map(summary(authors)).fillna("")
    for s in done:
        calls = load_setup(s, n_x, snps)
        rows, calls, _ = against_authors(calls, authors, samples)
        parents_rows[s["name"]], calls = against_parents(calls, parents, meta, snps)
        scores += [dict(setup=s["name"], **r) for r in rows]
        all_calls.append(calls.assign(setup=s["name"]))
        per_embryo[s["name"]] = per_embryo["sample"].map(summary(calls)).fillna("")
        print(f"{s['name']:22s} {len(calls):4d} calls")

    scores = pd.DataFrame(scores)
    info = {s["name"]: dict(input=s["input"], caller=s["caller"]) for s in done} | {"Authors": dict(input="", caller="CNAT 4.0.1")}
    tag = lambda t: t.assign(input=t["setup"].map(lambda n: info[n]["input"]),  # noqa: E731
                             caller=t["setup"].map(lambda n: info[n]["caller"]))
    folder = at(I.TABLE_FOLDER)
    folder.mkdir(parents=True, exist_ok=True)

    summary_table = wide_summary(done, scores, parents_rows)
    summary_table.to_csv(folder / "summary.tsv", sep="\t", index=False, float_format="%.4g")
    tag(scores)[["setup", "input", "caller", "level", "metric", "x", "n", "value", "low", "high"]].to_csv(
        folder / "scores.tsv", sep="\t", index=False, float_format="%.4g")
    authors_rows = against_parents(authors.assign(overlap=np.nan),
                                   parents, meta, snps)[1]
    calls = pd.concat([authors_rows.assign(setup="Authors", size_mb=(authors_rows.end - authors_rows.start) / 1e6)]
                      + all_calls, ignore_index=True)
    offsets = genome_offsets(snps)
    calls["genome_start_mb"] = (calls["start"] / 1e6 + calls["chr"].map(offsets)).round(3)
    calls["genome_end_mb"] = (calls["end"] / 1e6 + calls["chr"].map(offsets)).round(3)
    cols = ["setup", "input", "caller", "sample", "chr", "start", "end", "genome_start_mb", "genome_end_mb",
            "size_mb", "direction", "label", "overlap", "on_parent", "breakpoint_mb"]
    tag(calls)[cols].round({"size_mb": 2, "overlap": 4, "breakpoint_mb": 2}).to_csv(
        folder / "calls.tsv", sep="\t", index=False)
    by_source = {"Authors": authors} | {c["setup"].iloc[0]: c for c in all_calls if len(c)}
    tag(per_chromosome(samples, by_source))[["setup", "input", "caller", "sample", "chr", "state"]].to_csv(folder / "per_chromosome.tsv", sep="\t", index=False)
    per_embryo.to_csv(folder / "per_embryo.tsv", sep="\t", index=False)

    if I.WRITE_PROFILES:
        write_profiles(samples, offsets, folder)
    at(I.REPORT).write_text(write_report(missing, summary_table, per_embryo, authors, len(samples)), encoding="utf-8")
    print(f"-> {I.REPORT}; {I.TABLE_FOLDER}/summary.tsv, scores.tsv, calls.tsv, per_chromosome.tsv, per_embryo.tsv"
          + (", profile_bins.tsv.gz, profile_snps.tsv.gz" if I.WRITE_PROFILES else ""))

if __name__ == "__main__":
    main()
