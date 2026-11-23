"""Step 0 (GSE44994 only) - write the paper's result of every embryo as a plain table.

Run (repo root, ml_env active):   python scripts/python/ml/00_paper_results_to_table.py
Reads:   data/processed/metadata.tsv   reported_result, e.g. "46,XX,+9(pter-q33),-12(pter-q23)"
         data/processed/lrr.parquet    only chr, pos, cytoband: to turn bands into base pairs
Writes:  data/ml/reported_abnormal.tsv  one row per changed region:
         sample, chr, start, end, direction (gain / loss), whole (whole chromosome?), copies

This is the only script with its own file names, because it is specific to this dataset.
For other data, write such a table yourself (or leave it out, see inputs.py).
"""
import re

import pandas as pd

import tables

METADATA = "data/processed/metadata.tsv"
SNPS = "data/processed/lrr.parquet"
OUTPUT = "reported_abnormal.tsv"


def band_range(snps: pd.DataFrame, band: str) -> tuple:
    """First and last position of the SNPs of one chromosome that lie in a cytoband
    ("q32" includes q32.1, q32.2 ...). `snps`: pos and band of that chromosome."""
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
        # band not on the array: use the next coarser band
        wanted = wanted.split(".")[0] if "." in wanted else wanted[:-1]


def regions_of_embryo(sample: str, text: str, snps: pd.DataFrame) -> list:
    body = re.sub(r"^([0-9]+,)?(XX|XY)", "", text)
    tokens = [m.group(0) for m in re.finditer(r"[+-]([0-9]{1,2}|X)(\([^)]*\)|[pq])?", body)]
    if len(tokens) != len(re.findall(r"[+-]", re.sub(r"\([^)]*\)", "", body))):
        raise ValueError(f"cannot parse the result of {sample}: {text}")
    rows = []
    for token in dict.fromkeys(tokens):              # each different token once, in order
        name = re.match(r"[+-]([0-9]{1,2}|X)", token).group(1)
        chromosome = 23 if name == "X" else int(name)
        where = token[1 + len(name):].strip("()")
        where = re.sub(r"([0-9])([pq]ter)", r"\1-\2", where)      # GEO has one "q32.1qter"
        here = snps[snps["chr"] == chromosome]
        if where == "":                              # whole chromosome
            start, end = here["pos"].min(), here["pos"].max()
        elif where in ("p", "q"):                    # whole arm
            arm = here[here["band"].str.startswith(where)]
            start, end = arm["pos"].min(), arm["pos"].max()
        else:                                        # from band to band
            ends = [p for band in where.split("-") for p in band_range(here, band)]
            start, end = min(ends), max(ends)
        gain = token.startswith("+")
        rows.append(dict(sample=sample, chr=chromosome, start=int(start), end=int(end),
                         direction="gain" if gain else "loss", whole=where == "",
                         copies=2 + (1 if gain else -1) * tokens.count(token), reported=token))
    return rows


def main():
    metadata = tables.read_table(METADATA)
    snps = tables.read_table(SNPS, ["chr", "pos", "cytoband"]).rename(columns={"cytoband": "band"})
    snps["band"] = snps["band"].astype(str)
    rows = [row for sample, text in zip(metadata["sample"], metadata["reported_result"])
            for row in regions_of_embryo(sample, text, snps)]
    out = pd.DataFrame(rows)
    out.to_csv(tables.output_path(OUTPUT), sep="\t", index=False)
    print(f"{len(out)} reported regions in {out['sample'].nunique()} of {len(metadata)} embryos "
          f"-> data/ml/{OUTPUT}")


if __name__ == "__main__":
    main()
