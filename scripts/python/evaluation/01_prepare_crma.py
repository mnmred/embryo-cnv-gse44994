"""Evaluation step 1 - put the CRMA v2 LRR/BAF in the form the end of the main pipeline reads.

The setup "CRMA LRR/BAF + the end of the pipeline" runs the unchanged steps 06b (ASPCF) and
06c (rules) on CRMA's LRR and BAF. Those steps read a folder with fixed file names, so this
script writes such a folder.

Run (repo root, bioenv active):   python scripts/python/evaluation/01_prepare_crma.py   (~1 min)
Reads:   data/ascat/crma/lrr.tsv.gz, baf.tsv.gz, arrays.tsv      (scripts/R/ascat/01_crma_lrr_baf.R)
         data/processed/probe_level/snp_annotation.tsv          (only the cytoband of every SNP)
Writes:  data/evaluation/crma_rules/
         lrr.tsv.gz, baf.tsv.gz   for step 06b   (BAF clipped to 0-1, as step 06a does)
         lrr.parquet              for step 06c   (probe | chr | pos | cytoband | one column per array)
         arrays.tsv               sample, n_X, mapd, lrr_sd
         calibration.tsv          k and alpha, measured on the X of the one-X arrays, as step 04 does

SNPs dropped:
  - the 172 without a cytoband (they are not in the pipeline's SNP table either): step 06c
    groups SNPs by chromosome arm, and one SNP without an arm makes the noise of every
    embryo NaN, so nothing would be called
  - the 10 with a missing BAF
  - the 2 pairs that share a position (step 06c finds SNP rows by position)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "scripts" / "python"))
from snparray import lrr as L, sex  # noqa: E402

CRMA = REPO / "data" / "ascat" / "crma"
ANNOTATION = REPO / "data" / "processed" / "probe_level" / "snp_annotation.tsv"
OUT = REPO / "data" / "evaluation" / "crma_rules"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    lrr = pd.read_csv(CRMA / "lrr.tsv.gz", sep="\t")
    baf = pd.read_csv(CRMA / "baf.tsv.gz", sep="\t")
    arrays = pd.read_csv(CRMA / "arrays.tsv", sep="\t").set_index("sample")
    samples = [c for c in lrr.columns if c not in ("probe", "chr", "pos")]
    if not lrr["probe"].equals(baf["probe"]):
        raise ValueError("lrr.tsv.gz and baf.tsv.gz are not in the same SNP order")

    band = pd.read_csv(ANNOTATION, sep="\t", usecols=["man_fsetid", "cytoband"]).set_index("man_fsetid")["cytoband"]
    lrr.insert(3, "cytoband", band.reindex(lrr["probe"]).values)
    keep = lrr["cytoband"].notna() & ~lrr.duplicated(["chr", "pos"], keep=False)
    keep &= lrr[samples].notna().all(axis=1) & baf[samples].notna().all(axis=1)
    lrr, baf = lrr[keep].reset_index(drop=True), baf[keep].reset_index(drop=True)
    lrr["pos"] = baf["pos"] = lrr["pos"].astype(np.int64)
    print(f"{len(keep)} SNPs, {(~keep).sum()} dropped (no cytoband, missing BAF or shared position), {len(lrr)} kept")

    # the three text tables of step 06a, and the parquet table of step 04
    lrr.drop(columns="cytoband").to_csv(OUT / "lrr.tsv.gz", sep="\t", index=False, float_format="%.4f")
    baf[samples] = baf[samples].clip(0, 1)
    baf.to_csv(OUT / "baf.tsv.gz", sep="\t", index=False, float_format="%.4f")
    lrr.astype({s: np.float32 for s in samples}).to_parquet(OUT / "lrr.parquet")

    # noise and calibration, with the functions of step 04
    values = lrr[samples].to_numpy(np.float64)
    chrom, pos = lrr["chr"].to_numpy(), lrr["pos"].to_numpy()
    n_x = arrays.loc[samples, "n_X"].to_numpy(int)
    out = pd.DataFrame({"n_X": n_x,
                        "mapd": [L.mapd(values[chrom <= 22, j]) for j in range(len(samples))]},
                       index=pd.Index(samples, name="sample"))
    out["lrr_sd"] = L.mapd_to_sd(out["mapd"])
    out.to_csv(OUT / "arrays.tsv", sep="\t")
    k = L.one_copy_shift(values, sex.x_not_on_y_mask(chrom, pos), n_x)
    alpha = L.compression_alpha(k)
    pd.Series({"k": k, "alpha": alpha}).to_csv(OUT / "calibration.tsv", sep="\t", header=False)

    print(f"noise (median): MAPD {out['mapd'].median():.3f}; one lost copy shifts LRR by k = {k:.3f}, alpha = {alpha:.3f}")
    print(f"-> {OUT.relative_to(REPO)}/")


if __name__ == "__main__":
    main()
