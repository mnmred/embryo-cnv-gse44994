"""Step 04 - log R ratio (LRR) for every SNP and array, plus sex and calibration.

Run:     python scripts/python/04_compute_lrr.py
Reads:   data/processed/signal_total.parquet, probe_level/snp_annotation.tsv, metadata.tsv
Writes:  data/processed/lrr.parquet       probe | chr | pos | cytoband | one LRR column per array
         data/processed/arrays.tsv        per array: X copies, reported sex, noise before/after correction
         data/processed/calibration.tsv   k (LRR shift of one lost copy) and alpha (compression)

Order of operations (functions in snparray/sex.py, lrr.py):
    count_x_copies -> reference_signal -> log_r_ratio -> correct_waves -> one_copy_shift
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from snparray import io, lrr, paths, sex


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=paths.PROCESSED)
    layout = paths.Layout(ap.parse_args().out_dir)

    snps = io.load_snp_annotation(layout)
    total = io.load_signals(layout, snps, kinds=("total",))["total"]
    samples = list(total.columns); T = total.values
    chrom, pos = snps.chr.values, snps.pos.values
    autosomal, is_x = chrom <= 22, chrom == 23 #Is it on chrom 1-22 or X?
    x_not_on_y = sex.x_not_on_y_mask(chrom, pos) #Is it on X and outside the two tips that also exist on Y? then it is a probe to count x copies
    arm = chrom.astype(str) + snps.cytoband.str[0].values

    # sex: one or two X chromosomes (from X signal, not from genotype calls)
    n_x = sex.count_x_copies(T, chrom, pos)

    # LRR = signal - reference, then remove amplification bias
    reference = lrr.reference_signal(T, is_x, n_x == 2)
    raw = lrr.log_r_ratio(T, reference)
    # Three facts per SNP that can bias brightness for technical reasons
    neighbourhood_gc = lrr.regional_gc(snps.gc.values, chrom, pos)   # G/C share of the surrounding region
    fragment_length  = snps.frag_len.values                           # length of the DNA piece the SNP is on
    probe_gc         = snps.gc.values                                 # G/C share of the probe itself
    covariates = [neighbourhood_gc, fragment_length, probe_gc]
    # correct the bias linked to those three facts, per embryo
    values = lrr.correct_waves(raw, covariates, autosomal)

    # calibration: how far does LRR move for one lost copy?
    k = lrr.one_copy_shift(values, x_not_on_y, n_x) # thanks to known X copy number, we can see how much LRR moves for one lost copy
    alpha = lrr.compression_alpha(k) # share of the LRR shift that is actually observed (the rest is compressed by the amplification)

    meta = pd.read_csv(layout.metadata, sep="\t").set_index("sample")
    arrays = pd.DataFrame({
        "gsm": meta.loc[samples, "gsm"].values,
        "n_X": n_x,
        "reported_sex": meta.loc[samples, "reported_result"].str.extract(r"(XX|XY)")[0].values,
        "X_median_lrr": np.median(values[x_not_on_y], 0),
        "mapd_before": [lrr.mapd(raw[autosomal, j]) for j in range(len(samples))],
        "mapd": [lrr.mapd(values[autosomal, j]) for j in range(len(samples))],
        "wave_before": [lrr.wave_sd(raw[:, j], chrom, pos, arm) for j in range(len(samples))],
        "wave": [lrr.wave_sd(values[:, j], chrom, pos, arm) for j in range(len(samples))],
    }, index=pd.Index(samples, name="sample"))
    arrays["lrr_sd"] = lrr.mapd_to_sd(arrays.mapd)

    io.save_snp_matrix(values, snps, samples, layout.lrr)
    arrays.to_csv(layout.arrays, sep="\t")
    pd.Series({"k": k, "alpha": alpha}).to_csv(layout.calibration, sep="\t", header=False)

    # evaluation against authors reported sex
    agree = (arrays.reported_sex == np.where(arrays.n_X == 2, "XX", "XY")).sum()
    print(f"LRR: {values.shape[0]} SNPs x {values.shape[1]} arrays -> {layout.lrr}")
    print(f"sex from X dosage agrees with the reported sex for {agree}/{len(arrays)} arrays")
    print(f"noise (median): MAPD {arrays.mapd.median():.3f}, wave SD {arrays.wave_before.median():.3f}"
          f" -> {arrays.wave.median():.3f} after correction")
    print(f"one lost copy shifts LRR by k = {k:.3f} (theory 1.0); alpha = {alpha:.3f}")
    print("expected LRR per copy number:",
          ", ".join(f"{c} cop{'y' if c == 1 else 'ies'} {lrr.expected_lrr(c, alpha):+.3f}" for c in (1, 2, 3, 4)))


if __name__ == "__main__":
    main()
