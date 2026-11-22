"""Step 06c - give every ASPCF segment of step 06b a copy number (normal, full or mosaic change).

Run:     python scripts/python/06c_call_copy_number.py         (~10 s)
Reads:   data/processed/segments_raw_ASPCF.tsv (step 06b), lrr.parquet, arrays.tsv,
         calibration.tsv (step 04)
Writes:  data/processed/segments_ASPCF.tsv   every joined segment of every array: position,
                                             mean LRR, change, z, final copy number and kind
         data/processed/cnv_calls_ASPCF.tsv  only the reported (abnormal) segments
         data/processed/call_noise.tsv       per array: sigma (50-SNP averages), tau (arms)

Order of operations (functions in snparray/calling.py):
    bin_noise -> between_arm_sd   (once per array)
    segment_rows -> call_segments (per array)
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from snparray import calling, paths

METHOD = "ASPCF"          


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=paths.PROCESSED)
    args = ap.parse_args()
    layout = paths.Layout(args.out_dir)

    lrr_table = pd.read_parquet(layout.lrr)
    arrays = pd.read_csv(layout.arrays, sep="\t", index_col=0)
    alpha = pd.read_csv(layout.calibration, sep="\t", header=None, index_col=0)[1]["alpha"]
    samples = list(arrays.index)
    chrom, pos = lrr_table.chr.values, lrr_table.pos.values
    cytoband = lrr_table.cytoband.astype(str).values
    arm = chrom.astype(str) + pd.Series(cytoband).str[0].values

    # noise of every array: clipped LRR (extreme SNPs pulled in to 4 SD), sigma and tau
    clipped, noise = {}, {}
    for s in samples:
        sd = arrays.loc[s, "lrr_sd"]
        clipped[s] = np.clip(lrr_table[s].values.astype(np.float64), -4 * sd, 4 * sd)
        sigma, bin_means, bin_arm = calling.bin_noise(clipped[s], chrom, arm)
        noise[s] = dict(sigma_bin=sigma, tau_arm=calling.between_arm_sd(bin_means, bin_arm, alpha))
    noise = pd.DataFrame(noise).T
    noise.index.name = "sample"
    noise.to_csv(layout.call_noise, sep="\t")
    print(f"noise (median): sigma {noise.sigma_bin.median():.3f}, tau {noise.tau_arm.median():.3f}")

    # copy number of every segment
    raw = calling.segment_rows(pd.read_csv(layout.segments_raw(METHOD), sep="\t"), chrom, pos)
    segments_per_array = []
    for s in samples:
        called = calling.call_segments(
            segments=raw[raw["sample"] == s],
            y=clipped[s],
            pos=pos,
            cytoband=cytoband,
            sigma=noise.loc[s, "sigma_bin"],
            tau=noise.loc[s, "tau_arm"],
            alpha=alpha,
            n_x=arrays.loc[s, "n_X"],
        )
        segments_per_array.append(called.assign(sample=s))
    segments = pd.concat(segments_per_array, ignore_index=True)

    segments = segments[["sample"] + [c for c in segments.columns if c != "sample"]]
    calls = segments[segments.reported].drop(columns="reported")
    segments.to_csv(layout.segments(METHOD), sep="\t", index=False, float_format="%.4g")
    calls.to_csv(layout.cnv_calls(METHOD), sep="\t", index=False, float_format="%.4g")
    print(f"{len(raw)} ASPCF pieces -> {len(segments)} joined segments -> {layout.segments(METHOD)}")
    print(f"{len(calls)} reported ({(calls.kind == 'full').sum()} full, {(calls.kind == 'mosaic').sum()} mosaic)"
          f" in {calls['sample'].nunique()} of {len(samples)} arrays -> {layout.cnv_calls(METHOD)}")


if __name__ == "__main__":
    main()
