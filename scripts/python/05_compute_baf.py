"""Step 05 - B-allele frequency (BAF) for every SNP and array.

Run:     python scripts/python/05_compute_baf.py
Reads:   data/processed/signal_A.parquet, signal_B.parquet, genotype_calls.parquet,
         arrays.tsv (X copies, from step 04)
Writes:  data/processed/baf.parquet        probe | chr | pos | cytoband | one BAF column per array
         data/processed/snp_info.parquet   per SNP: AA/AB/BB cluster centres, counts, population B frequency
         data/processed/baf_qc.tsv         per array: agreement of BAF with the deposited genotype calls

Order of operations (functions in snparray/baf.py):
    polar_coordinates -> compute_baf (genotype_clusters + theta_to_baf, leave-one-out)
    -> population_b_frequency
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from snparray import baf, io, paths


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=paths.PROCESSED)
    layout = paths.Layout(ap.parse_args().out_dir)

    snps = io.load_snp_annotation(layout)
    signals = io.load_signals(layout, snps, kinds=("A", "B"))  # load the allele-specific signals from step 03
    samples = list(signals["A"].columns) # the 47 embryo names (the A and B tables have the same columns)
    genotypes = io.load_genotypes(layout, snps, samples)  # loads the GTYPE's output as training labels for the BAF clusters
    n_x = pd.read_csv(layout.arrays, sep="\t", index_col=0).loc[samples, "n_X"].values
    autosomal, is_x = snps.chr.values <= 22, snps.chr.values == 23

    # angle theta gives the A/B mix, whatever thebrightness
    theta, _ = baf.polar_coordinates(signals["A"].values, signals["B"].values)

    # we will ignore the NO CALLs (-1) and the X SNPs in embryos with one X copy (n_x == 1)
    use = genotypes >= 0
    use[np.ix_(is_x, n_x == 1)] = False
    result = baf.compute_baf(theta, genotypes, use, autosomal)

    snp_info = snps[["fsetid", "probe", "chr", "pos"]].assign(
        theta_AA=result["centres"][:, 0], 
        theta_AB=result["centres"][:, 1], 
        theta_BB=result["centres"][:, 2],
        n_AA=result["counts"][:, 0], 
        n_AB=result["counts"][:, 1], 
        n_BB=result["counts"][:, 2]
        )
    b = result["baf"]
    baf_qc = pd.DataFrame({
        "call_rate": (genotypes[autosomal] >= 0).mean(0),
        # concordance: fraction of calls that match the genotype calls (ignoring NO CALLs), this isnt accuracy but agreement!
        "concordance_with_calls": result["concordance"],
        "calls_do_not_match_cel": result["mismatch"],
        "baf_sd_of_AB_calls": [np.std(b[autosomal & (genotypes[:, j] == 1), j]) for j in range(len(samples))],
    }, index=pd.Index(samples, name="sample"))

    io.save_snp_matrix(b, snps, samples, layout.baf)
    snp_info.to_parquet(layout.snp_info)
    baf_qc.to_csv(layout.baf_qc, sep="\t")

    print(f"BAF: {b.shape[0]} SNPs x {b.shape[1]} arrays -> {layout.baf}")
    print(f"SNPs with all three genotype clusters observed: {(result['counts'] >= 3).all(1).mean():.1%}")
    print(f"BAF vs deposited calls: median concordance {np.median(result['concordance']):.3f}")
    bad = [s for s, m in zip(samples, result["mismatch"]) if m]
    print(f"call sets that do not belong to their CEL file (excluded from clustering): {bad or 'none'}")


if __name__ == "__main__":
    main()
