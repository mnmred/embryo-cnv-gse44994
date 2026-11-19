"""Step 01 - sample sheet and the authors' genotype calls from the GEO series matrix.

Run (repo root, bioenv active):   python scripts/python/01_parse_metadata.py
Reads:   data/meta/GSE44994_series_matrix.txt
Writes:  data/processed/metadata.tsv          one row per embryo (47)
         data/processed/genotype_calls.parquet AA/AB/BB calls per SNP and array
"""
import argparse
from pathlib import Path

from snparray import io, paths


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--series-matrix", type=Path, default=paths.SERIES_MATRIX)
    ap.add_argument("--out-dir", type=Path, default=paths.PROCESSED)
    args = ap.parse_args()
    layout = paths.Layout(args.out_dir)
    layout.root.mkdir(parents=True, exist_ok=True)

    meta = io.read_sample_sheet(args.series_matrix)
    meta.to_csv(layout.metadata, sep="\t", index=False)
    print(f"{len(meta)} arrays -> {layout.metadata}")
    print(meta["reported_class"].value_counts().to_string())

    gt = io.read_genotype_calls(args.series_matrix)
    gt.to_parquet(layout.genotype_calls)
    print(f"genotype calls {gt.shape} (SNPs x arrays), median call rate "
          f"{(gt >= 0).mean().median():.2f} -> {layout.genotype_calls}")


if __name__ == "__main__":
    main()
