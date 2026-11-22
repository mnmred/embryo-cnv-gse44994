"""Step 06a - write LRR and BAF as gzipped text tables for the R segmentation (step 06b).

Run:     python scripts/python/06a_export_for_r.py       (~1 min)
Reads:   data/processed/lrr.parquet (step 04), baf.parquet (step 05)
Writes:  data/processed/lrr.tsv.gz, baf.tsv.gz
         probe | chr | pos | one column per array, values rounded to 4 decimals
         (BAF clipped to 0-1, which copynumber::aspcf expects)

Why: R in bioenv has no parquet reader (the R package arrow is not installed because versions need to be downgraded for compatibility), so the
two tables are exported once as text, which R's data.table::fread reads quickly.
"""
import argparse
from pathlib import Path

import pandas as pd

from snparray import paths


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=paths.PROCESSED)
    layout = paths.Layout(ap.parse_args().out_dir)

    for source, target in [(layout.lrr, layout.lrr_tsv), (layout.baf, layout.baf_tsv)]:
        table = pd.read_parquet(source).drop(columns="cytoband")
        if source == layout.baf:                  # copynumber::aspcf expects BAF within 0-1
            values = table.columns[3:]
            table[values] = table[values].clip(0, 1)
        table.to_csv(target, sep="\t", index=False, float_format="%.4f")
        print(f"{source.name}: {table.shape[0]} SNPs x {table.shape[1] - 3} arrays -> {target}")


if __name__ == "__main__":
    main()
