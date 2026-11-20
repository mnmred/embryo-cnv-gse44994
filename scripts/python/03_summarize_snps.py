"""Step 03 - one signal per SNP from its ~12 probes (median polish).

Run:     python scripts/python/03_summarize_snps.py           (~20 s; first run compiles numba)
Reads:   data/processed/probe_level/   (step 02)
Writes:  data/processed/signal_total.parquet   all probes of a SNP     -> used for LRR
         data/processed/signal_A.parquet       allele-A probes         -> used for BAF
         data/processed/signal_B.parquet       allele-B probes         -> used for BAF
         (rows = SNP id "fsetid", columns = GSM id, values = log2 signal)

NOTE on DESIGN:
Step 2 and 3 are the steps done by oligo::snprma, but here we do it in Python with numba for speed.
snprma itself is not used because its preprocessCore dependency crashes in bioenv when it starts threads 
(fixable by reinstalling preprocessCore with --disable-threading, but a lockfile can't capture that).
Two deliberate differences from snprma's output:
     - the two DNA strands are pooled, so A and B each use ~6 probes instead
       of ~3, which means less noise on this amplified few-cell data
     - an extra 'total' signal (all probes of a SNP) is produced for LRR
"""
import argparse
from pathlib import Path

from snparray import io, paths, summarize


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=paths.PROCESSED)
    layout = paths.Layout(ap.parse_args().out_dir)

    x, probes, gsm = io.load_probe_intensities(layout)
    print(f"{x.shape[0]} probes x {x.shape[1]} arrays")
    signals = summarize.summarize_snps(x, probes, gsm)
    for kind, table in signals.items(): # total, A, B
        table.to_parquet(layout.signal(kind))
        print(f"signal_{kind}: {table.shape[0]} SNPs -> {layout.signal(kind)}")


if __name__ == "__main__":
    main()
