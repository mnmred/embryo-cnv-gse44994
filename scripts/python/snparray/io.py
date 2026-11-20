"""
Reading and writing the pipeline's files.
"""
import csv
from pathlib import Path

import numpy as np
import pandas as pd

from .paths import Layout


# series matrix
def read_sample_sheet(series_matrix: Path) -> pd.DataFrame:
    """One row per array: GSM id, sample name, couple id, parent's karyotype,
    the authors' result for the embryo."""
    rows = {}
    with open(series_matrix, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            # parse the line into a key and values, then store in rows
            key, *vals = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
            key = key[len("!Sample_"):]
            if key == "characteristics_ch1":                 # "tag: value" fields
                tag = vals[0].split(":", 1)[0].strip()
                rows[tag] = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
            elif key in ("title", "geo_accession"):
                rows[key] = vals
    meta = pd.DataFrame(rows).rename(columns={
        "geo_accession": "gsm", "title": "sample",
        "karyotype of patient donor": "donor_karyotype",
        "karyotype classification": "donor_class",
        "analysis result": "reported_result",
        "result classification": "reported_class"})
   # extract the couple id from the sample name (e.g. SNP018ZXP2 -> 18) 
    meta["couple"] = meta["sample"].str.extract(r"^SNP(\d+)")[0]
    return meta[["gsm", "sample", "couple", "donor_karyotype", "donor_class",
                 "reported_result", "reported_class"]]


def read_genotype_calls(series_matrix: Path) -> pd.DataFrame:
    """The authors' genotype calls (GTYPE 4.1): rows = SNP ids, columns = GSM ids,
    values 0 = AA, 1 = AB, 2 = BB, -1 = no call."""
    # Local dictoinary to convert the authors' string calls to integers for storage in int8.
    GENOTYPE_CODES = {"AA": 0, "AB": 1, "BB": 2, "No Call": -1}
    # Skip all lines starting with ! to see the table 
    gt = pd.read_csv(series_matrix, sep="\t", comment="!", index_col=0, low_memory=False)
    gt = gt.replace(GENOTYPE_CODES).astype("int8")
    gt.index.name = "probe"
    return gt

# step 03
def load_probe_intensities(layout: Layout):
    """Normalised log2 probe intensities written by scripts/R/02_normalize_cel.R.

    Returns (x, probes, gsm):
      x       float32 array, probes x arrays
      probes  table with fid, fsetid (SNP id), allele (0 = A, 1 = B), strand
      gsm     GSM id of each array (column order of x)
    """
    probes = pd.read_csv(layout.probe_table, sep="\t")
    cels = layout.cel_order.read_text().split()
    x = np.fromfile(layout.probe_intensities, dtype=np.float32).reshape(len(cels), len(probes)).T
    gsm = [c.split("_")[0] for c in cels]
    return np.ascontiguousarray(x), probes, gsm

# step 4
def load_snp_annotation(layout: Layout) -> pd.DataFrame:
    """SNPs on chromosomes 1-22 and X, sorted by chromosome and position.

    Columns: fsetid (internal id), probe (SNP name, e.g. SNP_A-1780419), chr (X = 23),
    pos (hg19 bp), cytoband, frag_len (PCR fragment length), gc (probe GC fraction).
    Every SNP-level matrix in this pipeline uses this row order."""
    ann = pd.read_csv(layout.snp_annotation, sep="\t", dtype={"chrom": str})
    ann = ann[ann.chrom.notna()].copy()                       # drop control probes
    ann["chr"] = ann.chrom.replace({"X": "23"}).astype(int)
    ann = ann.rename(columns={"man_fsetid": "probe", "physical_pos": "pos"})
    ann = ann.sort_values(["chr", "pos"]).reset_index(drop=True)
    return ann[["fsetid", "probe", "chr", "pos", "cytoband", "frag_len", "gc"]]


def sample_names(layout: Layout) -> dict:
    """GSM id -> sample name (e.g. GSM1095396 -> SNP018ZXP2)."""
    meta = pd.read_csv(layout.metadata, sep="\t")
    return dict(zip(meta.gsm, meta["sample"]))


def load_signals(layout: Layout, snps: pd.DataFrame, kinds=("total", "A", "B")) -> dict:
    """Per-SNP log2 signals from step 03, rows in `snps` order, columns renamed
    from GSM ids to sample names."""
    names = sample_names(layout)
    out = {}
    for kind in kinds:
        s = pd.read_parquet(layout.signal(kind)).loc[snps.fsetid.values]
        out[kind] = s.rename(columns=names).reset_index(drop=True)
    return out

def save_snp_matrix(values: np.ndarray, snps: pd.DataFrame, samples: list, path: Path) -> None:
    """Save a SNP x array matrix (LRR or BAF) with the annotation columns in front.

    Written twice: `path` as parquet (fast, small, for Python) and the same name
    with .csv (readable anywhere: R, Excel; 4 decimals is far below the noise)."""
    ANNOTATION_COLUMNS = ["probe", "chr", "pos", "cytoband"]
    table = pd.concat([snps[ANNOTATION_COLUMNS].reset_index(drop=True),
                       pd.DataFrame(values.astype(np.float32), columns=samples)], axis=1)
    table.to_parquet(path)
    table.to_csv(path.with_suffix(".csv"), index=False, float_format="%.4f")



