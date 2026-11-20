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
