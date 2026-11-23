"""Read the tables named in inputs.py into one simple form (class Data).

    read_table   any .parquet / .csv / tab-separated file -> pandas table
    load         everything the ML caller needs, checked and with standard column names
    x_only       which SNPs lie on X outside the regions that are also on Y
"""
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

import inputs

REPO = Path(__file__).resolve().parents[3]          # scripts/python/ml/tables.py -> repo


def read_table(path, columns=None) -> pd.DataFrame:
    """Read a table; the format follows from the file name (see inputs.py)."""
    path = REPO / path
    if not path.exists():
        raise FileNotFoundError(f"{path} does not exist; check scripts/python/ml/inputs.py")
    name = path.name.lower().removesuffix(".gz")
    if name.endswith(".parquet"):
        return pd.read_parquet(path, columns=columns)
    table = pd.read_csv(path, sep="," if name.endswith(".csv") else "\t")
    return table if columns is None else table[columns]


def output_path(file_name: str) -> Path:
    folder = REPO / inputs.OUTPUT_FOLDER
    folder.mkdir(parents=True, exist_ok=True)
    return folder / file_name


@dataclass
class Data:
    samples: list             # sample names, in the order of the LRR columns
    n_x: np.ndarray           # X copies of every sample (1 or 2)
    chrom: np.ndarray         # chromosome of every SNP
    pos: np.ndarray           # position (bp) of every SNP
    lrr: np.ndarray           # SNPs x samples
    segments: pd.DataFrame    # sample, chr, start, end, row_start, row_end (SNP rows; end exclusive)
    known_abnormal: set       # (sample, chromosome) pairs that are not normal examples


def load() -> Data:
    sample_table = read_table(inputs.SAMPLES_TABLE)
    samples = list(sample_table[inputs.SAMPLE_NAME_COLUMN].astype(str))
    n_x = sample_table[inputs.SAMPLE_X_COPIES_COLUMN].to_numpy(int)

    lrr_table = read_table(inputs.LRR_TABLE,
                           [inputs.LRR_CHROMOSOME_COLUMN, inputs.LRR_POSITION_COLUMN] + samples)
    chrom = lrr_table[inputs.LRR_CHROMOSOME_COLUMN].to_numpy(int)
    pos = lrr_table[inputs.LRR_POSITION_COLUMN].to_numpy(np.int64)
    lrr = lrr_table[samples].to_numpy(np.float32, copy=True)     # own copy: it is clipped later
    genome_order = chrom * 10_000_000_000 + pos     # one increasing number per SNP
    if (np.diff(genome_order) < 0).any():
        raise ValueError("the LRR table must be sorted by chromosome and position")
    if np.isnan(lrr).any():
        raise ValueError("the LRR table has missing values; remove those SNPs first")

    segments = read_table(inputs.SEGMENTS_TABLE).rename(columns={
        inputs.SEGMENT_SAMPLE_COLUMN: "sample", inputs.SEGMENT_CHROMOSOME_COLUMN: "chr",
        inputs.SEGMENT_START_COLUMN: "start", inputs.SEGMENT_END_COLUMN: "end"})
    segments = segments[["sample", "chr", "start", "end"]].copy()
    segments["sample"] = segments["sample"].astype(str)
    unknown = set(segments["sample"]) - set(samples)
    if unknown:
        raise ValueError(f"segments of samples that are not in the samples table: {sorted(unknown)}")
    # SNP rows of every segment: first SNP at or after `start`, last SNP at or before `end`
    chr_of_segment = segments["chr"].to_numpy(int) * 10_000_000_000
    segments["row_start"] = np.searchsorted(genome_order, chr_of_segment + segments["start"].to_numpy(np.int64), "left")
    segments["row_end"] = np.searchsorted(genome_order, chr_of_segment + segments["end"].to_numpy(np.int64), "right")

    known_abnormal = set()
    for path in inputs.KNOWN_ABNORMAL_TABLES:
        table = read_table(path)
        known_abnormal |= set(zip(table[inputs.SEGMENT_SAMPLE_COLUMN].astype(str),
                                  table[inputs.SEGMENT_CHROMOSOME_COLUMN].astype(int)))

    return Data(samples, n_x, chrom, pos, lrr, segments, known_abnormal)


def x_only(data: Data) -> np.ndarray:
    """True for SNPs on X outside the pseudo-autosomal regions (one copy in males)."""
    return ((data.chrom == inputs.X_CHROMOSOME)
            & (data.pos > inputs.X_PAR1_END) & (data.pos < inputs.X_PAR2_START))
