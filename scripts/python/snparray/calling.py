"""Copy numbers for segments cut by an R package (step 06b: copynumber::aspcf).

The package only cuts each chromosome into pieces and reports each piece's mean LRR;
it does not say how many copies a piece has. This module does that. Nothing here is
specific to ASPCF, so segments of another method (e.g. DNAcopy's CBS) could be
called the same way.

    segment_rows    segment bp coordinates -> SNP rows of lrr.parquet
    bin_noise       per array: scatter of 50-SNP averages (sigma)
    between_arm_sd  per array: scatter of whole-arm averages (tau)
    call_segments   join neighbouring pieces that change copy number the same way,
                    then test every joined piece and give it a copy number

A piece is reported when all three hold:
    |z| >= 3            z = (median LRR - normal LRR) / sqrt(sigma^2 * 50 / n_snps + tau^2)
                        sigma: noise of 50-SNP averages; tau: amplification moves whole arms
                        a little, which does not average away over many SNPs
    |change| >= 0.4     copies gained (+) or lost (-) per cell, from the median LRR
                        (lrr.copies_from_lrr); >= 0.8 is "full", 0.4-0.8 "mosaic"
    size >= 2.36 Mb     for full calls (the authors' minimum); 10 Mb for mosaic calls
Pieces with |change| < 0.4, and pieces that fail a test, keep the normal copy number.
The normal copy number is 2, except on X in a one-X (male) array: 1, but 2 in the
pseudo-autosomal regions at the two ends of X, which the Y chromosome also carries.
"""
import numpy as np
import pandas as pd

from .lrr import copies_from_lrr, expected_lrr
from .sex import x_not_on_y_mask

BIN_SIZE = 50                                   # SNPs per average for sigma
MIN_ARM_BINS = 10                               # arms shorter than this are not used for tau
MIN_Z, MIN_CHANGE, FULL_CHANGE = 3.0, 0.4, 0.8
MIN_SIZE_FULL, MIN_SIZE_MOSAIC = 2.36e6, 10e6   # bp


def normal_copies(chrom: int, n_x: int, pos: np.ndarray) -> int:
    """Normal copy number of a piece with SNP positions `pos` on chromosome `chrom`:
    1 for X in a one-X (male) array, unless most of the piece lies in a pseudo-autosomal
    region (also on Y, so 2 copies); otherwise 2."""
    if chrom == 23 and n_x == 1:
        on_y_too = ~x_not_on_y_mask(np.full(len(pos), 23), pos)
        return 2 if on_y_too.mean() >= 0.5 else 1
    return 2


def segment_rows(segments: pd.DataFrame, chrom: np.ndarray, pos: np.ndarray) -> pd.DataFrame:
    """Add row_start / row_end (SNP rows; end exclusive) to segments with chr, start, end
    (bp of their first and last SNP, as the R packages report them)."""
    row = pd.Series(np.arange(len(chrom)), index=pd.MultiIndex.from_arrays([chrom, pos]))
    i = row.reindex(pd.MultiIndex.from_arrays([segments.chr, segments.start])).values
    j = row.reindex(pd.MultiIndex.from_arrays([segments.chr, segments.end])).values
    if np.isnan(i).any() or np.isnan(j).any():
        raise ValueError("segment ends do not match SNP positions in lrr.parquet")
    return segments.assign(row_start=i.astype(int), row_end=j.astype(int) + 1)


def _bins(chrom: np.ndarray, size: int):
    """Bin number of every SNP: `size` consecutive SNPs of one chromosome (the last bin of a
    chromosome takes the leftover SNPs), plus the first SNP row of every bin."""
    edges = np.r_[0, np.flatnonzero(np.diff(chrom)) + 1, len(chrom)]
    bin_of, first, n = np.empty(len(chrom), np.int64), [], 0
    for a, e in zip(edges[:-1], edges[1:]):
        k = max(1, (e - a) // size)
        bin_of[a:e] = n + np.minimum(np.arange(e - a) // size, k - 1)
        first.append(a + np.arange(k) * size)
        n += k
    return bin_of, np.concatenate(first)


def bin_noise(y: np.ndarray, chrom: np.ndarray, arm: np.ndarray, size: int = BIN_SIZE):
    """sigma: robust SD of the array's autosomal 50-SNP averages around their arm median.
    Also returns the averages and their arms, for between_arm_sd."""
    bin_of, first = _bins(chrom, size)
    means = np.bincount(bin_of, y) / np.bincount(bin_of)
    bin_arm, auto = arm[first], chrom[first] <= 22
    resid = means[auto] - pd.Series(means[auto]).groupby(bin_arm[auto]).transform("median").values
    return 1.4826 * np.median(np.abs(resid - np.median(resid))), means[auto], bin_arm[auto]


def between_arm_sd(bin_means: np.ndarray, bin_arm: np.ndarray, alpha: float) -> float:
    """tau: robust SD of the mean LRR of autosomal arms (>= 10 bins) whose mean does not
    itself show a change of >= 0.4 copies."""
    arms = pd.Series(bin_means).groupby(bin_arm).agg(["mean", "size"])
    arms = arms[arms["size"] >= MIN_ARM_BINS]
    normal = arms.loc[np.abs(copies_from_lrr(arms["mean"], alpha) - 2) < MIN_CHANGE, "mean"]
    return 1.4826 * np.median(np.abs(normal - np.median(normal)))


def call_segments(segments: pd.DataFrame, y: np.ndarray, pos: np.ndarray, cytoband: np.ndarray,
                  sigma: float, tau: float, alpha: float, n_x: int) -> pd.DataFrame:
    """segments: one array's pieces with chr, row_start, row_end; y: its (clipped) LRR.

    Neighbouring pieces on one chromosome are joined while they change copy number the
    same way (gain, loss, or less than 0.4 copies either way). Every joined piece is
    tested (see the module docstring). One row per joined piece:
    chr, start, end, size_mb, cytoband, n_snps, n_pieces, row_start, row_end,
    mean_lrr (holds the median), change, z, reported, copies (final), kind (normal / full / mosaic)."""
    def change_of(c, i, j):
        return float(copies_from_lrr(np.median(y[i:j]), alpha)) - normal_copies(c, n_x, pos[i:j])

    runs = []                                   # [chr, row_start, row_end, direction, n_pieces]
    for p in segments.sort_values("row_start").itertuples():
        d = change_of(p.chr, p.row_start, p.row_end)
        direction = int(np.sign(d)) if abs(d) >= MIN_CHANGE else 0
        if runs and runs[-1][0] == p.chr and runs[-1][3] == direction:
            runs[-1][2] = p.row_end
            runs[-1][4] += 1
        else:
            runs.append([p.chr, p.row_start, p.row_end, direction, 1])

    rows = []
    for c, i, j, direction, n_pieces in runs:
        normal = normal_copies(c, n_x, pos[i:j])
        level = float(np.median(y[i:j]))
        change = change_of(c, i, j)
        z = (level - float(expected_lrr(normal, alpha))) / np.sqrt(sigma ** 2 * BIN_SIZE / (j - i) + tau ** 2)
        size = pos[j - 1] - pos[i]
        kind = "full" if abs(change) >= FULL_CHANGE else "mosaic"
        reported = bool(direction != 0 and abs(z) >= MIN_Z and abs(change) >= MIN_CHANGE
                        and size >= (MIN_SIZE_MOSAIC if kind == "mosaic" else MIN_SIZE_FULL))
        if reported:
            copies = normal + np.sign(change) * (max(1, round(abs(change))) if kind == "full" else 0.5)
        else:
            copies, kind = normal, "normal"
        rows.append(dict(chr=int(c), start=int(pos[i]), end=int(pos[j - 1]), size_mb=size / 1e6,
                         cytoband=f"{cytoband[i]}-{cytoband[j - 1]}", n_snps=int(j - i), n_pieces=n_pieces,
                         row_start=int(i), row_end=int(j), mean_lrr=level, change=change, z=float(z),
                         reported=reported, copies=float(copies), kind=kind))
    return pd.DataFrame(rows)
