"""The numbers the model sees: a few summaries of the LRR of one segment.

The model never sees single SNPs. Every segment (real or simulated) becomes one row of
FEATURES, computed by segment_features from the segment's LRR values and from the
noise of the sample it belongs to (sample_noise).
"""
import numpy as np
import pandas as pd

import inputs
import settings

FEATURES = [
    "lrr_median",                  # level of the segment: 0 = two copies, lower = loss, higher = gain
    "lrr_median_in_noise_units",   # the same, divided by the noise expected for a segment this long
    "halves_difference",           # do the two halves of the segment agree on the level?
    "log10_snps",                  # length of the segment
    "segment_bin_noise",           # scatter of 50-SNP averages inside the segment
    "sample_snp_noise",            # the sample's noise between neighbouring SNPs
    "sample_bin_noise",            # the sample's scatter of 50-SNP averages (waves)
    "sample_chromosome_noise",     # the sample's scatter of whole-chromosome levels
]


def robust_sd(values: np.ndarray) -> float:
    """Standard deviation estimated from the median absolute deviation (ignores outliers)."""
    return float(1.4826 * np.median(np.abs(values - np.median(values))))


def bin_means(values: np.ndarray) -> np.ndarray:
    """Average of every 50 consecutive SNPs (leftover SNPs at the end are dropped)."""
    n_bins = len(values) // settings.BIN_SNPS
    return values[:n_bins * settings.BIN_SNPS].reshape(n_bins, settings.BIN_SNPS).mean(1)


def sample_noise(data) -> pd.DataFrame:
    """Three noise levels per sample, measured on the autosomes, and the LRR clipped in place:
    single SNPs further than 4 SD from 0 are pulled in, so they cannot dominate an average."""
    rows = []
    autosomes = [c for c in np.unique(data.chrom) if c != inputs.X_CHROMOSOME]
    for j, sample in enumerate(data.samples):
        column = data.lrr[:, j]
        autosomal = column[data.chrom != inputs.X_CHROMOSOME]
        snp_noise = float(np.median(np.abs(np.diff(autosomal))) / (np.sqrt(2) * 0.6745))
        np.clip(column, -settings.CLIP_AT_SD * snp_noise, settings.CLIP_AT_SD * snp_noise, out=column)

        bin_residuals, chromosome_levels = [], []
        for c in autosomes:
            bins = bin_means(column[data.chrom == c])
            bin_residuals.append(bins - np.median(bins))
            chromosome_levels.append(np.median(bins))
        rows.append(dict(sample=sample, sample_snp_noise=snp_noise,
                         sample_bin_noise=robust_sd(np.concatenate(bin_residuals)),
                         sample_chromosome_noise=robust_sd(np.array(chromosome_levels))))
    return pd.DataFrame(rows).set_index("sample")



def segment_features(values: np.ndarray, noise: dict) -> dict:
    """values: LRR of the SNPs of one segment (at least MIN_SEGMENT_SNPS of them);
    noise: that sample's row of sample_noise()."""
    n = len(values)
    level = float(np.median(values))
    bins = bin_means(values)
    expected_noise = noise["sample_bin_noise"] / np.sqrt(n / settings.BIN_SNPS)
    return dict(
        lrr_median=level,
        lrr_median_in_noise_units=level / expected_noise,
        halves_difference=abs(float(np.median(values[:n // 2]) - np.median(values[n // 2:]))),
        log10_snps=np.log10(n),
        segment_bin_noise=robust_sd(bins) if len(bins) >= 4 else noise["sample_bin_noise"],
        sample_snp_noise=noise["sample_snp_noise"],
        sample_bin_noise=noise["sample_bin_noise"],
        sample_chromosome_noise=noise["sample_chromosome_noise"],
    )
