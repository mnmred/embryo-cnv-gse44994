"""Log R ratio (LRR): how much DNA is at this SNP compared with a normal genome?

    LRR = log2( signal in this array / signal of a normal reference )

0 means the usual two copies. One copy lost should give log2(1/2) = -1 and one
gained log2(3/2) = +0.58 in theory, but arrays compress the scale

Steps (in the order step 04 calls them):
    reference_signal   per-SNP "normal" level = median over all arrays
    log_r_ratio        sample minus reference
    regional_gc        GC-content covariate for the wave correction
    correct_waves      remove amplification bias linked to GC / fragment length
    one_copy_shift     k, the LRR change for one lost copy
    compression_alpha  turn k into a scale factor alpha
    expected_lrr       expected LRR for any copy number (uses alpha)

Quality control - how noisy is one array's LRR?
    mapd               median absolute difference between neighbouring SNPs: SNP-level
                       noise, insensitive to real copy-number changes (rare steps)
    mapd_to_sd         MAPD -> standard deviation
    wave_sd            size of slow technical "waves" (from DNA amplification) at
                       megabase scale, with the SNP-level noise taken out
"""
import numpy as np
import pandas as pd


def reference_signal(total_log2: np.ndarray, is_x: np.ndarray, two_x: np.ndarray) -> np.ndarray:
    """Per-SNP reference = median signal over all arrays. Aneuploidies are rare
    and differ between embryos, so the median behaves like a normal genome.
    On X only arrays with two X copies are used."""
    ref = np.median(total_log2, 1)
    ref[is_x] = np.median(total_log2[is_x][:, two_x], 1)
    return ref


def log_r_ratio(total_log2: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """LRR before correction (signals are already log2)."""
    return total_log2 - reference[:, None]


def regional_gc(probe_gc: np.ndarray, chrom: np.ndarray, pos: np.ndarray,
                half_window: int = 500_000) -> np.ndarray:
    """Average probe GC fraction of all SNPs within +-half_window bp: a stand-in
    for the GC content of the surrounding genome (no genome sequence needed)."""
    out = np.empty(len(probe_gc))
    for c in np.unique(chrom):
        idx = np.flatnonzero(chrom == c); p = pos[idx]
        csum = np.r_[0, np.cumsum(probe_gc[idx])]
        lo = np.searchsorted(p, p - half_window); hi = np.searchsorted(p, p + half_window, "right")
        out[idx] = (csum[hi] - csum[lo]) / (hi - lo)
    return out


def correct_waves(lrr: np.ndarray, covariates: list, autosomal: np.ndarray,
                  n_bins: int = 25, rounds: int = 2) -> np.ndarray:
    """Remove technical bias per array: for each covariate (regional GC, PCR
    fragment length, probe GC) SNPs are split into n_bins quantile bins and each
    bin's median LRR (estimated on autosomes) is subtracted. Finally every array
    is centred so its autosomal median is 0."""
    bins = [pd.qcut(v, n_bins, labels=False, duplicates="drop") for v in covariates]
    out = lrr.copy()
    for j in range(out.shape[1]):
        v = out[:, j]
        for _ in range(rounds):
            for b in bins:
                v = v - pd.Series(v[autosomal]).groupby(b[autosomal]).median().reindex(b).values
        out[:, j] = v - np.median(v[autosomal])
    return out


def one_copy_shift(lrr: np.ndarray, x_nonpar: np.ndarray, n_x: np.ndarray) -> float:
    """k = how far LRR drops for one lost copy, measured on X of one-X arrays
    (the reference on X is built from two-X arrays). Uses no embryo labels."""
    return float(-np.median(np.median(lrr[x_nonpar], 0)[n_x == 1]))


def compression_alpha(k: float) -> float:
    """Arrays respond to copy number on a compressed scale. Model:
        observed ratio = 1 + alpha * (copies / 2 - 1)
    One lost copy gives log2(1 - alpha / 2) = -k, hence alpha = 2 (1 - 2^-k)."""
    return 2 * (1 - 2 ** -k)


def expected_lrr(copies, alpha: float):
    """Expected LRR for a copy number (scalar or array) under the compression model."""
    return np.log2(1 + alpha * (np.asarray(copies, float) / 2 - 1))


def copies_from_lrr(lrr_value, alpha: float):
    """Inverse of expected_lrr: estimated copy number from a (segment) mean LRR."""
    return 2 * (1 + (2 ** np.asarray(lrr_value, float) - 1) / alpha)


# --------------------------------------------------------------------- quality control

def mapd(values: np.ndarray) -> float:
    return float(np.median(np.abs(np.diff(values))))


def mapd_to_sd(mapd_value):
    """SD of Gaussian noise with this MAPD (neighbour differences have SD sqrt(2) sd)."""
    return mapd_value / (np.sqrt(2) * 0.6745)


def wave_sd(values: np.ndarray, chrom: np.ndarray, pos: np.ndarray, arm: np.ndarray,
            bin_bp: int = 1_000_000, min_snps: int = 20) -> float:
    """SD of the LRR waves between 1-Mb bins, without the SNP noise.
    """
    auto = chrom <= 22
    s = pd.Series(values[auto])
    s = s - s.groupby(arm[auto]).transform("median")           # ignore whole-arm copy changes
    bins = chrom[auto] * 10_000 + pos[auto] // bin_bp          # one id per chromosome and Mb
    n_snps = pd.Series(bins).value_counts()                    # SNPs per bin
    rank = s.groupby(bins).cumcount().values                   # 0, 1, 2 ... along each bin
    half = rank * 2 // n_snps.reindex(bins).values             # 0 = first half of the bin, 1 = second
    med = s.groupby([bins, half]).median().unstack()           # rows = bins, columns = half 0 / half 1
    med = med[n_snps.reindex(med.index).values >= min_snps]    # skip near-empty bins
    wave_variance = np.cov(med[0], med[1])[0, 1]
    return float(np.sqrt(max(wave_variance, 0.0)))
