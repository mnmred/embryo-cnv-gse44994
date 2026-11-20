"""
Females have two X chromosomes, males one (plus a Y, which this array does not
measure). One X copy means half the X DNA, so the X signal is lower than the
autosome signal by the same amount as any other single-copy loss.

"""
import numpy as np

# hg19 pseudo-autosomal regions: present on both X and Y, so they behave like
# two copies in males too and must be left out
HG19_PAR1_END = 2_699_520
HG19_PAR2_START = 154_931_044


def x_not_on_y_mask(chrom: np.ndarray, pos: np.ndarray) -> np.ndarray:
    """True for X-chromosome SNPs outside the pseudo-autosomal regions"""
    return (chrom == 23) & (pos > HG19_PAR1_END) & (pos < HG19_PAR2_START)


def x_dosage(total_log2: np.ndarray, chrom: np.ndarray, pos: np.ndarray) -> np.ndarray:
    """Median X signal minus median autosome signal, per array, after removing
    each SNP's typical level (its median over all arrays)"""
    d = total_log2 - np.median(total_log2, 1, keepdims=True)
    return np.median(d[x_not_on_y_mask(chrom, pos)], 0) - np.median(d[chrom <= 22], 0)


def count_x_copies(total_log2: np.ndarray, chrom: np.ndarray, pos: np.ndarray) -> np.ndarray:
    """1 or 2X copies per array: two-group k-means on the X dosage. This assumes we have both
     male and female, in our casewe do"""
    x_signal = x_dosage(total_log2, chrom, pos)
    low_centre = x_signal.min()     # the most male-looking array
    high_centre = x_signal.max()    # the most female-looking array

    for _ in range(20): # iterate until convergence, 20 is safe because there are only 2 groups
        distance_to_low = np.abs(x_signal - low_centre)
        distance_to_high = np.abs(x_signal - high_centre)
        in_high_group = distance_to_high < distance_to_low # true if closer to 'high centre'
        in_low_group  = distance_to_low <= distance_to_high
        # update centres to the mean of their group before the next iteration
        low_centre = x_signal[in_low_group].mean()
        high_centre = x_signal[in_high_group].mean()

    # Bright X = two copies (female), dim X = one copy (male)
    n_x_copies = np.where(in_high_group, 2, 1)
    return n_x_copies
