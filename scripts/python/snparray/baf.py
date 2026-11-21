"""B-allele frequency (BAF): which share of the DNA at a SNP carries allele B?

At each SNP a person has two copies, each carrying allele A or B:
    AA -> BAF 0      AB -> BAF 0.5      BB -> BAF 1
With three copies (trisomy) heterozygous SNPs sit at 1/3 or 2/3; with one copy
(monosomy) no SNP can be heterozygous, so only 0 and 1 remain.

How BAF is computed (PennCNV-Affy approach):
  1. polar_coordinates   A and B signals -> angle theta (0 = only A, 1 = only B)
  2. genotype_clusters   per SNP, where do AA / AB / BB arrays sit on theta?
                         (from the authors' genotype calls; rare clusters imputed)
  3. theta_to_baf        put this array's theta on the 0 / 0.5 / 1 scale between
                         those cluster centres
compute_baf runs all three, leave-one-out: an array's own call never helps
define the clusters its BAF is measured against (with only 3-5 arrays per
cluster, its own call would otherwise pin its BAF to exactly 0.5).

Caveat for this dataset: DNA amplified from a few cells is so noisy that
single-SNP BAF barely separates heterozygous SNPs; use it averaged over regions.
"""
import warnings

import numpy as np


def polar_coordinates(a_log2: np.ndarray, b_log2: np.ndarray):
    """theta = 2/pi * atan2(B, A) in [0, 1] and log2 R = log2(A + B), from log2 signals."""
    # undo the log2 to get the original signals, then compute polar coordinates
    a, b = 2.0 ** a_log2, 2.0 ** b_log2
    return 2 / np.pi * np.arctan2(b, a), np.log2(a + b)


def cluster_centres(values: np.ndarray, genotypes: np.ndarray, use: np.ndarray, min_n: int = 3):
    """Median of `values` (SNP x array) over the arrays called AA, AB, BB."""
    centres = np.full((len(values), 3), np.nan)
    counts = np.zeros((len(values), 3), int)
    # one passs per genotype, AA, AB, BB
    for g in range(3):
        m = (genotypes == g) & use
        counts[:, g] = m.sum(1)
        with warnings.catch_warnings():      # SNPs without any such call is NaN
            warnings.simplefilter("ignore", RuntimeWarning)
            centres[:, g] = np.nanmedian(np.where(m, values, np.nan), 1)
    centres[counts < min_n] = np.nan # use min 3 to have enough for a centre defined by the median
    return centres, counts

def impute_centres(centres: np.ndarray) -> np.ndarray:
    """Fill missing cluster centres by linear regression on the observed ones,
    learnt from SNPs where all three clusters were observed."""
    filled = centres.copy()
    observed = ~np.isnan(centres)  # SNP x 3: which of AA / AB / BB has a centre
    complete = observed.all(axis=1) # SNPs with all three centres: the training set

    # group the incomplete SNPs by what they miss, for example (True, True, False) means only BB missing
    missing_patterns = {tuple(row) for row in observed[~complete]}

    # one regression per pattern of missing clusters
    for pattern in missing_patterns:
        snps_with_pattern = (observed == pattern).all(axis=1)
        known = np.flatnonzero(pattern)                  # so [0, 1] is AA and AB observed
        missing = np.flatnonzero(~np.array(pattern))     # so [2] is BB missing

        # predict each missing cluster from the observed ones
        for target in missing:
            # no cluster observed: nothing to regress on, use the typical centre
            if known.size == 0:
                filled[snps_with_pattern, target] = np.median(centres[complete, target])
                continue

            # learn on complete SNPs: target centre = b0 + b1 * known centre (+ b2 * the other)
            X_train = np.column_stack([np.ones(complete.sum()),
                                       centres[np.ix_(complete, known)]])
            y_train = centres[complete, target]
            coefficients = np.linalg.lstsq(X_train, y_train, rcond=None)[0]

            # apply the learnt line to the SNPs with this pattern
            X_new = np.column_stack([np.ones(snps_with_pattern.sum()),
                                     centres[np.ix_(snps_with_pattern, known)]])
            filled[snps_with_pattern, target] = X_new @ coefficients
    return filled

def genotype_clusters(theta: np.ndarray, genotypes: np.ndarray, use: np.ndarray,
                      min_n: int = 3, min_gap: float = 0.05):
    """AA / AB / BB centres on the theta scale for every SNP. Centres that are not
    in the order AA < AB < BB (by at least min_gap) are replaced by the typical
    centres. Returns (centres: SNP x 3, counts: SNP x 3)."""
    centres, counts = cluster_centres(theta, genotypes, use, min_n)
    centres = impute_centres(centres)
    # centres too close or out of order, baf would explode or run backwards
    bad = ~((centres[:, 1] - centres[:, 0] > min_gap) & (centres[:, 2] - centres[:, 1] > min_gap))
    centres[bad] = np.median(centres[~bad], 0)
    return centres, counts


def theta_to_baf(theta: np.ndarray, centres: np.ndarray) -> np.ndarray:
    """Piecewise-linear map: AA centre -> 0, AB centre -> 0.5, BB centre -> 1.
    Not clipped, so values slightly below 0 or above 1 occur (noise)."""
    aa, ab, bb = centres[:, 0], centres[:, 1], centres[:, 2]
    if theta.ndim == 2:
        aa, ab, bb = aa[:, None], ab[:, None], bb[:, None]
    return np.where(theta < ab, 0.5 * (theta - aa) / (ab - aa), 0.5 + 0.5 * (theta - ab) / (bb - ab))


def genotypes_from_baf(baf: np.ndarray) -> np.ndarray:
    """Simple genotype read-out: < 0.15 AA, 0.35-0.65 AB, > 0.85 BB, else -1."""
    g = np.full(baf.shape, -1, np.int8)
    g[baf < 0.15] = 0; g[(baf > 0.35) & (baf < 0.65)] = 1; g[baf > 0.85] = 2
    return g


def call_concordance(baf: np.ndarray, genotypes: np.ndarray, rows: np.ndarray) -> np.ndarray:
    """Per array: fraction of SNPs (in `rows`) where the BAF read-out agrees with
    the deposited genotype call, among SNPs where both give a genotype."""
    own = genotypes_from_baf(baf)
    both = (genotypes >= 0) & (own >= 0) & rows[:, None]
    return ((own == genotypes) & both).sum(0) / both.sum(0)


def compute_baf(theta: np.ndarray, genotypes: np.ndarray, use: np.ndarray,
                autosomal: np.ndarray, min_concordance: float = 0.8) -> dict:
    """Full BAF computation for all arrays (SNP x array inputs).

    Pass 1: clusters from every usable call -> BAF -> arrays whose deposited calls
            agree with their own CEL data less than `min_concordance` are flagged
            (their calls belong to another sample) and excluded from clustering.
    Pass 2: leave-one-out clusters for each array -> final BAF.

    Returns dict: baf, concordance (per array), mismatch (bool per array),
    use (final mask of calls used), centres / counts (clusters from all used calls).
    """
    use = use.copy()
    centres, _ = genotype_clusters(theta, genotypes, use)
    mismatch = call_concordance(theta_to_baf(theta, centres), genotypes, autosomal) < min_concordance
    use[:, mismatch] = False
    baf = np.empty_like(theta)
    # leave-one-out: for each array, compute clusters from all other arrays and convert its theta to BAF
    for j in range(theta.shape[1]):
        loo = use.copy(); loo[:, j] = False
        centres_j, _ = genotype_clusters(theta, genotypes, loo)
        baf[:, j] = theta_to_baf(theta[:, j], centres_j)
    centres, counts = genotype_clusters(theta, genotypes, use)
    return dict(baf=baf, concordance=call_concordance(baf, genotypes, autosomal),
                mismatch=mismatch, use=use, centres=centres, counts=counts)


