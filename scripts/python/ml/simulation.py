"""Make training examples with a known answer.

Real labelled segments are too rare, so examples are made from the data itself:
take a stretch of a chromosome that is normal (2 copies), and move its LRR to where it
would be with another copy number. The real noise and waves of that embryo stay in.

    one_copy_shifts      how far LRR drops for one lost copy, measured per male on X
    lrr_shift            LRR change for any copy number (same formula as snparray/lrr.py)
    copies_from_lrr      the reverse: copy number that fits an LRR level
    normal_chromosomes   the chromosomes of each embryo that may serve as normal examples
    simulate             the examples: one row of features per simulated segment, plus the truth
"""
import numpy as np
import pandas as pd

import inputs
import settings
import tables
from features import segment_features


def one_copy_shifts(data) -> dict:
    """{sample: k} for every male. k = minus the median LRR of its X (outside the regions
    that are also on Y): males have one X, so this is the LRR of a one-copy loss."""
    on_x = tables.x_only(data)
    return {sample: float(-np.median(data.lrr[on_x, j]))
            for j, sample in enumerate(data.samples) if data.n_x[j] == 1}


def lrr_shift(copies, k: float):
    """Expected LRR of a region with `copies` copies. Arrays compress the scale:
    observed ratio = 1 + alpha * (copies / 2 - 1), with alpha set so that 1 copy gives -k."""
    alpha = 2 * (1 - 2 ** -k)
    return np.log2(1 + alpha * (np.asarray(copies, float) / 2 - 1))


def copies_from_lrr(lrr_value, k: float):
    alpha = 2 * (1 - 2 ** -k)
    return 2 * (1 + (2 ** np.asarray(lrr_value, float) - 1) / alpha)


def normal_chromosomes(data) -> dict:
    """{sample: [(chromosome, first SNP row, last SNP row + 1), ...]} of chromosomes with
    two normal copies: not listed as abnormal, and X only in samples with two X."""
    edges = np.r_[0, np.flatnonzero(np.diff(data.chrom)) + 1, len(data.chrom)]
    out = {}
    for sample, n_x in zip(data.samples, data.n_x):
        out[sample] = [(int(data.chrom[a]), int(a), int(b)) for a, b in zip(edges[:-1], edges[1:])
                       if (sample, int(data.chrom[a])) not in data.known_abnormal
                       and not (data.chrom[a] == inputs.X_CHROMOSOME and n_x != 2)
                       and b - a >= settings.MIN_SEGMENT_SNPS]
    return out


def _draw_state(rng) -> tuple:
    """A state and a true copy number for one simulated segment."""
    other = [s for s in settings.STATES if s != "normal"]
    if rng.random() < settings.SHARE_NORMAL:
        state = "normal"
    else:
        state = other[rng.integers(len(other))]
    low, high = settings.STATES[state]["simulate"]
    copies = rng.uniform(low, high)
    if state == "gain" and rng.random() < settings.SHARE_GAIN_WITH_4_COPIES:
        copies = 4.0
    return state, copies


def simulate(data, noise: pd.DataFrame, samples: list, shifts: list, per_embryo: int, rng) -> pd.DataFrame:
    """`per_embryo` simulated segments for every sample in `samples`.
    shifts: one-copy shifts (k) to draw from, one per simulated segment."""
    normal = normal_chromosomes(data)
    rows = []
    for sample in samples:
        j = data.samples.index(sample)
        chromosomes = normal[sample]
        if not chromosomes:
            continue
        lengths = np.array([b - a for _, a, b in chromosomes])
        sample_noise = noise.loc[sample].to_dict()
        for _ in range(per_embryo):
            # where: a random stretch of a normal chromosome (longer chromosomes more often)
            chromosome, a, b = chromosomes[rng.choice(len(chromosomes), p=lengths / lengths.sum())]
            if rng.random() < settings.SHARE_WHOLE_CHROMOSOME:
                n_snps = b - a
            else:                                   # any length, short ones as likely as long ones
                n_snps = int(np.exp(rng.uniform(np.log(settings.MIN_SEGMENT_SNPS), np.log(b - a))))
            start = a + rng.integers(b - a - n_snps + 1)
            # what: a state, and the LRR change that goes with it
            state, copies = _draw_state(rng)
            shift = float(lrr_shift(copies, rng.choice(shifts)))
            values = data.lrr[start:start + n_snps, j] + shift
            rows.append(dict(sample=sample, chr=chromosome, row_start=int(start), row_end=int(start + n_snps),
                             n_snps=n_snps, true_state=state, true_copies=copies, lrr_shift=shift,
                             **segment_features(values, sample_noise)))
    return pd.DataFrame(rows)
