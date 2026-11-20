"""Many probes -> one signal per SNP.

Probes differ in how strongly theybind, so their raw brightness is not directly comparable. 
Median polish (the summarisation step of RMA / oligo::snprma) fits

    log2 intensity[probe, array] = overall + probe effect + array effect + noise

robustly (with medians), and the array effect is the SNP's signal in that array.

    summarize_snps(x, probes) -> {"total": all probes,
                                  "A": allele-A probes, "B": allele-B probes}
"""
import numpy as np
import pandas as pd
from numba import njit, prange


@njit(cache=True) #make it C so it runs faster, and cache it so it doesn't have to be recompiled every time
def median_polish(z, max_iter=10, eps=0.01):
    """Tukey median polish of one probes x arrays block (log2 scale).
    Returns one value per array (overall level + array effect). Modifies z."""
    n_rows, n_cols = z.shape
    row = np.zeros(n_rows); col = np.zeros(n_cols); overall = 0.0; old = 0.0
    for _ in range(max_iter):
        # median of each row subtracted, row=probe
        for i in range(n_rows):
            d = np.median(z[i, :]); z[i, :] -= d; row[i] += d
        d = np.median(col); col -= d; overall += d
        # median of each column subtracted, column=embryo
        for j in range(n_cols):
            d = np.median(z[:, j]); z[:, j] -= d; col[j] += d
        # keep snp real brightness
        d = np.median(row); row -= d; overall += d
        new = np.abs(z).sum()
        if new == 0 or abs(new - old) < eps * new:
            break
        old = new
    return overall + col


@njit(parallel=True, cache=True)
def _polish_groups(x, order, starts):
    n_groups = starts.size - 1
    out = np.empty((n_groups, x.shape[1]))
    for g in prange(n_groups):
        idx = order[starts[g]:starts[g + 1]]
        out[g] = median_polish(x[idx].astype(np.float64))
    return out


def summarize_groups(x: np.ndarray, keys: list):
    """Median-polish every group of probe rows that share the same key(s).
    Returns (values: groups x arrays, group_keys: groups x len(keys))."""
    order = np.lexsort(keys[::-1]).astype(np.int64) # sort
    k = np.stack(keys, 1)[order]
    starts = np.r_[0, np.flatnonzero((k[1:] != k[:-1]).any(1)) + 1, len(order)] # and record where each group starts
    return _polish_groups(x, order, starts), k[starts[:-1]]


def summarize_snps(x: np.ndarray, probes: pd.DataFrame, columns: list) -> dict:
    """Per-SNP log2 signals: 'total' (every probe of the SNP, used for LRR) and
    'A' / 'B' (allele-specific probes, both strands pooled, used for BAF).
    Each value is a DataFrame indexed by fsetid with one column per array."""
    out = {}
    values, keys = summarize_groups(x, [probes.fsetid.values]) #group by SNP
    out["total"] = pd.DataFrame(values, index=keys[:, 0], columns=columns)
    values, keys = summarize_groups(x, [probes.allele.values, probes.fsetid.values]) #group by allele and SNP
    for allele, name in [(0, "A"), (1, "B")]: #split results into allele A and B
        m = keys[:, 0] == allele
        out[name] = pd.DataFrame(values[m], index=keys[m, 1], columns=columns)
    for df in out.values():
        df.index.name = "fsetid"
    return out
