"""Score a list of calls against reported results.

A reported region is "found" when a call of the same sample, on the same chromosome and
in the same direction (gain / loss) shares at least 1 bp with it; a call is "confirmed"
when it shares at least 1 bp with such a reported region. The strict version asks for a
reciprocal overlap of at least 50 %.

Both tables need the columns: sample, chr, start, end, direction ("gain" / "loss").
"""
import numpy as np
import pandas as pd

MIN_RECIPROCAL_OVERLAP = 0.5


def best_overlap(these: pd.DataFrame, those: pd.DataFrame) -> np.ndarray:
    """For every row of `these`: the best reciprocal overlap (0-1) with a row of `those`
    of the same sample, chromosome and direction."""
    out = np.zeros(len(these))
    for i, row in enumerate(these.itertuples()):
        same = those[(those["sample"] == row.sample) & (those["chr"] == row.chr)
                     & (those["direction"] == row.direction)]
        if len(same):
            shared = np.minimum(row.end, same["end"]) - np.maximum(row.start, same["start"]) + 1
            longest = np.maximum(row.end - row.start, same["end"] - same["start"]) + 1
            out[i] = max(0.0, float((shared / longest).max()))
    return out


def score(calls: pd.DataFrame, reported: pd.DataFrame, all_samples: list):
    """Returns (metrics, reported + overlap / found, calls + overlap / confirmed)."""
    reported, calls = reported.copy(), calls.copy()
    reported["overlap"] = best_overlap(reported, calls)
    reported["found"] = reported["overlap"] > 0
    calls["overlap"] = best_overlap(calls, reported)
    calls["confirmed"] = calls["overlap"] > 0

    abnormal = np.isin(all_samples, reported["sample"])     # embryos with a reported change
    called = np.isin(all_samples, calls["sample"])          # embryos with at least one call
    strict = f"per region (overlap >= {MIN_RECIPROCAL_OVERLAP:.0%})"
    rows = [
        ("per embryo", "sensitivity", called[abnormal]),
        ("per embryo", "specificity", ~called[~abnormal]),
        ("per region (any overlap)", "recall", reported["found"]),
        ("per region (any overlap)", "precision", calls["confirmed"]),
        (strict, "recall", reported["overlap"] >= MIN_RECIPROCAL_OVERLAP),
        (strict, "precision", calls["overlap"] >= MIN_RECIPROCAL_OVERLAP),
    ]
    metrics = pd.DataFrame([dict(group=g, name=n, x=int(np.sum(hit)), n=len(hit)) for g, n, hit in rows])
    metrics["value"] = metrics["x"] / metrics["n"]
    return metrics, reported, calls


def state_metrics(truth: pd.Series, answer: pd.Series) -> dict:
    """How well `answer` (a state per segment) matches the simulated truth."""
    abnormal = truth != "normal"
    return {
        "exact state right": float((answer == truth).mean()),
        "abnormal segments reported": float((answer[abnormal] != "normal").mean()),
        "normal segments reported (false alarms)": float((answer[~abnormal] != "normal").mean()),
        "n": len(truth),
    }
