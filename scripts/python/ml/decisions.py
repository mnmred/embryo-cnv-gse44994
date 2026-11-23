"""From the model's probabilities to final calls.

    decide            probability of every state -> one final state per segment
    add_copies        final state -> copy number, and whether that is abnormal
    join_neighbours   abnormal segments next to each other with the same state -> one region
"""
import numpy as np
import pandas as pd

import settings

STATE_NAMES = list(settings.STATES)


def decide(p: pd.DataFrame, normal_state: pd.Series, min_probability: float = settings.MIN_PROBABILITY) -> pd.DataFrame:
    """Final state of every row. A segment keeps its normal state unless the model gives
    the other states together at least `min_probability`; then the most probable of them wins."""
    p_normal_state = p.to_numpy()[np.arange(len(p)), [STATE_NAMES.index(s) for s in normal_state]]
    p_abnormal = 1 - p_normal_state
    others = p.copy()
    for state in STATE_NAMES:
        others.loc[normal_state == state, state] = -1.0          # the normal state cannot win here
    state = np.where(p_abnormal >= min_probability, others.idxmax(axis=1), normal_state)
    return pd.DataFrame(dict(state=state, probability_abnormal=p_abnormal), index=p.index)


def add_copies(table: pd.DataFrame) -> pd.DataFrame:
    """Adds copies, normal_copies and abnormal to segments that have a state, a normal_state
    and copies_from_lrr. A gain counts as 4 copies when the LRR level fits 3.5 copies or more."""
    table["copies"] = table["state"].map(lambda s: settings.STATES[s]["copies"])
    table.loc[(table["state"] == "gain") & (table["copies_from_lrr"] >= 3.5), "copies"] = 4.0
    table["normal_copies"] = table["normal_state"].map(lambda s: settings.STATES[s]["copies"])
    table["abnormal"] = table["state"] != table["normal_state"]
    return table


def join_neighbours(called: pd.DataFrame) -> pd.DataFrame:
    """Abnormal segments only; neighbours on one chromosome with the same state become one row."""
    rows = []
    for (sample, chromosome), group in called.groupby(["sample", "chr"], sort=False):
        previous_state = None
        for seg in group.sort_values("row_start").itertuples():
            if seg.state == previous_state and seg.abnormal:
                last = rows[-1]
                last.update(end=seg.end, n_snps=last["n_snps"] + seg.n_snps, n_pieces=last["n_pieces"] + 1,
                            probability_abnormal=min(last["probability_abnormal"], seg.probability_abnormal),
                            copies=max(last["copies"], seg.copies) if "gain" in seg.state else last["copies"])
            elif seg.abnormal:
                rows.append(dict(sample=sample, chr=chromosome, start=seg.start, end=seg.end, n_snps=seg.n_snps,
                                 n_pieces=1, state=seg.state, copies=seg.copies,
                                 normal_copies=seg.normal_copies,
                                 probability_abnormal=seg.probability_abnormal))
            previous_state = seg.state
    calls = pd.DataFrame(rows)
    calls.insert(4, "size_mb", (calls["end"] - calls["start"]) / 1e6)
    calls["direction"] = np.where(calls["copies"] > calls["normal_copies"], "gain", "loss")
    return calls
