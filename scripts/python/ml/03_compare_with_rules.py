"""Step 3 (this repository only) - the ML caller against the rules of step 06c.

Run (repo root, ml_env active):   python scripts/python/ml/03_compare_with_rules.py
Reads:   data/ml: simulated_test.tsv, cnv_calls_ML.tsv, reported_abnormal.tsv
         data/processed (through snparray.paths, as the main pipeline does):
         lrr.parquet, arrays.tsv, call_noise.tsv, calibration.tsv, cnv_calls_ASPCF.tsv
Writes:  data/ml/comparison_with_rules.txt
         data/ml/simulated_test_rules.tsv   the rules' answer on every simulated segment

Two comparisons:
    1. the same simulated segments, called by the rules (snparray.calling.call_segments,
       imported unchanged) and by the model
    2. the real calls of both, scored against the paper's results in the same way

Unlike steps 01 and 02 this script knows the file layout of the main pipeline, so it
only works in this repository.
"""
import sys

import numpy as np
import pandas as pd

import scoring
import settings
import tables

sys.path.insert(0, str(tables.REPO / "scripts" / "python"))
from snparray import calling, paths     # noqa: E402  (the main pipeline's code)

LONG = 1000      # SNPs, as in 02_evaluate.py
lines = []


def say(text=""):
    print(text)
    lines.append(str(text))


def rules_on_simulated(test: pd.DataFrame, layout) -> pd.Series:
    """State that step 06c gives each simulated segment: the segment's real LRR (clipped as
    in 06c) plus the simulated shift goes through call_segments as one piece."""
    samples = sorted(test["sample"].unique())
    lrr_table = pd.read_parquet(layout.lrr, columns=["chr", "pos", "cytoband"] + samples)
    arrays = pd.read_csv(layout.arrays, sep="\t", index_col=0)
    noise = pd.read_csv(layout.call_noise, sep="\t", index_col=0)
    alpha = pd.read_csv(layout.calibration, sep="\t", header=None, index_col=0)[1]["alpha"]
    pos, cytoband = lrr_table["pos"].to_numpy(), lrr_table["cytoband"].astype(str).to_numpy()

    states = pd.Series(index=test.index, dtype=object)
    for sample, rows in test.groupby("sample"):
        limit = 4 * arrays.loc[sample, "lrr_sd"]
        clipped = np.clip(lrr_table[sample].to_numpy(np.float64), -limit, limit)
        for row in rows.itertuples():
            a, b = row.row_start, row.row_end
            piece = pd.DataFrame([dict(chr=row.chr, row_start=0, row_end=b - a)])
            call = calling.call_segments(piece, clipped[a:b] + row.lrr_shift, pos[a:b], cytoband[a:b],
                                         sigma=noise.loc[sample, "sigma_bin"], tau=noise.loc[sample, "tau_arm"],
                                         alpha=alpha, n_x=arrays.loc[sample, "n_X"]).iloc[0]
            if call["kind"] == "normal":
                states[row.Index] = "normal"
            else:
                direction = "loss" if call["change"] < 0 else "gain"
                states[row.Index] = direction if call["kind"] == "full" else f"mosaic {direction}"
    return states


def main():
    layout = paths.Layout()

    say("1. THE SAME SIMULATED SEGMENTS (left-out embryos)")
    test = pd.read_csv(tables.output_path("simulated_test.tsv"), sep="\t")
    test["state_rules"] = rules_on_simulated(test, layout)
    test[["sample", "row_start", "row_end", "state_rules"]].to_csv(
        tables.output_path("simulated_test_rules.tsv"), sep="\t", index=False)
    groups = {"all": test, f"under {LONG} SNPs": test[test.n_snps < LONG],
              f"{LONG} SNPs or more": test[test.n_snps >= LONG]}
    rows = [dict(method=name, segments=group, **scoring.state_metrics(part["true_state"], part[column]))
            for name, column in [("rules of step 06c", "state_rules"),
                                 ("ML: boosted trees", "state_boosted_trees"),
                                 ("ML: logistic regression", "state_logistic_regression")]
            for group, part in groups.items()]
    say(pd.DataFrame(rows).to_string(index=False, float_format="%.3f"))
    say()
    say("   Truth (rows) against the rules of step 06c (columns):")
    states = list(settings.STATES)
    say(pd.crosstab(test["true_state"], test["state_rules"]).reindex(index=states, columns=states,
                                                                    fill_value=0).to_string())
    say()

    say("2. THE REAL CALLS AGAINST THE PAPER'S RESULTS")
    reported = pd.read_csv(tables.output_path("reported_abnormal.tsv"), sep="\t")
    samples = list(pd.read_csv(layout.arrays, sep="\t", index_col=0).index)
    rule_calls = pd.read_csv(layout.cnv_calls("ASPCF"), sep="\t")
    rule_calls["direction"] = np.where(rule_calls["change"] > 0, "gain", "loss")
    ml_calls = pd.read_csv(tables.output_path("cnv_calls_ML.tsv"), sep="\t")
    both = None
    for name, calls in [("rules of step 06c", rule_calls), (f"ML: {settings.MODEL_FOR_CALLS}", ml_calls)]:
        metrics, _, _ = scoring.score(calls, reported, samples)
        metrics[name] = metrics["x"].astype(str) + " / " + metrics["n"].astype(str) \
            + metrics["value"].map(lambda v: f"  ({v:.0%})")
        both = metrics[["group", "name", name]] if both is None else both.assign(**{name: metrics[name]})
    say(both.to_string(index=False))

    tables.output_path("comparison_with_rules.txt").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
