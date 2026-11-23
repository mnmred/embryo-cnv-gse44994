"""Step 2 - how good are the calls of step 1?

Run (repo root, ml_env active):   python scripts/python/ml/02_evaluate.py
Reads:   segments_ML.tsv, cnv_calls_ML.tsv, simulated_test.tsv (step 1) and, from inputs.py,
         the samples table and the reported results (if any)
Writes:  evaluation.txt        everything printed here
         evaluation.tsv        the scores against the reported results, one row per score
         cnv_calls_ML.tsv      gets two columns: overlap, confirmed (by a reported result)

Three tests, from most artificial to most real:
    1. simulated segments of left-out embryos   the truth is known exactly, but simulated
    2. the X chromosome                         real: one copy in males, two in females
    3. the reported results                     real, but only as good as the report
"""
import pandas as pd

import inputs
import scoring
import settings
import tables
from decisions import STATE_NAMES, add_copies, decide, join_neighbours

MODELS = ["boosted_trees", "logistic_regression"]
LONG = 1000                     # SNPs: segments are also scored separately above / below this
OTHER_MIN_PROBABILITIES = [0.5, 0.9, 0.95, 0.99]     # what if settings.MIN_PROBABILITY were different?
lines = []


def say(text=""):
    print(text)
    lines.append(str(text))


def simulated_test(test: pd.DataFrame):
    say("1. SIMULATED SEGMENTS OF LEFT-OUT EMBRYOS")
    say(f"   {len(test)} segments, each called by a model that never saw its embryo.")
    groups = {"all": test, f"under {LONG} SNPs": test[test.n_snps < LONG],
              f"{LONG} SNPs or more": test[test.n_snps >= LONG]}
    rows = []
    for model in MODELS:
        for group, part in groups.items():
            rows.append(dict(model=model, segments=group,
                             **scoring.state_metrics(part["true_state"], part[f"state_{model}"])))
    say(pd.DataFrame(rows).to_string(index=False, float_format="%.3f"))
    say()
    say(f"   Truth (rows) against the final state of {settings.MODEL_FOR_CALLS} (columns):")
    states = list(settings.STATES)
    table = pd.crosstab(test["true_state"], test[f"state_{settings.MODEL_FOR_CALLS}"])
    say(table.reindex(index=states, columns=states, fill_value=0).to_string())
    say()


def x_chromosome_test(segments: pd.DataFrame):
    say("2. THE X CHROMOSOME (real data)")
    say("   Share of the X SNPs per most probable state. Males have one X: the right answer")
    say("   is 'loss'. Females have two: the right answer is 'normal'.")
    samples = tables.read_table(inputs.SAMPLES_TABLE)
    n_x = dict(zip(samples[inputs.SAMPLE_NAME_COLUMN].astype(str), samples[inputs.SAMPLE_X_COPIES_COLUMN]))
    on_x = segments[(segments["chr"] == inputs.X_CHROMOSOME) & (segments["state"] != "too short")].copy()
    on_x["sex"] = on_x["sample"].map(lambda s: "males" if n_x[s] == 1 else "females")
    on_x = on_x[(on_x["sex"] == "females") | (on_x["normal_state"] == "loss")]     # not the ends also on Y
    rows = []
    for model in MODELS:
        for sex, part in on_x.groupby("sex"):
            share = part.groupby(f"most_probable_{model}")["n_snps"].sum() / part["n_snps"].sum()
            rows.append(dict(model=model, sex=sex, embryos=part["sample"].nunique(),
                             **{s: share.get(s, 0.0) for s in settings.STATES}))
    say(pd.DataFrame(rows).to_string(index=False, float_format="%.3f"))
    say()


def reported_results_test(calls: pd.DataFrame, segments: pd.DataFrame):
    say("3. THE REPORTED RESULTS (real data)")
    if inputs.REPORTED_RESULTS_TABLE is None:
        say("   no reported results given in inputs.py")
        return
    reported = tables.read_table(inputs.REPORTED_RESULTS_TABLE)
    samples = list(tables.read_table(inputs.SAMPLES_TABLE)[inputs.SAMPLE_NAME_COLUMN].astype(str))
    metrics, reported, calls = scoring.score(calls.drop(columns=["overlap", "confirmed"], errors="ignore"),
                                             reported, samples)
    say(metrics.to_string(index=False, float_format="%.3f"))
    say()
    say("   Reported but not found:")
    missed = reported.loc[~reported["found"], ["sample", "chr", "start", "end", "direction"]]
    say(missed.to_string(index=False) if len(missed) else "   none")
    say()
    say("   Called but not reported:")
    extra = calls.loc[~calls["confirmed"], ["sample", "chr", "start", "end", "size_mb", "state",
                                             "probability_abnormal"]]
    say(extra.to_string(index=False) if len(extra) else "   none")
    say()
    say("   The same scores if MIN_PROBABILITY (settings.py) had another value. Do not pick the")
    say("   value that looks best here and then quote these scores: that would be tuning on the answer.")
    called = segments[segments["state"] != "too short"].copy()
    p = called[[f"p_{s.replace(' ', '_')}" for s in STATE_NAMES]].set_axis(STATE_NAMES, axis=1)
    rows = []
    for min_probability in OTHER_MIN_PROBABILITIES:
        called["state"] = decide(p, called["normal_state"], min_probability)["state"]
        other, _, _ = scoring.score(join_neighbours(add_copies(called)), reported, samples)
        value = other.set_index(["group", "name"])["value"]
        rows.append(dict(min_probability=min_probability, calls=int(other.loc[3, "n"]),
                         recall=value["per region (any overlap)", "recall"],
                         precision=value["per region (any overlap)", "precision"],
                         embryo_specificity=value["per embryo", "specificity"]))
    say(pd.DataFrame(rows).to_string(index=False, float_format="%.3f"))
    metrics.to_csv(tables.output_path("evaluation.tsv"), sep="\t", index=False, float_format="%.4g")
    calls.to_csv(tables.output_path("cnv_calls_ML.tsv"), sep="\t", index=False, float_format="%.4g")


def main():
    simulated_test(pd.read_csv(tables.output_path("simulated_test.tsv"), sep="\t"))
    segments = pd.read_csv(tables.output_path("segments_ML.tsv"), sep="\t")
    x_chromosome_test(segments)
    reported_results_test(pd.read_csv(tables.output_path("cnv_calls_ML.tsv"), sep="\t"), segments)
    tables.output_path("evaluation.txt").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
