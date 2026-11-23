"""Step 1 - train the model on simulated segments and give every real segment a state.

Run (repo root, ml_env active):   python scripts/python/ml/01_train_and_call.py     (~20 s)
Reads:   the tables named in inputs.py
Writes (to the output folder of inputs.py):
    segments_ML.tsv        every real segment: features, probability of each state, final state
    cnv_calls_ML.tsv       only the abnormal ones, neighbours with the same state joined
    simulated_test.tsv     simulated segments of left-out embryos, with the truth and the
                           answer of each model (read by 02_evaluate.py)

How an embryo gets its calls (no embryo is ever called by a model that has seen it):
    the embryos are split in FOLDS groups; for each group
      1. simulate training segments from the OTHER embryos           (simulation.py)
      2. fit the models on them                                      (make_models)
      3. simulate test segments from the embryos of this group
      4. let the models answer on those and on the group's real segments
"""
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

import inputs
import settings
import simulation
import tables
from decisions import STATE_NAMES, add_copies, decide, join_neighbours
from features import FEATURES, sample_noise, segment_features


def make_models() -> dict:
    return {
        "boosted_trees": HistGradientBoostingClassifier(random_state=settings.SEED),
        "logistic_regression": make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000)),
    }


def real_segment_table(data, noise) -> pd.DataFrame:
    """The real segments with their features and the state that is normal for them:
    "normal" (2 copies), or "loss" (1 copy) on the X of a male, outside the regions also on Y."""
    on_x_only = tables.x_only(data)
    rows = []
    for seg in data.segments.itertuples():
        j = data.samples.index(seg.sample)
        a, b = seg.row_start, seg.row_end
        male_x = data.n_x[j] == 1 and b > a and on_x_only[a:b].mean() >= 0.5
        row = dict(sample=seg.sample, chr=seg.chr, start=seg.start, end=seg.end,
                   size_mb=(seg.end - seg.start) / 1e6, row_start=a, row_end=b, n_snps=b - a,
                   normal_state="loss" if male_x else "normal")
        if b - a >= settings.MIN_SEGMENT_SNPS:
            row |= segment_features(data.lrr[a:b, j], noise.loc[seg.sample].to_dict())
        rows.append(row)
    return pd.DataFrame(rows)


def probabilities(model, table: pd.DataFrame) -> pd.DataFrame:
    """Probability of every state for every row, columns in the order of settings.STATES."""
    p = pd.DataFrame(model.predict_proba(table[FEATURES]), columns=model.classes_, index=table.index)
    return p.reindex(columns=STATE_NAMES, fill_value=0.0)


def main():
    rng = np.random.default_rng(settings.SEED)
    data = tables.load()
    noise = sample_noise(data)                      # also clips extreme SNPs
    shift_of_male = simulation.one_copy_shifts(data)
    all_shifts = list(shift_of_male.values()) or [settings.ONE_COPY_SHIFT_IF_NO_MALES]
    typical_shift = float(np.median(all_shifts))
    print(f"{len(data.samples)} samples, {len(data.segments)} segments, {len(shift_of_male)} males; "
          f"one lost copy lowers LRR by {typical_shift:.3f} (median of the males)")

    real = real_segment_table(data, noise)
    callable_ = real["n_snps"] >= settings.MIN_SEGMENT_SNPS
    fold_of_sample = dict(zip(rng.permutation(data.samples), np.arange(len(data.samples)) % settings.FOLDS))

    real_parts, test_parts = [], []
    for fold in range(settings.FOLDS):
        test_samples = [s for s in data.samples if fold_of_sample[s] == fold]
        train_samples = [s for s in data.samples if fold_of_sample[s] != fold]
        # the one-copy shift is learnt from the training males only
        train_shifts = [shift_of_male[s] for s in train_samples if s in shift_of_male] \
            or [settings.ONE_COPY_SHIFT_IF_NO_MALES]
        train = simulation.simulate(data, noise, train_samples, train_shifts,
                                    settings.SIMULATED_PER_EMBRYO_TRAIN, rng)
        test = simulation.simulate(data, noise, test_samples, all_shifts,
                                   settings.SIMULATED_PER_EMBRYO_TEST, rng)
        real_here = real[callable_ & real["sample"].isin(test_samples)].copy()
        test["fold"], real_here["fold"] = fold, fold

        for name, model in make_models().items():
            model.fit(train[FEATURES], train["true_state"])
            for table in (test, real_here):
                p = probabilities(model, table)
                table[f"most_probable_{name}"] = p.idxmax(axis=1)
                decided = decide(p, table.get("normal_state", pd.Series("normal", index=table.index)))
                table[f"state_{name}"] = decided["state"]
                table[f"probability_abnormal_{name}"] = decided["probability_abnormal"]
                if name == settings.MODEL_FOR_CALLS and table is real_here:
                    for state in STATE_NAMES:
                        real_here[f"p_{state.replace(' ', '_')}"] = p[state]
        real_parts.append(real_here)
        test_parts.append(test)
        print(f"group {fold + 1} of {settings.FOLDS}: trained on {len(train)} simulated segments "
              f"of {len(train_samples)} embryos, called {len(test_samples)} embryos")

    # ---- the real segments: final state, copy number, abnormal or not
    called = pd.concat(real_parts)
    called["state"] = called[f"state_{settings.MODEL_FOR_CALLS}"]
    called["probability_abnormal"] = called[f"probability_abnormal_{settings.MODEL_FOR_CALLS}"]
    called["copies_from_lrr"] = simulation.copies_from_lrr(called["lrr_median"], typical_shift)
    called = add_copies(called)

    too_short = real[~callable_].assign(state="too short", abnormal=False)
    segments = pd.concat([called, too_short]).sort_values(["sample", "row_start"])
    first = ["sample", "chr", "start", "end", "size_mb", "n_snps", "state", "copies", "normal_copies",
             "abnormal", "probability_abnormal", "copies_from_lrr"]
    segments = segments[first + [c for c in segments.columns if c not in first]]
    calls = join_neighbours(called)

    segments.to_csv(tables.output_path("segments_ML.tsv"), sep="\t", index=False, float_format="%.4g")
    calls.to_csv(tables.output_path("cnv_calls_ML.tsv"), sep="\t", index=False, float_format="%.4g")
    pd.concat(test_parts).to_csv(tables.output_path("simulated_test.tsv"), sep="\t", index=False,
                                 float_format="%.5g")
    print(f"{len(segments)} segments ({len(too_short)} too short to call) -> {inputs.OUTPUT_FOLDER}/segments_ML.tsv")
    print(f"{len(calls)} abnormal regions in {calls['sample'].nunique()} of {len(data.samples)} samples "
          f"-> {inputs.OUTPUT_FOLDER}/cnv_calls_ML.tsv")
    print(calls["state"].value_counts().to_string())


if __name__ == "__main__":
    main()
