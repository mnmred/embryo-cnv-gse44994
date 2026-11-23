"""THE FILE TO EDIT to change what 02_evaluate.py compares.

Every path is relative to the repository root.
"""

# ------------------------------------------------------------------ 1. the setups compared
# One entry per setup. `name` is its label in the report; `input` (the LRR/BAF used) and `caller` describe it.
# `reader` says how to read its calls:
#   rules   a cnv_calls_ASPCF.tsv of step 06c        (direction = sign of `change`)
#   ascat   a segments_ASCAT_<setting>.tsv           (every segment with its copy number)
#   ml      a cnv_calls_ML.tsv of the ML caller      (column `direction`: gain / loss)
# `run` is what to run when the file does not exist yet (shown in the report).
ASCAT_SETTING = "calibrated"     # "default", "fixed" or "calibrated" (see scripts/R/ascat/02_run_ascat.R)

SETUPS = [
    dict(name="Pipeline end to end", input="pipeline", caller="ASPCF + rules (06c)", reader="rules",
         calls="data/processed/cnv_calls_ASPCF.tsv",
         run="the main pipeline, steps 01 to 06c"),
    dict(name="CRMA + pipeline end", input="CRMA v2", caller="ASPCF + rules (06c)", reader="rules",
         calls="data/evaluation/crma_rules/cnv_calls_ASPCF.tsv",
         run="steps 1 to 3 of scripts/python/evaluation/README.md"),
    dict(name="CRMA + ASCAT", input="CRMA v2", caller="ASCAT", reader="ascat",
         calls=f"data/ascat/from_crma/segments_ASCAT_{ASCAT_SETTING}.tsv",
         run="scripts/R/ascat/02_run_ascat.R data/ascat/crma data/ascat/from_crma"),
    dict(name="Pipeline + ASCAT", input="pipeline", caller="ASCAT", reader="ascat",
         calls=f"data/ascat/from_pipeline/segments_ASCAT_{ASCAT_SETTING}.tsv",
         run="scripts/R/ascat/02_run_ascat.R data/processed data/ascat/from_pipeline"),
    dict(name="CRMA + GBoost", input="CRMA v2", caller="ASPCF + GBoost", reader="ml",
         calls="data/evaluation/crma_ml/cnv_calls_ML.tsv",
         run="steps 1 to 4 of scripts/python/evaluation/README.md (step 4: edit ml/inputs.py)"),
    dict(name="Pipeline + GBoost", input="pipeline", caller="ASPCF + GBoost", reader="ml",
         calls="data/ml/cnv_calls_ML.tsv",
         run="scripts/python/ml/01_train_and_call.py"),
]

# ------------------------------------------------------------------ 2. the references
METADATA = "data/processed/metadata.tsv"            # reported_result (authors), donor_karyotype (parent), couple
AUTHORS = "data/ml/reported_abnormal.tsv"            # the authors' results as regions (ml/00_paper_results_to_table.py)
SNPS = "data/processed/lrr.parquet"                  # only chr, pos, cytoband: to place bands and chromosome ends
ARRAYS = "data/processed/arrays.tsv"                 # n_X: one or two X per embryo

# ------------------------------------------------------------------ 3. choices
MIN_CALL_MB = 2.36             # calls shorter than this are dropped for every setup (the authors' smallest change)
MIN_RECIPROCAL_OVERLAP = 0.5   # the strict event match
END_TOLERANCE_MB = 1.0         # an edge this close to the first / last SNP of a chromosome is a chromosome end
PERMUTATIONS = 10_000          # chance level of the parent test: karyotypes shuffled among couples
SEED = 1

# ------------------------------------------------------------------ 4. LRR and BAF for figures
# The LRR and BAF of every `input` of SETUPS, written as profile_bins.tsv.gz and profile_snps.tsv.gz.
# BAF is clipped to 0-1, as step 06a does. Set WRITE_PROFILES = False to skip them (they take ~1 min).
WRITE_PROFILES = True
PROFILES = {
    "pipeline": dict(lrr="data/processed/lrr.parquet", baf="data/processed/baf.parquet"),
    "CRMA v2": dict(lrr="data/evaluation/crma_rules/lrr.parquet", baf="data/evaluation/crma_rules/baf.tsv.gz"),
}
PROFILE_BIN_SNPS = 50          # SNPs per bin (as the noise bins of step 06c)
PROFILE_EVERY_SNP = 10         # profile_snps keeps one SNP in this many
HET_BAF = (0.2, 0.8)           # a SNP whose BAF is in this range looks heterozygous

# ------------------------------------------------------------------ 5. outputs
REPORT = "results/evaluation.md"
TABLE_FOLDER = "data/evaluation"     # the tables for figures (see the README)
