"""THE FILE TO EDIT when you want to run the ML caller on other data.

Every path is relative to the repository root. Tables may be .parquet, .tsv,
.csv or .txt (also gzipped, e.g. .tsv.gz); the format is read from the file name.
The steps 01 and 02 take all their file names from here.
"""

# ------------------------------------------------------------------ 1. input tables

# LRR: one row per SNP, one column per sample. Rows sorted by chromosome and position.
LRR_TABLE = "data/processed/lrr.parquet"
LRR_CHROMOSOME_COLUMN = "chr"        # numbers; X uses the number in X_CHROMOSOME below
LRR_POSITION_COLUMN = "pos"          # base pairs

# Segments: one row per piece that a segmentation method cut. The ML model gives each a state.
SEGMENTS_TABLE = "data/processed/segments_raw_ASPCF.tsv"
SEGMENT_SAMPLE_COLUMN = "sample"
SEGMENT_CHROMOSOME_COLUMN = "chr"
SEGMENT_START_COLUMN = "start"       # base pairs of the first SNP in the piece
SEGMENT_END_COLUMN = "end"           # base pairs of the last SNP in the piece

# Samples: one row per sample. Its names must be the sample columns of LRR_TABLE.
SAMPLES_TABLE = "data/processed/arrays.tsv"
SAMPLE_NAME_COLUMN = "sample"
SAMPLE_X_COPIES_COLUMN = "n_X"       # 1 = one X (male), 2 = two X (female)

# Chromosomes that may NOT be used as "normal" examples for training.
# Any number of tables; each needs a sample column and a chromosome column
# (named as in the segments table). Use [] if you have none.
KNOWN_ABNORMAL_TABLES = [
    "data/processed/cnv_calls_ASPCF.tsv",     # calls of the rule-based step 06c
    "data/ml/reported_abnormal.tsv",          # the paper's results (made by 00_paper_results_to_table.py)
]

# Optional: results reported by someone else, to compare the ML calls with.
# Columns: sample, chr, direction ("gain" or "loss"). Use None if you have none.
REPORTED_RESULTS_TABLE = "data/ml/reported_abnormal.tsv"

# ------------------------------------------------------------------ 2. about the genome

X_CHROMOSOME = 23
# Ends of X that are also on Y (pseudo-autosomal regions), so males have two copies there.
# These are hg19 positions; change them if your positions use another genome build.
X_PAR1_END = 2_699_520
X_PAR2_START = 154_931_044

# ------------------------------------------------------------------ 3. output folder

OUTPUT_FOLDER = "data/ml"
