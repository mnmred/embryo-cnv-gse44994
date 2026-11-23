"""Choices of the ML caller that are not file names (those are in inputs.py)."""

# ------------------------------------------------------------------ copy-number states
# The model gives every segment one of these states. For each state: the true copy
# numbers that are simulated for it, and the copy number written to the output.
#   normal        2 copies in every cell
#   loss / gain   one copy lost / gained in every cell (gain: sometimes two gained)
#   mosaic ...    the change is in 40-80 % of the cells, so 0.4-0.8 copies on average
#                 (the same definition as step 06c of the main pipeline)
STATES = {
    #  state          lowest, highest simulated copies    reported copies
    "loss":          dict(simulate=(1.0, 1.0),           copies=1.0),
    "mosaic loss":   dict(simulate=(1.2, 1.6),           copies=1.5),
    "normal":        dict(simulate=(2.0, 2.0),           copies=2.0),
    "mosaic gain":   dict(simulate=(2.4, 2.8),           copies=2.5),
    "gain":          dict(simulate=(3.0, 3.0),           copies=3.0),
}
SHARE_NORMAL = 0.5            # half of the simulated segments stay normal, the rest is split evenly
SHARE_GAIN_WITH_4_COPIES = 0.15   # part of the simulated gains have 4 copies instead of 3

# ------------------------------------------------------------------ simulation
SIMULATED_PER_EMBRYO_TRAIN = 300   # simulated segments per embryo to train on
SIMULATED_PER_EMBRYO_TEST = 100    # simulated segments per left-out embryo to test on
MIN_SEGMENT_SNPS = 100             # shorter segments are not simulated and not called
SHARE_WHOLE_CHROMOSOME = 0.3       # part of the simulated segments that cover a whole chromosome

# Used only if the data has no male: LRR drop of one lost copy (k of the main pipeline).
ONE_COPY_SHIFT_IF_NO_MALES = 0.28

# ------------------------------------------------------------------ noise
BIN_SNPS = 50                 # LRR is averaged per 50 neighbouring SNPs to measure noise
CLIP_AT_SD = 4.0              # single SNPs further than this from 0 are pulled in

# ------------------------------------------------------------------ model and testing
MODEL_FOR_CALLS = "boosted_trees"  # "boosted_trees" or "logistic_regression"
FOLDS = 5                     # embryos are split in 5 groups; each group is called by a
                              # model trained on the other 4 (it never saw those embryos)
MIN_PROBABILITY = 0.9         # an abnormal state is reported only if the model is this sure
SEED = 1
