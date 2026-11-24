# Copy-number changes in human embryos from SNP arrays (GSE44994)

This repo finds copy-number changes in 47 embryo biopsies, starting from the raw CEL files
(Affymetrix Mapping 250K Nsp), and checks the calls against the authors' results and the parents' karyotypes.

## Tasks and code

| Task | Method (packages) | Code |
|---|---|---|
| LRR and BAF from CEL files | Quantile normalisation and median polish (`oligo`), then LRR and BAF written in Python | `scripts/R/02_*`, `scripts/python/01_*`–`05_*` |
| Segments and copy numbers | ASPCF (`copynumber`), then rules calibrated on the male X | `scripts/python/06a_*`, `scripts/R/06b_*`, `scripts/python/06c_*` |
| Quality check | Against the authors' results and the parents' karyotypes | `scripts/python/evaluation/` |
| Existing methods | CRMA v2 (`aroma.affymetrix`) and ASCAT | `scripts/R/ascat/` |
| Something extra | Gradient boosting (`scikit-learn`) trained on simulated segments | `scripts/python/ml/` |

## Environment

You need conda. From the repo folder, run:

```bash
conda install -n base -c conda-forge conda-lock
conda-lock install -n bioenv conda-lock.yml
conda env create -f ascat_env.yml
conda env create -f ml_env.yml
```

## Data

1. Download these two files from [GEO GSE44994](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44994):
   - `GSE44994_RAW.tar` (supplementary file, ~1.1 GB): the raw CEL files and the authors' CHP files
   - `GSE44994_series_matrix.txt.gz` (series matrix): sample metadata and karyotypes
2. Put both files in `data/raw/`.
3. Unpack them:

```bash
conda activate bioenv
python scripts/setup_data.py
```

## Run

Run everything from the repo root, in this order.

**1. Main pipeline** (`bioenv`; step 02 needs ~10 GB RAM)

```bash
python scripts/python/01_parse_metadata.py
Rscript scripts/R/02_normalize_cel.R
python scripts/python/03_summarize_snps.py
python scripts/python/04_compute_lrr.py
python scripts/python/05_compute_baf.py
python scripts/python/06a_export_for_r.py
Rscript scripts/R/06b_segment.R
python scripts/python/06c_call_copy_number.py
```

**2. CRMA v2 and ASCAT** (`ascat_env`)

```bash
Rscript scripts/R/ascat/00_install_sfit.R
Rscript scripts/R/ascat/01_crma_lrr_baf.R
Rscript scripts/R/ascat/02_run_ascat.R data/processed data/ascat/from_pipeline
Rscript scripts/R/ascat/02_run_ascat.R data/ascat/crma data/ascat/from_crma
```

**3. Gradient boosting** (`ml_env`)

```bash
python scripts/python/ml/00_paper_results_to_table.py
python scripts/python/ml/01_train_and_call.py
python scripts/python/ml/02_evaluate.py
python scripts/python/ml/03_compare_with_rules.py
```

**4. Evaluation** (`bioenv`). This runs steps 06b and 06c again, now on CRMA's LRR and BAF:

```bash
python scripts/python/evaluation/01_prepare_crma.py
Rscript scripts/R/06b_segment.R data/evaluation/crma_rules
python scripts/python/06c_call_copy_number.py --out-dir data/evaluation/crma_rules
python scripts/python/evaluation/02_evaluate.py
```

To include gradient boosting on CRMA, do this before running `02_evaluate.py`:

1. Point the file paths in `scripts/python/ml/inputs.py` to `data/evaluation/crma_rules`.
2. Set `OUTPUT_FOLDER` to `data/evaluation/crma_ml`.
3. Run `01_train_and_call.py` again.

**5. Figures** (`bioenv`). Written to `results/figures/evaluation/`:

```bash
python scripts/python/evaluation/03_figures.py
Rscript scripts/R/07_evaluation_figures.R
```

Data: Xiong et al. 2014, [Genomics Data 2:92–95](https://pmc.ncbi.nlm.nih.gov/articles/PMC4535754/).

AI Use: AI code assistance was used to accelerate the coding process. Any generated code has been carefully reviewed and adapted. The author takes full responsibility for the code.
