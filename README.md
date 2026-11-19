## How to set up the environment

Requires conda. From the repo folder:

```bash
conda install -n base -c conda-forge conda-lock
conda-lock install -n bioenv conda-lock.yml
conda activate bioenv
```

## Data setup

1. **Download the data** from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44994
   - `GSE44994_RAW.tar` (Supplementary file, select all files): raw CEL + CHP files, ~1.1 GB
   - `GSE44994_series_matrix.txt.gz` (Series Matrix File): sample metadata and karyotypes

2. **Put both files in `data/raw/`**

3. **Unpack the data in the environment:**
   ```bash
   conda activate bioenv
   python scripts/setup_data.py
   ```

This creates:

| Path | Contents |
|---|---|
| `data/raw/cel/` | Raw intensity files, input for the analysis |
| `data/raw/chp/` | Authors' Affymetrix output, used for comparison |
| `data/meta/samples.csv` | One row per sample: GSM ID, sample name, file paths |
| `data/meta/GSE44994_series_matrix.txt` | Sample metadata (karyotypes, descriptions) |
