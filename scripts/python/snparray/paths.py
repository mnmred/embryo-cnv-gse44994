"""Where every file lives.

All step scripts read inputs and write outputs through `Layout`, so the folder
structure is defined here and nowhere else.

    data/raw/cel/                        CEL files (created by scripts/setup_data.py)
    data/meta/GSE44994_series_matrix.txt sample sheet + the authors' genotype calls
    data/processed/                      everything this pipeline writes (default Layout)
"""
from dataclasses import dataclass
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]          # scripts/python/snparray/paths.py -> repo
CEL_DIR = REPO / "data" / "raw" / "cel"
SERIES_MATRIX = REPO / "data" / "meta" / "GSE44994_series_matrix.txt"
PROCESSED = REPO / "data" / "processed"


@dataclass(frozen=True)
class Layout:
    """File names inside one output folder (default: data/processed)."""

    root: Path = PROCESSED

    # step 01 - sample sheet and the authors' genotype calls
    @property
    def metadata(self) -> Path:
        return self.root / "metadata.tsv"

    @property
    def genotype_calls(self) -> Path:
        return self.root / "genotype_calls.parquet"

   
