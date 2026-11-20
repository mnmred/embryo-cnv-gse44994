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
SERIES_MATRIX = REPO / "data" / "meta" / "GSE44994_series_matrix.txt"
PROCESSED = REPO / "data" / "processed"


@dataclass(frozen=True)
class Layout:
    """File names inside one output folder (default: data/processed)."""

    root: Path = PROCESSED

    @property
    def metadata(self) -> Path:
        return self.root / "metadata.tsv"

    @property
    def genotype_calls(self) -> Path:
        return self.root / "genotype_calls.parquet"
    
    @property
    def probe_dir(self) -> Path:
        return self.root / "probe_level"

    @property
    def probe_intensities(self) -> Path:
        return self.probe_dir / "pm_log2.f32"

    @property
    def probe_table(self) -> Path:
        return self.probe_dir / "pm_probes.tsv"


    @property
    def cel_order(self) -> Path:
        return self.probe_dir / "cel_order.txt"

    def signal(self, kind: str) -> Path:
        """kind = 'total' (all probes), 'A' or 'B' (allele-specific probes)."""
        return self.root / f"signal_{kind}.parquet"

    # step 04 - LRR -----------------------------------------------------------------
    @property
    def lrr(self) -> Path:
        return self.root / "lrr.parquet"

    @property
    def arrays(self) -> Path:
        """one row per array: sex call and LRR noise measures"""
        return self.root / "arrays.tsv"

    @property
    def calibration(self) -> Path:
        return self.root / "calibration.tsv"

    @property
    def snp_annotation(self) -> Path:
        return self.probe_dir / "snp_annotation.tsv"

   
