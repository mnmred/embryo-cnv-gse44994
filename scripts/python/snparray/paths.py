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

    # step 04 - LRR 
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

    @property
    def baf(self) -> Path:
        return self.root / "baf.parquet"

    @property
    def snp_info(self) -> Path:
        """one row per SNP: genotype-cluster centres and population B-allele frequency"""
        return self.root / "snp_info.parquet"

    @property
    def baf_qc(self) -> Path:
        return self.root / "baf_qc.tsv"

    # step 06 
    @property
    def lrr_tsv(self) -> Path:
        """lrr.parquet as a gzipped text table, for R (step 06a)"""
        return self.root / "lrr.tsv.gz"

    @property
    def baf_tsv(self) -> Path:
        """baf.parquet as a gzipped text table, for R (step 06a)"""
        return self.root / "baf.tsv.gz"

    def segments_raw(self, method: str) -> Path:
        """pieces cut by an R method (step 06b: ASPCF), before copy numbers"""
        return self.root / f"segments_raw_{method}.tsv"

    def segments(self, method: str) -> Path:
        """every segment of every array with its tests and final copy number
        (method = ASPCF from step 06c, or HMM from the alternative step)"""
        return self.root / f"segments_{method}.tsv"

    def cnv_calls(self, method: str) -> Path:
        """only the reported (abnormal) segments of one method"""
        return self.root / f"cnv_calls_{method}.tsv"

    @property
    def call_noise(self) -> Path:
        """per array: noise of 50-SNP bins (sigma) and of whole arms (tau), step 06c"""
        return self.root / "call_noise.tsv"

