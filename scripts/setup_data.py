#!/usr/bin/env python3
"""
Unpack the GSE44994 raw archive (GSE44994_RAW.tar) into a clean project layout.
Safe to re-run: files that already exist are skipped.
"""
import argparse
import csv
import gzip
import re
import shutil
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
CEL_DIR = RAW / "cel"
CHP_DIR = RAW / "chp"
META = ROOT / "data" / "meta"
TMP = RAW / "_extract_tmp"


def find_archive(given):
    if given:
        p = Path(given).expanduser().resolve()
        if not p.exists():
            sys.exit(f"Archive not found: {p}")
        return p
    candidates = sorted(RAW.glob("*RAW*.tar"))
    if not candidates:
        sys.exit(f"No .tar found in {RAW}. Put GSE44994_RAW.tar there or pass its path.")
    return candidates[0]


def safe_extract_tar(path, dest):
    with tarfile.open(path) as tf:
        try:
            tf.extractall(dest, filter="data")  # Python >= 3.11.4: blocks path traversal
        except TypeError:  # older Python without the filter argument
            for m in tf.getmembers():
                target = (dest / m.name).resolve()
                if not str(target).startswith(str(dest.resolve())) or not (m.isfile() or m.isdir()):
                    raise RuntimeError(f"Refusing unsafe archive member: {m.name}")
            tf.extractall(dest)


def gunzip_to(src, dst):
    with gzip.open(src, "rb") as fin, open(dst, "wb") as fout:
        shutil.copyfileobj(fin, fout, length=1 << 20)


def place(src, kind, keep_gz):
    """Move one CEL/CHP file (gzipped or not) into its final folder."""
    out_dir = CEL_DIR if kind == "CEL" else CHP_DIR
    is_gz = src.name.lower().endswith(".gz")
    if is_gz and kind == "CEL" and keep_gz:
        dst = out_dir / src.name
        if not dst.exists():
            shutil.move(str(src), dst)
        return dst
    dst = out_dir / (src.name[:-3] if is_gz else src.name)
    if dst.exists():
        return dst
    if is_gz:
        gunzip_to(src, dst)
    else:
        shutil.move(str(src), dst)
    return dst


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("archive", nargs="?", help="path to GSE44994_RAW.tar (default: search data/raw/)")
    ap.add_argument("--keep-gz", action="store_true", help="keep CEL files gzipped (oligo reads .CEL.gz)")
    args = ap.parse_args()

    for d in (CEL_DIR, CHP_DIR, META):
        d.mkdir(parents=True, exist_ok=True)

    archive = find_archive(args.archive)
    print(f"Unpacking {archive.name} ...")
    if TMP.exists():
        shutil.rmtree(TMP)
    TMP.mkdir(parents=True)
    safe_extract_tar(archive, TMP)

    samples = {}
    n = {"CEL": 0, "CHP": 0}
    files = sorted(p for p in TMP.rglob("*") if p.is_file())
    for i, p in enumerate(files, 1):
        m = re.search(r"\.(CEL|CHP)(\.gz)?$", p.name, re.I)
        if not m:
            print(f"\n  skipping unrecognised file: {p.name}")
            continue
        kind = m.group(1).upper()
        dst = place(p, kind, args.keep_gz)
        n[kind] += 1
        gsm = re.match(r"(GSM\d+)", dst.name)
        gsm = gsm.group(1) if gsm else dst.name.split(".")[0]
        sname = re.sub(r"^GSM\d+_", "", dst.name).split(".")[0]
        rec = samples.setdefault(gsm, {"gsm": gsm, "sample": sname, "cel": "", "chp": ""})
        rec[kind.lower()] = str(dst.relative_to(ROOT))
        print(f"\r  {i}/{len(files)} files processed", end="", flush=True)
    print()
    shutil.rmtree(TMP)

    # Series matrix (downloaded separately from GEO) -> data/meta/
    for sm in RAW.glob("*series_matrix*"):
        if "Zone.Identifier" in sm.name:
            continue
        dst = META / (sm.name[:-3] if sm.name.endswith(".gz") else sm.name)
        if not dst.exists():
            if sm.name.endswith(".gz"):
                gunzip_to(sm, dst)
            else:
                shutil.copy2(sm, dst)
            print(f"Series matrix -> {dst.relative_to(ROOT)}")

    manifest = META / "samples.csv"
    with open(manifest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["gsm", "sample", "cel", "chp"])
        w.writeheader()
        for k in sorted(samples):
            w.writerow(samples[k])

    missing_cel = [s["gsm"] for s in samples.values() if not s["cel"]]
    print(f"\nDone: {n['CEL']} CEL and {n['CHP']} CHP files for {len(samples)} samples.")
    print(f"Manifest: {manifest.relative_to(ROOT)}")
    if missing_cel:
        print(f"WARNING: no CEL file for {len(missing_cel)} samples: {', '.join(missing_cel[:10])}")
    if not any(META.glob("*series_matrix*")):
        print("Reminder: also download GSE44994_series_matrix.txt.gz from GEO into data/raw/ and re-run.")


if __name__ == "__main__":
    main()