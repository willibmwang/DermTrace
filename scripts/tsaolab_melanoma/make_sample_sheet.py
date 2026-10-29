"""Sample sheet for the Tsao lab Greek melanoma FFPE RNA-seq (core project 20088-11).

The core delivered one folder of paired FASTQs (one lane per sample) plus md5.txt and no sample key, so the
sheet is built from md5.txt: every sample the core listed gets a row, whether or not its FASTQs are here yet.
Clinical annotation (patient, site, stage, BRAF/NRAS status from the pathology report) is not in the delivery and
has to be added later from the lab's own key.

Usage: python -I make_sample_sheet.py <raw_data_dir> <out.csv>
"""
import re
import sys
from pathlib import Path

import pandas as pd

PROJECT = "20088-11"
PAT = re.compile(r"^(?P<core_id>20088XR-11-(?P<n>\d+))_S(?P<s>\d+)_(?P<lane>L\d{3})_R(?P<read>[12])_001\.fastq\.gz$")


def main(raw, out):
    md5 = {}
    for line in (raw / "md5.txt").read_text().splitlines():
        if line.strip():
            h, f = line.split()
            md5[f] = h
    rows = {}
    for f, h in sorted(md5.items()):
        m = PAT.match(f)
        if not m:
            sys.exit(f"unexpected file name in md5.txt: {f}")
        n = int(m["n"])
        r = rows.setdefault(n, dict(
            sample_id=f"GM{n:02d}", core_id=m["core_id"], project=PROJECT, species="human",
            tissue="melanoma", specimen="FFPE", disease_label="melanoma",
            lane=m["lane"], index_s=f"S{m['s']}", r1="", r2="", md5_r1="", md5_r2="",
        ))
        assert r["lane"] == m["lane"] and r["index_s"] == f"S{m['s']}", f
        r["r" + m["read"]] = f
        r["md5_r" + m["read"]] = h
    df = pd.DataFrame([rows[k] for k in sorted(rows)])
    assert (df.r1 != "").all() and (df.r2 != "").all(), "sample without both reads in md5.txt"
    df.to_csv(out, index=False)
    present = sum((raw / f).exists() for f in df.r1) + sum((raw / f).exists() for f in df.r2)
    print(f"wrote {out}: {len(df)} samples (core numbers {df.core_id.str.split('-').str[-1].tolist()}), "
          f"{present}/{2 * len(df)} FASTQs present")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
