"""Build sample sheets (CSV) for the Quintara RNA-seq orders from Chelvi's xlsx submission forms.

Usage: python -I make_sample_sheets.py <datasets_dir> <out_dir>
"""
import re
import sys
from pathlib import Path

import pandas as pd

# Every cell line was run through the same 4 oxygen arms, in this order, in triplicate.
OXYGEN_ARMS = ["20%", "20%->1%", "1%", "1%->20%"]


def oxygen_from_name(name):
    pcts = re.findall(r"(\d+)\s*%", name)
    return "->".join(f"{p}%" for p in pcts)


def sheet_0013831748(datasets):
    f = datasets / "0013831748-2024 RNAseq-Chelvi - Quintara -Differential Analysis.xlsx"
    key = pd.read_excel(f, sheet_name="Sample-Key")
    key.columns = ["sample_name", "qb_id", "qb_sample", "species", "group"]
    rows = []
    for r in key.itertuples():
        n = int(r.qb_id.split("-")[1])
        submitted_as = "HL-1" if n <= 12 else "NIH3T3" if n <= 24 else "BJ"
        # RNA-seq shows samples 13-24 are not mouse NIH3T3: ~93% of reads map to human, they are female,
        # fibroblast genes are off, and every sample carries adenovirus-5 E1A/E1B transcripts (HEK293 signature).
        cell_line = "HEK293-like" if submitted_as == "NIH3T3" else submitted_as
        rows.append(dict(
            sample_id=f"Q{n}",
            qb_id=r.qb_id,
            sample_name=r.sample_name.strip(),
            project="0013831748",
            species="human" if cell_line == "HEK293-like" else r.species.strip().lower(),
            submitted_as=submitted_as,
            cell_line=cell_line,
            cell_type={"HL-1": "cardiomyocyte", "HEK293-like": "HEK293_like", "BJ": "fibroblast"}[cell_line],
            disease_label="normal",
            group=r.group.replace("Group_", ""),
            oxygen=OXYGEN_ARMS[(n - 1) // 3 % 4],
            collection_date="2024-07-30",
        ))
    return pd.DataFrame(rows)


def sheet_1030390(datasets):
    f = datasets / "1030390-2025 ATCC & AC16_RNAseq-6172025-Chelvi with Groups.xlsx"
    raw = pd.read_excel(f, header=None, skiprows=1)
    raw = raw[raw[0].astype(str).str.startswith("Group")]
    rows = []
    for i, r in enumerate(raw.itertuples(index=False)):
        n = i + 1  # column 2 has a typo (17 entered as 27); row order is authoritative
        name = str(r[3]).strip()
        cell_line = "ATCC_keloid" if n <= 12 else "AC16"
        rows.append(dict(
            sample_id=f"A{n}",
            qb_id=f"1030390-{n}",  # assumed: same Quintara folder naming as 0013831748
            sample_name=name,
            project="1030390",
            species="human",
            submitted_as=cell_line,
            cell_line=cell_line,
            cell_type="fibroblast" if n <= 12 else "cardiomyocyte",
            disease_label="keloid" if n <= 12 else "normal",
            group=str(r[0]).replace("Group ", ""),
            oxygen=OXYGEN_ARMS[(n - 1) // 3 % 4],
            collection_date=pd.Timestamp(r[15]).date().isoformat(),
        ))
    df = pd.DataFrame(rows)
    # sanity check: the oxygen arm parsed from the name must match the position-based arm
    assert (df.sample_name.map(oxygen_from_name) == df.oxygen).all(), df
    return df


if __name__ == "__main__":
    datasets, out = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    for fn in (sheet_0013831748, sheet_1030390):
        df = fn(datasets)
        if fn is sheet_0013831748:
            assert (df.sample_name.map(oxygen_from_name) == df.oxygen).all(), df
        path = out / f"samples_{df.project[0]}.csv"
        df.to_csv(path, index=False)
        print(f"wrote {path} ({len(df)} samples)")
