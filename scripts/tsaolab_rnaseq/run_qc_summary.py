"""One-row-per-sample sequencing/quantification QC from the fastp and salmon outputs of run_quant.sh.

Usage: python -I run_qc_summary.py <out.tsv> <sample_sheet.csv>=<quant_dir>=<raw_data_dir> [...]
"""
import gzip
import json
import sys
from pathlib import Path

import pandas as pd

INSTRUMENTS = {"AV": "Element_AVITI", "NB": "Illumina_NextSeq", "NS": "Illumina_NextSeq"}


def sequencing_run(raw, qb):
    """Instrument and run from the first read header of each R1 FASTQ (all lanes must agree)."""
    runs = set()
    for f in sorted((raw / qb).glob(f"{qb}_S*_R1_001.fastq.gz")):
        with gzip.open(f, "rt") as fh:
            inst, run = fh.readline()[1:].split(":")[:2]
        runs.add((INSTRUMENTS.get(inst[:2], inst), f"{inst}:{run}"))
    assert len(runs) == 1, (qb, runs)
    inst, run = runs.pop()
    return dict(instrument=inst, sequencing_run=run)


def fastp_totals(qdir, sid):
    before = after = adapter = q30_bases = bases = 0
    dup = []
    for j in sorted((qdir / "fastp").glob(f"{sid}_L00*.json")):
        d = json.loads(j.read_text())
        s = d["summary"]
        before += s["before_filtering"]["total_reads"]
        after += s["after_filtering"]["total_reads"]
        q30_bases += s["after_filtering"]["q30_bases"]
        bases += s["after_filtering"]["total_bases"]
        adapter += d.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0)
        dup.append(d.get("duplication", {}).get("rate", float("nan")))
    return dict(
        fastp_lanes=len(dup),
        reads_raw=before,
        reads_after_trim=after,
        pct_reads_kept=round(after / before * 100, 2),
        pct_adapter_trimmed=round(adapter / before * 100, 2),
        pct_q30_after_trim=round(q30_bases / bases * 100, 2),
        fastp_duplication_rate=round(sum(dup) / len(dup), 4),
    )


def main(out, runs):
    rows = []
    for sheet, qdir, raw in runs:
        for r in pd.read_csv(sheet, dtype=str).itertuples():
            sq = qdir / "quant" / r.sample_id
            meta = json.loads((sq / "aux_info" / "meta_info.json").read_text())
            lib = json.loads((sq / "lib_format_counts.json").read_text())
            rows.append(dict(
                sample_id=r.sample_id, qb_id=r.qb_id, species=r.species, cell_line=r.cell_line,
                **sequencing_run(raw, r.qb_id),
                **fastp_totals(qdir, r.sample_id),
                salmon_fragments=meta["num_processed"],
                salmon_mapped=meta["num_mapped"],
                mapping_rate_pct=round(meta["percent_mapped"], 2),
                library_type=lib["expected_format"],
                compatible_fragment_ratio=round(lib["compatible_fragment_ratio"], 4),
                frag_length_mean=round(meta.get("frag_length_mean", float("nan")), 1),
                salmon_index=Path(json.loads((sq / "cmd_info.json").read_text())["index"]).parent.name,
            ))
    df = pd.DataFrame(rows)
    df.to_csv(out, sep="\t", index=False)
    print(df.drop(columns=["qb_id"]).to_string(index=False))


if __name__ == "__main__":
    main(Path(sys.argv[1]), [tuple(map(Path, a.split("=", 2))) for a in sys.argv[2:]])
