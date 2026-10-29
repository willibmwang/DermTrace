"""One row per sample of sequencing / trimming / quantification / driver-alignment QC for the 20088-11 run.

Usage: python -I run_qc_summary.py <out.tsv> <sample_sheet.csv> <quant_out_dir> <raw_data_dir>
"""
import gzip
import json
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

INSTRUMENTS = {"LH": "Illumina_NovaSeqX", "A0": "Illumina_NovaSeq", "NB": "Illumina_NextSeq", "AV": "Element_AVITI"}
SAMTOOLS = Path(__import__("os").environ.get("VARCALL_ENV", "/opt/anaconda3/envs/varcall")) / "bin" / "samtools"


def main(out, sheet, qdir, raw):
    rows = []
    for r in pd.read_csv(sheet, dtype=str).itertuples():
        sq = qdir / "quant" / r.sample_id
        if not (sq / "quant.sf").exists():
            continue
        with gzip.open(raw / r.r1, "rt") as fh:
            inst, run, flowcell, lane = fh.readline()[1:].split(":")[:4]
        fp = json.loads((qdir / "fastp" / f"{r.sample_id}.json").read_text())
        s = fp["summary"]
        meta = json.loads((sq / "aux_info" / "meta_info.json").read_text())
        lib = json.loads((sq / "lib_format_counts.json").read_text())
        bam = qdir / "driver_bam" / f"{r.sample_id}.bam"
        driver_reads = int(subprocess.run([SAMTOOLS, "view", "-c", bam], capture_output=True, text=True).stdout.strip()) if bam.exists() else 0
        rows.append(dict(
            sample_id=r.sample_id, core_id=r.core_id, instrument=INSTRUMENTS.get(inst[:2], inst),
            sequencing_run=f"{inst}:{run}", flowcell=flowcell, lane=f"L{int(lane):03d}",
            reads_raw=s["before_filtering"]["total_reads"], reads_after_trim=s["after_filtering"]["total_reads"],
            pct_reads_kept=round(s["after_filtering"]["total_reads"] / s["before_filtering"]["total_reads"] * 100, 2),
            pct_adapter_trimmed=round(fp.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0) / s["before_filtering"]["total_reads"] * 100, 2),
            mean_read_length_after_trim=round((s["after_filtering"]["read1_mean_length"] + s["after_filtering"]["read2_mean_length"]) / 2, 1),
            insert_size_peak=fp.get("insert_size", {}).get("peak", float("nan")),
            pct_q30_after_trim=round(s["after_filtering"]["q30_rate"] * 100, 2),
            fastp_duplication_rate=round(fp.get("duplication", {}).get("rate", float("nan")), 4),
            salmon_fragments=meta["num_processed"], salmon_mapped=meta["num_mapped"],
            mapping_rate_pct=round(meta["percent_mapped"], 2), library_type=lib["expected_format"],
            compatible_fragment_ratio=round(lib["compatible_fragment_ratio"], 4),
            frag_length_mean=round(meta.get("frag_length_mean", float("nan")), 1),
            driver_ref_reads=driver_reads,
        ))
    df = pd.DataFrame(rows)
    df.to_csv(out, sep="\t", index=False)
    print(df.drop(columns=["core_id", "flowcell"]).to_string(index=False))


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4]))
