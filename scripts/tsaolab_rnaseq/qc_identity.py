"""Per-sample cell-identity QC: sex-chromosome expression and HEK293-like admixture.

Combines the adenovirus-5 screen (screen_ad5.sh) with gene TPMs from salmon and writes one row per sample:

  xist_tpm, chrY_tpm           XIST and summed chrY genes (human samples only)
  ad5_e1_per_million           Ad5 E1A/E1B read pairs per million pairs screened (HEK293 signature)
  hek_admixture_est_pct        E1 rate / median E1 rate of the HEK293-like samples (non-HEK human samples)
  sex_inferred                 male | female | no_Y_no_XIST (male with Y loss, or XIST-eroded female) | mixed
  qc_flag, qc_note             "suspect" + reason, or "ok"

Usage: python -I qc_identity.py <out.csv> <sample_sheet.csv>=<quant_dir>=<ad5_screen.tsv> [...]
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python -I does not add the script dir
from aggregate import FASTA, tx2gene  # noqa: E402

CHRY = ["RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y", "ZFY"]
ADMIX_SUSPECT_PCT = 2.0  # flag non-HEK samples with more HEK-like material than this


def main(out, runs):
    t2sym = tx2gene(FASTA["human"]).set_index("tx").gene_symbol
    rows = []
    for sheet, qdir, screen in runs:
        s = pd.read_csv(sheet, dtype=str)
        ad5 = pd.read_csv(screen, sep="\t", dtype={"sample_id": str}).set_index("sample_id")
        for r in s.itertuples():
            row = dict(sample_id=r.sample_id, qb_id=r.qb_id, species=r.species, cell_line=r.cell_line)
            a = ad5.loc[r.sample_id]
            row["ad5_e1_per_million"] = round(a.ad5_e1_pairs / a.pairs_screened * 1e6, 1)
            qf = Path(qdir) / "quant" / r.sample_id / "quant.sf"
            if r.species == "human" and qf.exists():
                q = pd.read_csv(qf, sep="\t", index_col=0)
                g = q.TPM.groupby(t2sym.loc[q.index].values).sum()
                row["xist_tpm"] = round(g.get("XIST", 0.0), 1)
                row["chrY_tpm"] = round(sum(g.get(k, 0.0) for k in CHRY), 1)
            rows.append(row)
    df = pd.DataFrame(rows)

    hek = df.cell_line == "HEK293-like"
    hek_rate = df.loc[hek, "ad5_e1_per_million"].median()
    human = df.species == "human"
    df.loc[human & ~hek, "hek_admixture_est_pct"] = (df.ad5_e1_per_million / hek_rate * 100).round(1)

    def sex(r):
        if pd.isna(r.get("xist_tpm")):
            return "not_assessed"
        x, y = r.xist_tpm >= 10, r.chrY_tpm >= 20
        return "mixed" if x and y else "female" if x else "male" if y else "no_Y_no_XIST"
    df["sex_inferred"] = df.apply(sex, axis=1)

    # A cell line should be uniformly chrY+ or chrY-. XIST is unreliable in cultured lines (erosion), so
    # line consistency is judged on chrY only, and a split line is noted on every one of its samples.
    df["chrY_present"] = df.chrY_tpm >= 20
    has_y = df[df.chrY_tpm.notna()].groupby("cell_line").chrY_present.agg(["sum", "count"])
    split = has_y[(has_y["sum"] > 0) & (has_y["sum"] < has_y["count"])]
    # trace male-cell admixture in a chrY- line, scaled to clean male BJ (chrY+, no Ad5)
    male_ref = df[(df.cell_line == "BJ") & (df.ad5_e1_per_million < 50)].chrY_tpm.median()

    def flag(r):
        notes = []
        if r.cell_line == "HEK293-like":
            notes.append("submitted as NIH3T3 (mouse); RNA is human, female, Ad5 E1+ (HEK293-like)")
        if r.get("hek_admixture_est_pct", 0) > ADMIX_SUSPECT_PCT:
            notes.append(f"~{r.hek_admixture_est_pct:.0f}% HEK293-like admixture (Ad5 E1 reads, XIST in male line)")
        if r.cell_line in split.index:
            s = split.loc[r.cell_line]
            notes.append(f"line split by chrY expression ({int(s['sum'])} chrY+ / {int(s['count'] - s['sum'])} chrY- "
                         f"samples), this sample chrY{'+' if r.chrY_present else '-'}: Y loss or different stock")
        elif r.cell_line in has_y.index and has_y.loc[r.cell_line, "sum"] == 0 and r.chrY_tpm >= 5:
            notes.append(f"trace chrY in chrY- line ({r.chrY_tpm:.1f} TPM, ~{r.chrY_tpm / male_ref * 100:.1f}% "
                         "male-cell admixture)")
        suspect = r.get("hek_admixture_est_pct", 0) > ADMIX_SUSPECT_PCT
        return pd.Series(dict(qc_flag="suspect" if suspect else "ok", qc_note="; ".join(notes)))
    df[["qc_flag", "qc_note"]] = df.apply(flag, axis=1)

    df.to_csv(out, index=False)
    print(f"HEK293-like Ad5 E1 median: {hek_rate:.0f} per million pairs")
    print(df.drop(columns=["species"]).to_string(index=False))


if __name__ == "__main__":
    main(Path(sys.argv[1]), [tuple(a.split("=", 2)) for a in sys.argv[2:]])
