# Tsao lab Greek melanoma FFPE RNA-seq (`tsaolab_melanoma`)

**Status (2026-10-29): processing in progress.** Expression matrices, QC and variant tables are added as they finish.

Bulk RNA-seq of 48 formalin-fixed, paraffin-embedded (FFPE) melanoma specimens from the Tsao lab's Greek cohort,
sequenced as core project `20088-11` (sample IDs `20088XR-11-01` … `-50`; 45 and 49 were not delivered).
Paired-end 2×151 on an Illumina NovaSeq X (`LH00150`), one lane per sample, spread over runs 1058 and 1062 and
lanes L001/L004/L006/L008. The raw FASTQs (about 259 GB, 96 files) are on the lab Google Drive in
`20088-11 greek melanoma fastq.gz` with the core's `md5.txt` — ask William for access.

**No sample key came with the delivery.** `sample_sheets/samples_20088-11.csv` therefore carries only what the file
names say (core ID, lane, index). Patient, primary site, stage, treatment and the pathology BRAF/NRAS status have to be
joined on `core_id` from the lab's own list before any clinical analysis.

## What is being produced (`scripts/tsaolab_melanoma/`)

1. `make_sample_sheet.py` — one row per sample from `md5.txt`.
2. `run_quant.sh` — per sample: `fastp` (adapter auto-detection, reads shorter than 25 nt dropped; FFPE inserts peak
   near 120 nt, so most pairs are adapter-trimmed), `salmon quant` 1.10.3 against GENCODE v46
   (`-l A --validateMappings --gcBias --seqBias`), and `minimap2 -x sr` of the same trimmed reads to a
   **melanoma driver transcript reference** (`qc/driver_ref.tsv`: MANE Select transcript of 46 driver genes, 34
   paralog decoys and all their pseudogenes, 0.47 Mb).
3. `variants.py` — joint `bcftools mpileup` over those alignments, allele counts per sample, codon-level annotation,
   hotspot calls (BRAF 600/601, NRAS/KRAS/HRAS 12/13/61, KIT, GNAQ/GNA11 209, MAP2K1, RAC1 P29, CTNNB1, IDH1 R132,
   SF3B1 R625, …) and any coding change in the tumour suppressors (NF1, TP53, CDKN2A, PTEN, ARID2, RB1, BAP1, …).
4. `run_qc_summary.py`, `build_dermtrace.py` — read-level QC and the DermTrace-format TPM/count matrices and metadata.
5. Donor check with `scripts/tsaolab_rnaseq/genotype_check.sh` + `genotype.py` (1,500 well-expressed transcripts,
   2M read pairs per sample) to find specimens that come from the same patient.

The same approach gave the SNP genotypes used to re-identify the oxygen-tension cell lines
(`data/processed/tsaolab_o2/qc/genotype_match.tsv`); for FFPE tumours the point is somatic hotspots, which RNA can
show only where the gene is expressed. Promoter mutations (TERT) are invisible, and low-fraction C>T/G>A and A>G
calls are flagged as possible deamination / RNA editing rather than reported as somatic.
