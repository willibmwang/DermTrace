# r3 candidate: stricter validated-pairs parser

## Status

| Version | What it is | Qualification |
|---|---|---|
| **r2** (`rqc-src-20261124-r2-08dd537f`) | Frozen quantification code | **Real-tool qualified on Rivanna** for GM01: job 23402243, all 14 provisional-v1 gates PASS. See `QUALIFICATION_GM01_r2.md`. Its original text is preserved, with a supersession note that separates the r2 deployment identity from the r3 repository identity |
| **r3** (`rqc-src-20261125-r3-4432a801`) | r2 plus a stricter `read_validated_pairs` (same rules as `run_fastqc_all.sbatch`) and its tests | **Parser-specific qualification PASSED on Rivanna** (user-reported, 2026-11-24). No new GM01 quantification was needed, for the reasons in "Effect on quantification". See `QUALIFICATION_STATUS.md` |

## What changed from r2 to r3

- **Changed in `rivanna_qc_tools.py`:** only `read_validated_pairs`, plus
  its new helper `vp_is_header` and the constants `VP_*`.
  - Comparing the two files definition by definition, the other 48
    top-level definitions are identical.
  - The parser is called in exactly one place, `cmd_resolve`
    (`resolve-specimen`).
- **Changed in `tests/run_local_tests.sh`:**
  - One fixture now moves both mates together. r3 requires R1 and R2 in
    the same folder.
  - A new 48-specimen eligibility section was added.
- **Changed in `README.md`:** one paragraph on the r3 manifest rule.
- **Unchanged:** `quantify_specimen.sbatch`, `lib_rivanna_qc.sh`, the
  other scripts and `config.example.env`.
- **Added later, documentation only:**
  - a supersession note in `QUALIFICATION_GM01_r2.md`; its original text is
    preserved;
  - `QUALIFICATION_STATUS.md`.

## How the parsers behave differently

Both parsers were run on the same manifest variants: the historical
layout, 43 rows from 11-01 to 11-42 plus 11-50. "Rows" are the
`(label, R1, R2)` tuples returned.

| Manifest variant | r2 | r3 |
|---|---|---|
| Historical: header `sample/r1_path/r2_path`, core-ID labels, absolute paths, LF | accept, 43 rows | accept, **same 43 rows** |
| Headerless; header `sample/r1/r2`; `SAMPLE/R1/R2`; labels `11-NN` or `GMNN`; no final newline | accept | accept, same rows |
| Header whose first word is not `sample`, `sample_id`, `specimen` or `label` | accept (line 1 skipped) | **reject** |
| A 4th column | accept (ignored) | **reject** |
| Relative paths (basenames only) | accept | **reject** (paths must be absolute) |
| Label padded with spaces | accept (spaces stripped) | **reject** |
| R2 in another folder, not R1's mate name, or R1 and R2 swapped | accept (the specimen is then blocked by the sheet match) | **reject** |
| Duplicate label, duplicate row, or the same specimen listed twice | accept | **reject** |
| Label naming a different specimen (`11-09` on 11-01's files) | accept (only that specimen is blocked) | **reject** |
| CRLF line endings, blank lines, padded paths, `.fq.gz`, rows with 2 fields | reject | reject |

**Two structural differences:**

1. **r3 rejects the whole manifest.** In r2, a bad row affected only
   that specimen. In r3, any defect anywhere blocks every specimen.
2. **r3 returns paths exactly as written.** r2 stripped surrounding
   whitespace; r3 refuses it instead. For valid manifests the output
   is identical.

## Eligibility decisions

Every specimen in the committed sample sheet (48) was run through the full
`resolve-specimen` command under both versions. The manifest replica used the
real filenames from the sheet.

| Manifest | Eligible under r2 | Eligible under r3 | Decisions changed | Resolved output where both eligible |
|---|---|---|---|---|
| Historical layout, core-ID labels | 42 | **42** | 0 | byte-identical |
| Headerless | 42 | 42 | 0 | byte-identical |
| `11-NN` labels | 42 | 42 | 0 | byte-identical |
| `GMNN` labels | 42 | 42 | 0 | byte-identical |
| One padded label (row for 11-21) | 42 | **0** | 42 | n/a |
| 11-37 listed twice | 42 | **0** | 42 | n/a |

In every clean layout, GM37 (11-37) and the corrupt 11-43, 11-44, 11-46,
11-47 and 11-48 are blocked by the built-in exclusions. That happens before
the manifest is read, so it is the same in r2 and r3. A corrupt specimen
whose exclusion is lifted is still blocked: it is not a validated pair.

## Effect on quantification

- `quantify_specimen.sbatch` uses only the 8 tab-separated fields that
  `resolve-specimen` prints:
  - sample ID and core ID;
  - R1 and R2 file names (taken from the **sample sheet**, then joined to
    `RAW_DIR`);
  - MD5s (from the **sample sheet**);
  - the validated-pairs label;
  - any lifted exclusions.
- The manifest only decides eligibility, supplies the label, and is used
  to check that its absolute paths are the same files as `RAW_DIR`.
- So when r2 and r3 make the same eligibility decision and print the same
  output, the inputs, fastp and Salmon arguments, index and container are
  all identical.
- Provenance does record a different `tools_sha256`, so r3 runs remain
  distinguishable from r2 runs.
- **Therefore GM01 is not re-quantified for r3.** This conclusion depends
  on the real Rivanna manifest giving identical `resolve-specimen` output
  under both versions, which qualification step 2 below checks directly.

## Local test results (developer machine, stub tools)

All results below are on the integrated r3 tree: the r2 files, the
legacy files from 412ad24, and the changes above.

| Suite | Python 3.9.6 | Python 3.11.15 |
|---|---|---|
| `tests/run_local_tests.sh` (quantification) | 188 passed, 0 failed | 188 passed, 0 failed |
| `tests/test_legacy_qc_scripts.py` (legacy QC, including `ParserParity`) | 48 OK | 48 OK |
| Legacy suite under a simulated Lmod login environment (`BASH_ENV`, an exported `module`) | 48 OK | — |

For comparison, the r2 tree fails 16 `ParserParity` sub-tests.

## Required r3 qualification on Rivanna (parser-specific)

1. Both suites pass with Rivanna's Python.
2. Using the **real** `validated_pairs.tsv`, sample sheet and `RAW_DIR`
   (with `--raw-dir`), every sheet specimen is resolved under both r2 and
   r3. To pass:
   - the r3 parser accepts the manifest;
   - each specimen gets the same exit code under both;
   - eligible specimens get byte-identical output under both;
   - exactly 42 are eligible;
   - exactly GM37, GM43, GM44, GM46, GM47 and GM48 are blocked.
3. Nothing is quantified. The r2 GM01 qualification and the cross-specimen
   runs (made with r2) stand as they are.

Step 2 passed on Rivanna; see `QUALIFICATION_STATUS.md`. The executable and
test files in this repository are byte-identical to the qualified r3 bundle.
Only the documentation files `R3_PARSER_CANDIDATE.md`, `README.md` and
`QUALIFICATION_STATUS.md` were updated afterwards, to record the result.
