# Archived processed data

These files were once committed under `data/processed/`, then removed from git. They were restored byte for byte from the last commit that contained them. Each subfolder is named after the date the files were removed and keeps their original layout.

| Folder | Removed in | Files |
|---|---|---|
| `2026-07-09/` | `e625f2f` ("stop tracking local data matrices") | `bulk_rnaseq/`, `microarray/` (including `signatures/`), `scrna_spatial/` and `spheroid/` matrices and tables |
| `2026-06-11/` | `f7e2dca` ("mvp rework") | `choi_extraction_summary.csv`, an earlier version of `spheroid/choi_extraction_summary.csv`. The other files removed that day were moved unchanged into `spheroid/` and `microarray/signatures/`, so they appear in `2026-07-09/` or in the current dataset folders. |

**These are not the versions used in the KeloidBench paper.** The bulk, microarray and single-cell matrices were rebuilt after July 9. Their SHA-256 hashes differ from those frozen in [`keloidbench/protocols/reviewer_revision_v1.json`](../../../keloidbench/protocols/reviewer_revision_v1.json). They are kept for provenance only. Use the current files in the parent dataset folders for analysis.
