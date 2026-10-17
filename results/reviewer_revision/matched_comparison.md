# Matched marker/program comparison

Revision-time nested LOSO on all eight development studies. Each fixed arm uses the same classifier family, C grid, training-only preprocessing, inner-OOF calibration, and threshold selection. Historical ensemble results are a separate evidence object.

| Arm | Mean primary macro F1 | Worst fold | Mean sample macro F1 |
|---|---:|---:|---:|
| inner_selected_marker | 0.8295 | 0.1818 | 0.8187 |
| inner_selected_program | 0.8548 | 0.4156 | 0.8226 |
| marker_ACTA2 | 0.3558 | 0.1429 | 0.3617 |
| marker_COL1A1 | 0.9144 | 0.4444 | 0.8648 |
| marker_COL3A1 | 0.8608 | 0.1818 | 0.8339 |
| marker_CTHRC1 | 0.8087 | 0.1818 | 0.8180 |
| marker_FN1 | 0.8181 | 0.3250 | 0.8074 |
| marker_POSTN | 0.8260 | 0.3250 | 0.8235 |
| program_composition_only | 0.7503 | 0.3250 | 0.7286 |
| program_fibrosis_only | 0.8023 | 0.1818 | 0.7839 |
| program_fused_multiview | 0.8356 | 0.3250 | 0.7978 |
| program_rank_programs_only | 0.8965 | 0.4156 | 0.8642 |
| program_scar_discriminative | 0.8419 | 0.4156 | 0.8041 |

Inner-selected program minus inner-selected marker: **+0.0253**; descriptive paired study bootstrap interval **[-0.1081, +0.1480]**.

No universal superiority claim is made. The development cohorts were already inspected; these results require independent confirmation. Missingness and per-study predictions are published alongside this report.
