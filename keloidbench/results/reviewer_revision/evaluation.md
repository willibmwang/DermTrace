# Reproduced eight-study evaluation

Keloid is code **0** in the saved predictions. Full coverage is primary.

| Analysis | Studies | Profiles | Mean macro F1 | Worst fold |
|---|---:|---:|---:|---:|
| Primary, all studies | 8 | 101 | 0.8260 | 0.1818 |
| Post hoc sensitivity | 7 | 92 | 0.9180 | 0.7619 |

The primary result fails the existing 0.85 mean / 0.75 worst-fold gates. The subset is outcome-selected.

The historical prediction cross-check reversed class encoding. That cross-check is now corrected; the main historical score, calculated from the training metrics, is unchanged.

See `all_studies.csv` for per-study/sample metrics and `selective_by_study.csv` for coverage by class. Selective metrics use sample grain and are separate from the donor/sample primary metric.

Uncertainty is descriptive: eight studies and repeated development do not establish prospective generalization.
