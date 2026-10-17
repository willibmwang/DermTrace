# Reproducing the reviewer revision

Use Python 3.10. The published artifacts support offline reproduction; training and GPU generation are separate commands.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-revision.txt
python -m unittest discover -s tests -v
python scripts/revision_report.py
python scripts/run_marker_benchmark.py --report-only
python scripts/generate_revision_docs.py
python scripts/generate_revision_figures.py
python scripts/export_static_site_data.py
python scripts/replay_narrative_audit.py results/reviewer_revision/narrative_audit/deterministic_v1/responses.jsonl
```

`revision_report.py` recomputes the historical eight-study evaluation from saved predictions and verifies label encoding and the original mean. `--report-only` reproduces all matched comparison summaries from the saved per-profile probabilities. Figures have SVG, PNG, and PDF versions. The generated README block and manuscript results use the same JSON artifacts as the dashboard and browser export.

## Rerun the matched comparison

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python scripts/run_marker_benchmark.py
```

This small linear-model benchmark runs on CPU. It verifies frozen input checksums and uses all eight outer folds, six single markers, five fixed program views, and inner-selected marker/program families. `protocols/reviewer_revision_v1.json` records all splits, preprocessing, selection, weighting, calibration, thresholds, source hashes, and the limits of the fixed processed inputs. The protocol was specified during revision after historical outcomes were known. The new comparison changes the historical train/validation reservation consistently for every arm; the historical ensemble is a separate evidence object.

Inputs required to rerun the matched comparison are committed: `cohort_manifest.csv`, `marker_inputs.csv`, and the five feature-view Parquet files. Large raw expression matrices are not required. Rebuilding the frozen marker inputs from source requires the source matrices named and hashed in the protocol; their omission from Git is intentional. Earlier sourcing/preprocessing scripts and the cohort registry retain source provenance.

## Narrative experiments

The 101 frozen packets contain the live Stage-B predictions and derived evidence. Model generation never changes classification. The conservative automated checker is a screening tool, and independent human assessment remains necessary.

For an HTTP model endpoint, configure `KELOIDBENCH_LLM_API_URL`, `KELOIDBENCH_LLM_MODEL`, and optionally `KELOIDBENCH_LLM_API_KEY` outside Git. Then run:

```bash
python scripts/audit_copilot_faithfulness.py --arm llm --repeats 3 --output results/reviewer_revision/narrative_audit/new_endpoint_run
```

For the pinned local Qwen3-1.7B model, use a CUDA environment with PyTorch 2.6.0 and Transformers 5.2.0, and the snapshot `70d244cc86ccca08cf5af4e1e306ecf908b1ad5e`:

```bash
python scripts/audit_copilot_faithfulness.py --arm llm --repeats 3 \
  --local-model /path/to/pinned/snapshot --batch-size 16 \
  --output results/reviewer_revision/narrative_audit/new_gpu_run
```

The Yale Slurm helper is `scripts/run_narrative_audit_gpu.sbatch`; set `KB_PYTHON` and `KB_MODEL` for another environment. Model weights are cached outside the repository. `protocols/narrative_local_v1.json` fixes the model and greedy decoding settings before generation. Raw narratives, failed checks, request payloads, prompt hashes, model metadata, and repeat identifiers are retained. Repeated greedy outputs are not independent samples. Output directories cannot be overwritten.

Replay any completed run without loading a model:

```bash
python scripts/replay_narrative_audit.py results/reviewer_revision/narrative_audit/qwen3_1_7b_v1/responses.jsonl
```

Run and replay the v2 audit (`protocols/narrative_local_v2.json`; Qwen3-8B snapshot `b968826`; set `KB_MODEL` to its local path):

```bash
sbatch scripts/run_narrative_audit_v2_gpu.sbatch
python scripts/replay_narrative_audit.py results/reviewer_revision/narrative_audit/qwen3_8b_v2/responses.jsonl
```

Regenerate the post hoc failure diagnostic for the frozen local run. It is exploratory, written after inspecting outputs, and does not alter the run or validator:

```bash
python scripts/diagnose_narrative_audit.py results/reviewer_revision/narrative_audit/qwen3_1_7b_v1 --output /tmp/narrative_diagnostics
diff /tmp/narrative_diagnostics/summary.json results/reviewer_revision/narrative_audit/qwen3_1_7b_v1_diagnostics/summary.json
```

Each run supplies two blank rating forms. Follow the [human-review rubric](results/reviewer_revision/narrative_audit/README.md). Do not report unfilled forms as completed human evaluation.

## Publication and provenance

The current primary estimate is the eight-study historical ensemble. The seven-study result is an outcome-selected sensitivity analysis. The marker benchmark, live Stage-B model, and historical ensemble have different model artifacts and must not share performance claims.

The former GSE185309 lockbox was evaluated in an earlier experiment (`results/breakthrough_sprint/breakthrough_lockbox_once.json`). No lockbox is newly scored for this revision. Prior access prevents describing GSE185309 as untouched.

`website/` is the static GitHub Pages source. The browser makes deterministic summaries and scores uploaded matrices locally. The accepted poster is preserved as originally supplied; the current README and revision reports supersede its older headline claims. No new clinical validation or independent prospective validation is claimed.
