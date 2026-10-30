"""Processed-data completion adapters and explicit information boundaries.

No source tables are written. Source and sample grouping are caller supplied;
no assertion that numerical normalization harmonizes different assay types.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
from typing import Callable, Iterator

import numpy as np
import pandas as pd
from .gse163973_reader import SAMPLE_ALIASES as ALIASES, iter_matrices

Verify = Callable[[Path, str], None]
PARTITIONS = frozenset({'train', 'validation', 'test', 'quarantine', 'auxiliary'})


class SplitAssignments:
    """Whole biological groups have one role, even across source samples."""
    def __init__(self, rows):
        self.rows = {}
        self.metadata = {}
        group_role = {}
        for row in rows:
            source, sample, group, role = (row.get(k) for k in ('source_id', 'sample_id', 'group_id', 'partition'))
            if not all(isinstance(x, str) and x.strip() for x in (source, sample, group, role)) or role not in PARTITIONS:
                raise ValueError('Invalid split assignment')
            key = (source, sample)
            if key in self.rows:
                raise ValueError('Duplicate source sample assignment')
            if group in group_role and group_role[group] != role:
                raise ValueError('Biological group crosses partitions')
            self.rows[key] = (group, role)
            self.metadata[key] = dict(row)
            group_role[group] = role
        if not self.rows:
            raise ValueError('Empty split assignments')

    @classmethod
    def read(cls, path, *, metadata_paths=None):
        """Read a pinned split and verify each supplied metadata input.

        Metadata names in the manifest are portable labels, not local paths.
        Callers must explicitly bind each name to the corresponding file.
        """
        path = Path(path)
        manifest = json.loads((path.parent/'manifest.json').read_text())
        matches = [x for x in manifest['files'] if x['file'] == path.name]
        if len(matches) != 1:
            raise ValueError('Split is absent or duplicated in its companion manifest')
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != matches[0]['sha256']:
            raise ValueError('Split file checksum mismatch')
        payload = json.loads(raw)
        rows = payload.get('assignments')
        if not isinstance(rows, list) or not rows:
            raise ValueError('Split must contain sample assignments')
        source_metadata = {'GSE163973': 'skin', 'GSE144424': 'cardiac',
                           'GSE95280_recount3': 'melanoma'}
        sources = {row.get('source_id') for row in rows if isinstance(row, dict)}
        if len(sources) == 0 or len([r for r in rows if isinstance(r, dict)]) != len(rows) or not sources <= source_metadata.keys():
            raise ValueError('Unknown source in split assignment')
        required_metadata = {source_metadata[source] for source in sources}
        global_metadata = manifest['input_metadata_sha256']
        expected_metadata = matches[0].get('input_metadata_sha256', global_metadata)
        if not isinstance(expected_metadata, dict) or not isinstance(global_metadata, dict):
            raise ValueError('Split metadata hashes must be a mapping')
        if not expected_metadata or any(global_metadata.get(k) != v for k, v in expected_metadata.items()):
            raise ValueError('Split metadata hashes disagree with companion manifest')
        declared_metadata = payload.get('input_metadata_sha256')
        if 'input_metadata_sha256' in matches[0]:
            if declared_metadata != expected_metadata or set(expected_metadata) != required_metadata:
                raise ValueError('Split metadata hashes do not cover its assigned sources')
        elif (declared_metadata is not None and declared_metadata != expected_metadata) or not required_metadata <= expected_metadata.keys():
            raise ValueError('Legacy split metadata hashes do not cover its assigned sources')
        paths = {} if metadata_paths is None else dict(metadata_paths)
        if set(paths) != set(expected_metadata):
            raise ValueError('Provide a path for each pinned split metadata input')
        for name, expected in expected_metadata.items():
            with Path(paths[name]).open('rb') as stream:
                if hashlib.file_digest(stream, 'sha256').hexdigest() != expected:
                    raise ValueError('Split source metadata checksum mismatch: '+name)
        return cls(rows)

    def get(self, source, sample):
        try:
            return self.rows[(source, sample)]
        except KeyError as e:
            raise ValueError(f'Missing split for {source} / {sample}') from e

    def check_skin_pairs(self):
        for suffix in '123':
            left = self.rows.get(('GSE163973', 'KF'+suffix))
            right = self.rows.get(('GSE163973', 'NF'+suffix))
            if left is not None and right is not None and left != right:
                raise ValueError('Conservative paired skin samples must share group and partition')


@dataclass(frozen=True)
class Profile:
    uid: str
    source_id: str
    sample_id: str
    group_id: str
    partition: str
    unit: str
    genes: tuple[str, ...]
    values: np.ndarray
    measured: np.ndarray
    cell_type: str | None = None
    audit_flags: tuple[str, ...] = ()

    def __post_init__(self):
        n = len(self.genes)
        if not self.uid or not self.group_id or self.partition not in PARTITIONS:
            raise ValueError('Profile identity, group and partition are required')
        if not n or len(set(self.genes)) != n or any(not x for x in self.genes):
            raise ValueError('Gene IDs must be unique and nonempty')
        x = np.asarray(self.values)
        measured = np.asarray(self.measured)
        if x.shape != (n,) or measured.shape != (n,) or measured.dtype != np.bool_:
            raise ValueError('Profile expression and availability axes do not agree')
        if not np.isfinite(x).all() or (x[measured] < 0).any() or (x[~measured] != 0).any():
            raise ValueError('Use finite nonnegative values; missing genes must have zero fill')


def align(profile, panel):
    """Exact Ensembl join only: unmeasured is separate from measured zero."""
    panel = tuple(panel)
    if not panel or len(set(panel)) != len(panel) or any(not isinstance(g, str) or not g.startswith('ENSG') for g in panel):
        raise ValueError('Expected a unique human Ensembl gene panel')
    index = pd.Index(profile.genes).get_indexer(panel)
    found = index >= 0
    out = np.zeros(len(panel), dtype=float)
    measured = np.zeros(len(panel), dtype=bool)
    out[found] = profile.values[index[found]]
    measured[found] = profile.measured[index[found]]
    return Profile(profile.uid, profile.source_id, profile.sample_id, profile.group_id,
                   profile.partition, profile.unit, panel, out, measured, profile.cell_type, profile.audit_flags)


class AuthorAnnotatedCellPolicy:
    """Keep source-annotated cells; retain descriptive QC flags for audit only."""
    def __init__(self, qc_path, expected_sha256):
        with Path(qc_path).open('rb') as stream:
            if hashlib.file_digest(stream, 'sha256').hexdigest() != expected_sha256:
                raise ValueError('Cell QC metadata checksum mismatch')
        qc = pd.read_parquet(qc_path)
        required = {'sample', 'barcode', 'annotated', 'cell_type', 'flag_low_counts',
                    'flag_low_features', 'flag_high_mito', 'flag_any'}
        if not required <= set(qc.columns):
            raise ValueError('QC metadata lacks required fields')
        self.qc = qc.set_index(['sample', 'barcode'], verify_integrity=True)

    def __call__(self, sample, barcode, row):
        try:
            qc = self.qc.loc[(sample, barcode)]
        except KeyError as e:
            raise ValueError('Cell is missing from QC metadata') from e
        if pd.isna(qc.annotated) or bool(qc.annotated) != pd.notna(row.source_metadata_row):
            raise ValueError('QC and source annotation status disagree')
        if str(qc.cell_type) != str(row.cellType):
            raise ValueError('QC and source cell-type labels disagree')
        return bool(qc.annotated)

    def flags(self, sample, barcode):
        qc = self.qc.loc[(sample, barcode)]
        names = ('flag_low_counts', 'flag_low_features', 'flag_high_mito')
        if any(pd.isna(qc[name]) for name in names):
            raise ValueError('Annotated cell has missing QC flags')
        return tuple(name for name in names if bool(qc[name]))


def cardiac_profiles(root, split, verify: Verify, primary_only=True, include_flagged=False) -> Iterator[Profile]:
    root = Path(root)
    manifest = json.loads((root/'manifest.json').read_text())
    for name in ('counts.parquet', 'samples.csv'):
        verify(root/name, manifest['outputs'][name])
    counts = pd.read_parquet(root/'counts.parquet')
    meta = pd.read_csv(root/'samples.csv').set_index('sample_id', verify_integrity=True)
    if not counts.index.is_unique or not counts.columns.is_unique or set(counts.columns) != set(meta.index):
        raise ValueError('Cardiac genes or sample join is inconsistent')
    # Metadata's donor relationship must agree with caller-provided grouping.
    donor_groups = {}
    for sample in counts.columns:
        group, _ = split.get('GSE144424', str(sample))
        donor = meta.loc[sample, 'donor_id']
        if pd.isna(donor) or not str(donor):
            raise ValueError('Unknown cardiac donor')
        if donor in donor_groups and donor_groups[donor] != group:
            raise ValueError('Cardiac donor crosses groups')
        donor_groups[donor] = group
    ids = tuple(map(str, counts.index))
    if any(not x.startswith('ENSG') for x in ids):
        raise ValueError('Unexpected cardiac feature identifiers')
    for sample in counts.columns:
        eligible = meta.loc[sample, 'published_primary_eligible']
        if primary_only and (pd.isna(eligible) or eligible != True):
            continue
        if not include_flagged and meta.loc[sample, 'quality_flag'] != 'no_specific_paper_flag':
            continue
        group, role = split.get('GSE144424', str(sample))
        yield Profile('GSE144424:'+str(sample), 'GSE144424', str(sample), group, role,
                      'bulk_gene_featurecounts', ids, counts[sample].to_numpy(float),
                      np.ones(len(ids), dtype=bool), 'cardiomyocyte')


def skin_pseudobulk_profiles(root, split, verify: Verify, samples=None) -> Iterator[Profile]:
    root = Path(root)
    split.check_skin_pairs()
    m = json.loads((root/'recovery_manifest.json').read_text())
    for name in ('pseudobulk_sum_counts.parquet', 'metadata.parquet', 'gene_annotations.parquet'):
        verify(root/name, m['output_sha256'][name])
    counts = pd.read_parquet(root/'pseudobulk_sum_counts.parquet')
    meta = pd.read_parquet(root/'metadata.parquet').set_index('sample_id', verify_integrity=True)
    gene = pd.read_parquet(root/'gene_annotations.parquet')
    if not counts.index.is_unique or not counts.columns.is_unique or set(counts.index) != set(meta.index):
        raise ValueError('Pseudobulk sample join is inconsistent')
    selected = set(meta.source_donor_key.astype(str)) if samples is None else set(samples)
    if not selected or not selected <= set(meta.source_donor_key.astype(str)):
        raise ValueError('Unknown or empty skin sample selection')
    for sample in selected:
        split.get('GSE163973', sample)
    # Exclude aggregated symbols: these values sum multiple distinct source IDs.
    unique = gene[~gene['symbol_aggregates_multiple_features'].astype(bool)].copy()
    unique = unique[unique['source_gene_id'].astype(str).str.startswith('ENSG')]
    if unique.source_gene_id.duplicated().any() or unique.gene_symbol.duplicated().any():
        raise ValueError('Unique source mapping is not one-to-one')
    unique = unique[unique.gene_symbol.isin(counts.columns)]
    symbols, ids = list(unique.gene_symbol), tuple(unique.source_gene_id.astype(str))
    for uid, row in meta.iterrows():
        sample = str(row['source_donor_key'])
        if sample not in selected:
            continue
        group, role = split.get('GSE163973', sample)
        yield Profile(str(uid), 'GSE163973', sample, group, role, 'pseudobulk_umi_counts', ids,
                      counts.loc[uid, symbols].to_numpy(float), np.ones(len(ids), dtype=bool), str(row['cell_type']))


def melanoma_profiles(root, split, verify: Verify, read_annotations, clean_only=True) -> Iterator[Profile]:
    """Read recount3 coverage with exact, source-supported, one-to-one IDs.

    Ambiguous aggregates and repeated Ensembl IDs are not guessed or summed.
    Default eligibility uses the source's clean oxygen annotation; lineage
    grouping must be supplied in the split assignments.
    """
    root = Path(root)
    m = json.loads((root/'recount_source_manifest.json').read_text())
    for name in ('expression_native_symbol.parquet', 'metadata.parquet'):
        verify(root/name, m['output_sha256'][name])
    path = root/'expression_native_symbol.parquet'
    anno = read_annotations(path)
    if anno['gene_axis'] != 'columns' or anno['species'] != 'human':
        raise ValueError('Unexpected melanoma gene axis or species')
    candidates = [g for g in anno['genes'] if g['mapping_basis'] == 'source_annotation'
                  and len(g['ensembl_gene_ids']) == 1 and len(g['source_ensembl_ids']) == 1]
    multiplicity = pd.Series([g['ensembl_gene_ids'][0] for g in candidates]).value_counts()
    candidates = [g for g in candidates if multiplicity[g['ensembl_gene_ids'][0]] == 1]
    labels = [g['native_label'] for g in candidates]
    ids = tuple(g['ensembl_gene_ids'][0] for g in candidates)
    counts = pd.read_parquet(path)
    meta = pd.read_parquet(root/'metadata.parquet').set_index('sample_id', verify_integrity=True)
    if not counts.index.is_unique or not counts.columns.is_unique or set(counts.index) != set(meta.index):
        raise ValueError('Melanoma sample join is inconsistent')
    if not set(labels) <= set(counts.columns) or not ids:
        raise ValueError('Melanoma gene annotations do not match expression columns')
    families = {}
    for uid in counts.index:
        row = meta.loc[uid]
        eligible = row['clean_hypoxia_eligible']
        if clean_only and (pd.isna(eligible) or eligible != True):
            continue
        group, _ = split.get('GSE95280_recount3', str(uid))
        family = row['cell_line_family']
        if pd.isna(family) or not str(family):
            raise ValueError('Unknown melanoma lineage family')
        if family in families and families[family] != group:
            raise ValueError('Melanoma lineage crosses groups')
        families[family] = group
    for uid in counts.index:
        row = meta.loc[uid]
        if clean_only and (pd.isna(row['clean_hypoxia_eligible']) or row['clean_hypoxia_eligible'] != True):
            continue
        group, role = split.get('GSE95280_recount3', str(uid))
        yield Profile('GSE95280_recount3:'+str(uid), 'GSE95280_recount3', str(uid), group, role,
                      'recount3_gene_base_pair_coverage', ids, counts.loc[uid, labels].to_numpy(float),
                      np.ones(len(ids), dtype=bool), 'melanoma')


def skin_cells(root, split, verify: Verify, cell_eligible, samples=None) -> Iterator[Profile]:
    """Stream a single compressed sample matrix at a time; annotated barcodes only.

    A selected sample can be supplied to bound memory and startup time. The
    required callback applies the caller's declared per-cell QC selection.
    """
    root = Path(root)
    if not callable(cell_eligible):
        raise ValueError('A declared per-cell eligibility policy is required')
    split.check_skin_pairs()
    m = json.loads((root/'recovery_manifest.json').read_text())
    archive = root/'downloads/GSE163973_RAW.tar'
    verify(archive, m['source_files'][archive.name]['sha256'])
    for name in ('cell_barcode_crosswalk.parquet', 'gene_annotations.parquet'):
        verify(root/name, m['output_sha256'][name])
    cross = pd.read_parquet(root/'cell_barcode_crosswalk.parquet')
    annot = pd.read_parquet(root/'gene_annotations.parquet')
    expected_ids = annot.source_gene_id.astype(str).tolist()
    if len(expected_ids) != len(set(expected_ids)):
        raise ValueError('Source Ensembl IDs are not unique')
    selected = set(ALIASES.values()) if samples is None else set(samples)
    if not selected or not selected <= set(ALIASES.values()):
        raise ValueError('Unknown or empty source sample selection')
    for sample, features, barcodes, matrix in iter_matrices(archive, selected):
        group, role = split.get('GSE163973', sample)
        if features.iloc[:, 0].tolist() != expected_ids:
            raise ValueError('Cell matrix axes do not match source annotations')
        source_cross = cross[cross.source_donor_key.eq(sample)].set_index('matrix_barcode', verify_integrity=True)
        if set(source_cross.index) != set(barcodes):
            raise ValueError('Cell barcode crosswalk is inconsistent')
        ids = tuple(expected_ids)
        for j, barcode in enumerate(barcodes):
            row = source_cross.loc[barcode]
            if pd.isna(row.source_metadata_row):
                continue
            decision = cell_eligible(sample, barcode, row)
            if type(decision) not in (bool, np.bool_):
                raise ValueError('Cell eligibility policy must return a boolean')
            if not decision:
                continue
            flags = cell_eligible.flags(sample, barcode) if hasattr(cell_eligible, 'flags') else ()
            yield Profile('GSE163973:'+sample+':'+barcode, 'GSE163973', sample, group, role,
                          'single_cell_umi_counts', ids, matrix.getcol(j).toarray().ravel().astype(float),
                          np.ones(len(ids), dtype=bool), str(row.cellType), tuple(flags))
        del matrix


@dataclass(frozen=True)
class MaskPlan:
    context: np.ndarray
    reward: np.ndarray
    evaluation: np.ndarray
    anchors: np.ndarray

    def validate(self, measured):
        measured = np.asarray(measured)
        masks = [np.asarray(v) for v in (self.context, self.reward, self.evaluation, self.anchors)]
        if measured.dtype != np.bool_ or any(v.dtype != np.bool_ or v.shape != measured.shape for v in masks):
            raise ValueError('Masks must be boolean and match measured genes')
        c, r, e, a = masks
        if not c.any() or not r.any() or not e.any() or not a.any():
            raise ValueError('Context, reward, evaluation and anchors must be nonempty')
        if ((c & r) | (c & e) | (r & e)).any() or ((c | r | e) & ~measured).any() or (a & ~c).any():
            raise ValueError('Masks overlap, include missing genes, or anchors are outside context')


def select_training_anchors(training_profiles, panel, count):
    """Select anchors by donor-weighted detection using training profiles only."""
    panel = tuple(panel)
    if not isinstance(count, int) or isinstance(count, bool) or count < 1 or count > len(panel):
        raise ValueError('Invalid requested anchor count')
    if len(set(panel)) != len(panel):
        raise ValueError('Panel contains duplicate genes')
    # Accumulate one vector per biological group, not a dense cell-by-gene table.
    common = np.ones(len(panel), dtype=bool)
    detected_by_group, sizes = {}, {}
    family = None
    seen = set()
    for p in training_profiles:
        if p.partition != 'train' or p.uid in seen:
            raise ValueError('Unique training profiles required for anchor selection')
        if family is None:
            family = p.unit
        if p.unit != family:
            raise ValueError('Select anchors within one measurement family')
        seen.add(p.uid)
        p = align(p, panel)
        common &= p.measured
        if p.group_id not in detected_by_group:
            detected_by_group[p.group_id] = np.zeros(len(panel), dtype=np.int64)
            sizes[p.group_id] = 0
        detected_by_group[p.group_id] += (p.values > 0)
        sizes[p.group_id] += 1
    if not seen:
        raise ValueError('Anchor selection requires training profiles')
    detection = sum(v / sizes[k] for k,v in detected_by_group.items()) / len(sizes)
    candidates = [i for i in range(len(panel)) if common[i] and detection[i] > 0]
    candidates.sort(key=lambda i: (-detection[i], panel[i]))
    if len(candidates) < count:
        raise ValueError('Insufficient training-measured detectable anchors')
    # Preserve panel order so the mask and frozen normalizer have identical axes.
    chosen = set(candidates[:count])
    return tuple(panel[i] for i in range(len(panel)) if i in chosen)


def make_mask_plan(profile, anchor_ids, seed, reward_fraction=0.2, evaluation_fraction=0.2):
    """Deterministic random masks; the fixed anchor identities come from training."""
    if not isinstance(seed, int) or isinstance(seed, bool) or seed < 0:
        raise ValueError('A nonnegative integer seed is required')
    if not 0 < reward_fraction < 1 or not 0 < evaluation_fraction < 1 or reward_fraction + evaluation_fraction >= 1:
        raise ValueError('Invalid mask fractions')
    anchor_ids = tuple(anchor_ids)
    idx = pd.Index(profile.genes).get_indexer(anchor_ids)
    if len(idx) == 0 or len(set(anchor_ids)) != len(idx) or (idx < 0).any() or not profile.measured[idx].all():
        raise ValueError('Fixed anchors must be unique and measured')
    anchors = np.zeros(len(profile.genes), dtype=bool)
    anchors[idx] = True
    available = np.flatnonzero(profile.measured & ~anchors)
    n_reward = max(1, int(len(available)*reward_fraction))
    n_eval = max(1, int(len(available)*evaluation_fraction))
    if n_reward + n_eval > len(available):
        raise ValueError('Too few measured non-anchor genes for reward and evaluation')
    order = np.random.default_rng(seed).permutation(available)
    reward = np.zeros(len(profile.genes), dtype=bool); reward[order[:n_reward]] = True
    evaluation = np.zeros(len(profile.genes), dtype=bool); evaluation[order[n_reward:n_reward+n_eval]] = True
    context = profile.measured & ~reward & ~evaluation
    plan = MaskPlan(context, reward, evaluation, anchors)
    plan.validate(profile.measured)
    return plan


@dataclass(frozen=True)
class Observation:
    """Only context values are exposed. Zeros elsewhere require the masks."""
    genes: tuple[str, ...]
    context_values: np.ndarray
    measured: np.ndarray
    context: np.ndarray
    reward: np.ndarray
    evaluation: np.ndarray
    anchors: np.ndarray


class RewardService:
    """Feedback service owns reward values, and never contains evaluation values."""
    def __init__(self, values, max_queries):
        if not isinstance(max_queries, int) or isinstance(max_queries, bool) or max_queries < 1:
            raise ValueError('Reward query budget must be a positive integer')
        self.__reward = values.copy()
        self.__locked = False
        self.__max_queries = max_queries
        self.__queries_used = 0

    @property
    def max_queries(self):
        return self.__max_queries

    @property
    def queries_used(self):
        return self.__queries_used

    def reward_loss(self, predictions):
        if self.__locked:
            raise ValueError('Adaptation has ended')
        if self.__queries_used >= self.__max_queries:
            raise ValueError('Reward query budget exhausted')
        p = np.asarray(predictions, dtype=float)
        if p.shape != self.__reward.shape or not np.isfinite(p).all():
            raise ValueError('Supply one finite prediction per reward gene')
        self.__queries_used += 1
        return float(np.mean((p - self.__reward)**2))

    def stop(self):
        """End adaptation and return legitimately observed reward-gene values."""
        if self.__locked:
            raise ValueError('Adaptation already stopped')
        self.__locked = True
        return self.__reward.copy()

    @property
    def stopped(self):
        return self.__locked


class EvaluationService:
    """A separate evaluator must be kept outside the adaptation loop."""
    def __init__(self, values, reward_service):
        self.__evaluation = values.copy()
        self.__reward_service = reward_service
        self.__locked = False

    def evaluate_once(self, predictions):
        if self.__locked:
            raise ValueError('Evaluation already performed')
        if not self.__reward_service.stopped:
            raise ValueError('Stop adaptation before final evaluation')
        p = np.asarray(predictions, dtype=float)
        if p.shape != self.__evaluation.shape or not np.isfinite(p).all():
            raise ValueError('Supply one finite prediction per held-out gene')
        self.__locked = True
        return float(np.mean((p - self.__evaluation)**2))


@dataclass(frozen=True)
class VisibleAnchorTotalNormalizer:
    """Sparse-compatible candidate: size factor is sum of fixed visible anchors."""
    anchors: tuple[str, ...]
    scale: float = 10000.0

    @classmethod
    def fit(cls, training, anchors, groups, scale=10000.0):
        anchors = tuple(anchors)
        if not anchors or len(set(anchors)) != len(anchors) or not training.index.is_unique or not training.columns.is_unique:
            raise ValueError('Need unique training axes and anchor IDs')
        a = training.loc[list(anchors)].to_numpy(float)
        if not np.isfinite(a).all() or (a < 0).any() or not (a.sum(axis=0) > 0).all():
            raise ValueError('Every training profile must have a positive measured anchor total')
        if groups.reindex(training.columns).isna().any() or not np.isfinite(scale) or scale <= 0:
            raise ValueError('Training groups and a positive scale are required')
        return cls(anchors, float(scale))

    def factor(self, observed):
        a = observed.loc[list(self.anchors)].to_numpy(float)
        if not np.isfinite(a).all() or (a < 0).any() or not a.sum() > 0:
            raise ValueError('Visible anchor total must be finite and positive')
        return float(a.sum())

    def transform(self, values, factor):
        v = values.to_numpy(float)
        if not np.isfinite(v).all() or (v < 0).any() or not np.isfinite(factor) or factor <= 0:
            raise ValueError('Invalid expression or visible-anchor factor')
        return pd.Series(np.log1p(v * (self.scale / factor)) / np.log(2), index=values.index)


def make_episode(profile, plan, normalizer, *, max_reward_queries=64):
    """Scale using context anchors only and expose hidden values via loss APIs."""
    plan.validate(profile.measured)
    anchors = tuple(np.asarray(profile.genes)[plan.anchors])
    if len(anchors) != len(normalizer.anchors) or set(anchors) != set(normalizer.anchors):
        raise ValueError('Plan anchors differ from frozen training reference')
    visible = pd.Series(profile.values[plan.context], index=np.asarray(profile.genes)[plan.context])
    factor = normalizer.factor(visible)
    target = np.zeros(len(profile.genes))
    active = plan.context | plan.reward | plan.evaluation
    target[active] = normalizer.transform(pd.Series(profile.values[active], index=np.asarray(profile.genes)[active]), factor).to_numpy()
    context = np.zeros(len(profile.genes))
    context[plan.context] = target[plan.context]
    obs = Observation(profile.genes, context, profile.measured.copy(), plan.context.copy(),
                      plan.reward.copy(), plan.evaluation.copy(), plan.anchors.copy())
    reward_service = RewardService(target[plan.reward], max_reward_queries)
    return obs, reward_service, EvaluationService(target[plan.evaluation], reward_service)


def fit_normalizer(training_profiles, anchors, normalizer_class):
    """The only fitting entry point accepts one unit family and train profiles."""
    a = tuple(anchors)
    if not a or len(set(a)) != len(a):
        raise ValueError('Unique anchors are required')
    columns = {}
    groups = {}
    family = None
    for p in training_profiles:
        if p.partition != 'train' or p.uid in columns:
            raise ValueError('Fit on unique training profiles only')
        if family is None:
            family = p.unit
        if p.unit != family:
            raise ValueError('Fit one declared measurement family at a time')
        idx = pd.Index(p.genes).get_indexer(a)
        if (idx < 0).any() or not p.measured[idx].all():
            raise ValueError('Training anchors must all be measured')
        columns[p.uid] = p.values[idx]
        groups[p.uid] = p.group_id
    if not columns:
        raise ValueError('Fit normalization on training profiles only')
    return normalizer_class.fit(pd.DataFrame(columns, index=a), a, pd.Series(groups))
