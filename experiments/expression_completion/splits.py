"""Make exploratory completion splits from supplied processed metadata.

The six skin labels are grouped by suffix as a precaution, not a verified
patient mapping. Never interpret mixed-skin folds as proven donor holdouts.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

SKIN = 'GSE163973'
CARDIAC = 'GSE144424'
MELANOMA = 'GSE95280_recount3'


def sha256(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def read_rows(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f))


def required(row, key):
    value = row.get(key)
    if not isinstance(value, str) or not value.strip() or value.strip().lower() in {'nan', 'none', 'unknown'}:
        raise ValueError(f'Missing {key}')
    return value.strip()


def boolean(row, key):
    value = required(row, key).lower()
    if value not in {'true', 'false'}:
        raise ValueError(f'{key} must be True or False')
    return value == 'true'


def collect_samples(skin_rows=None, cardiac_rows=None, melanoma_rows=None):
    samples = []
    if any(rows is not None and not rows for rows in (skin_rows, cardiac_rows, melanoma_rows)):
        raise ValueError('A supplied metadata source is empty')
    if skin_rows is not None:
        skin = sorted({required(r, 'source_donor_key') for r in skin_rows})
        if skin != ['KF1', 'KF2', 'KF3', 'NF1', 'NF2', 'NF3']:
            raise ValueError('Expected all six GSE163973 source sample labels')
        for sample in skin:
            samples.append(dict(source_id=SKIN, sample_id=sample,
                                group_id=f'{SKIN}:pair_block:{sample[-1]}',
                                primary_eligible=True, quality_flag='cell_qc_required'))
    for row in cardiac_rows or []:
        sample, donor = required(row, 'sample_id'), required(row, 'donor_id')
        quality = required(row, 'quality_flag')
        eligible = boolean(row, 'published_primary_eligible')
        if eligible and quality != 'no_specific_paper_flag':
            raise ValueError('Flagged cardiac sample cannot be primary eligible')
        samples.append(dict(source_id=CARDIAC, sample_id=sample,
                            group_id=f'{CARDIAC}:donor:{donor}',
                            primary_eligible=eligible, quality_flag=quality))
    for row in melanoma_rows or []:
        sample, family = required(row, 'sample_id'), required(row, 'cell_line_family')
        if family not in {'501MEL', 'IGR37_IGR39'}:
            raise ValueError('Unrecognized melanoma lineage; do not guess group')
        eligible = boolean(row, 'clean_hypoxia_eligible')
        samples.append(dict(source_id=MELANOMA, sample_id=sample,
                            group_id=f'{MELANOMA}:lineage:{family}',
                            primary_eligible=eligible,
                            quality_flag='clean_oxygen' if eligible else 'not_primary'))
    keys = [(r['source_id'], r['sample_id']) for r in samples]
    if len(keys) != len(set(keys)):
        raise ValueError('Repeated sample ID within source')
    return samples


def validate_assignments(rows):
    if not rows:
        raise ValueError('Empty assignment list')
    seen, group_roles = set(), {}
    for row in rows:
        key = (row['source_id'], row['sample_id'])
        group, role = row['group_id'], row['partition']
        if key in seen or not group or role not in {'train', 'validation', 'test'}:
            raise ValueError('Invalid or duplicate assignment')
        seen.add(key)
        if group in group_roles and group_roles[group] != role:
            raise ValueError('Group crosses partitions')
        group_roles[group] = role


def create_splits(samples):
    groups = {source: sorted({r['group_id'] for r in samples if r['source_id'] == source})
              for source in (SKIN, CARDIAC, MELANOMA)}
    if not any(groups.values()):
        raise ValueError('Supply at least one nonempty source')
    if groups[SKIN] and len(groups[SKIN]) != 3:
        raise ValueError('Skin requires three exploratory sample blocks')
    if groups[CARDIAC] and len(groups[CARDIAC]) < 3:
        raise ValueError('Cardiac requires at least three donor groups')
    if groups[MELANOMA] and len(groups[MELANOMA]) != 2:
        raise ValueError('This melanoma source requires its two known lineage groups')
    specs = []
    for source, name in ((SKIN, 'within_skin'), (CARDIAC, 'within_cardiac')):
        gs = groups[source]
        for i, test in enumerate(gs):
            val = gs[(i + 1) % len(gs)]
            roles = {g: 'test' if g == test else 'validation' if g == val else 'train' for g in gs}
            specs.append((f'{name}_holdout_{i+1:02d}', name, roles, None))
            if source == SKIN:
                specs.extend((f'within_{condition}_holdout_{i+1:02d}', f'within_{condition}', roles, prefix)
                             for condition, prefix in (('keloid', 'KF'), ('normal_scar', 'NF')))
    for src, target, name in ((SKIN, CARDIAC, 'skin_to_cardiac'),
                              (CARDIAC, SKIN, 'cardiac_to_skin'),
                              (SKIN, MELANOMA, 'skin_to_melanoma'),
                              (CARDIAC, MELANOMA, 'cardiac_to_melanoma'),
                              (MELANOMA, SKIN, 'melanoma_to_skin'),
                              (MELANOMA, CARDIAC, 'melanoma_to_cardiac')):
        if not groups[src] or not groups[target]:
            continue
        for i, validation in enumerate(groups[src]):
            roles = {g: 'validation' if g == validation else 'train' for g in groups[src]}
            roles.update({g: 'test' for g in groups[target]})
            specs.append((f'{name}_sourceval_{i+1:02d}', name, roles, None))
    if not specs:
        raise ValueError('Melanoma alone has only two lineages, insufficient for separate train, validation and test groups')
    result = {}
    for name, scenario, roles, prefix in specs:
        rows = [dict(row, partition=roles[row['group_id']]) for row in samples
                if row['group_id'] in roles and (prefix is None or row['sample_id'].startswith(prefix))]
        validate_assignments(rows)
        result[name] = {'schema_version': 1, 'name': name, 'scenario': scenario,
                        'status': 'exploratory; previously inspected cohorts',
                        'skin_group_caveat': 'Suffix blocks are sample groups, not proof of donor independence; unknown cross-condition pairings can cross roles.',
                        'assignments': rows}
    return result


def write_splits(splits, paths, output):
    metadata_hashes = {key: sha256(path) for key, path in paths.items()}
    source_keys = {SKIN: 'skin', CARDIAC: 'cardiac', MELANOMA: 'melanoma'}
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if any((output / f'{name}.json').exists() for name in splits) or (output / 'manifest.json').exists():
        raise ValueError('Refusing to overwrite existing split records; choose a new output directory')
    files = []
    for name, data in splits.items():
        source_ids = {row['source_id'] for row in data['assignments']}
        if not source_ids <= set(source_keys):
            raise ValueError('Unknown source in split assignment')
        needed = {source_keys[source] for source in source_ids}
        if not needed <= metadata_hashes.keys():
            raise ValueError('Missing metadata input for a source used by the split')
        input_hashes = {key: metadata_hashes[key] for key in sorted(needed)}
        data = {**data, 'input_metadata_sha256': input_hashes}
        path = output / f'{name}.json'
        path.write_text(json.dumps(data, indent=2) + '\n')
        files.append({'file': path.name, 'sha256': sha256(path), 'scenario': data['scenario'],
                      'input_metadata_sha256': input_hashes})
    manifest = {'schema_version': 1, 'files': files,
                'input_metadata_sha256': metadata_hashes,
                'total_files': len(files)}
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('skin', 'cardiac', 'melanoma'):
        parser.add_argument('--' + name, type=Path,
                            help='processed metadata for this source (optional)')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    paths = {name: getattr(args, name) for name in ('skin', 'cardiac', 'melanoma')
             if getattr(args, name) is not None}
    if not paths:
        parser.error('supply at least one processed metadata source')
    rows = {name: read_rows(path) for name, path in paths.items()}
    samples = collect_samples(rows.get('skin'), rows.get('cardiac'), rows.get('melanoma'))
    manifest = write_splits(create_splits(samples), paths, args.output)
    print(f'Wrote {manifest["total_files"]} exploratory split files.')


if __name__ == '__main__':
    main()
