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


def collect_samples(skin_rows, cardiac_rows, melanoma_rows):
    samples = []
    skin = sorted({required(r, 'source_donor_key') for r in skin_rows})
    if skin != ['KF1', 'KF2', 'KF3', 'NF1', 'NF2', 'NF3']:
        raise ValueError('Expected all six GSE163973 source sample labels')
    for sample in skin:
        samples.append(dict(source_id=SKIN, sample_id=sample,
                            group_id=f'{SKIN}:pair_block:{sample[-1]}',
                            primary_eligible=True, quality_flag='cell_qc_required'))
    for row in cardiac_rows:
        sample, donor = required(row, 'sample_id'), required(row, 'donor_id')
        quality = required(row, 'quality_flag')
        eligible = boolean(row, 'published_primary_eligible')
        if eligible and quality != 'no_specific_paper_flag':
            raise ValueError('Flagged cardiac sample cannot be primary eligible')
        samples.append(dict(source_id=CARDIAC, sample_id=sample,
                            group_id=f'{CARDIAC}:donor:{donor}',
                            primary_eligible=eligible, quality_flag=quality))
    for row in melanoma_rows:
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
    if len(groups[SKIN]) != 3 or len(groups[CARDIAC]) < 3 or len(groups[MELANOMA]) != 2:
        raise ValueError('Require three skin blocks, at least three cardiac donors, and two melanoma lineages')
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
        for i, validation in enumerate(groups[src]):
            roles = {g: 'validation' if g == validation else 'train' for g in groups[src]}
            roles.update({g: 'test' for g in groups[target]})
            specs.append((f'{name}_sourceval_{i+1:02d}', name, roles, None))
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
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if any((output / f'{name}.json').exists() for name in splits) or (output / 'manifest.json').exists():
        raise ValueError('Refusing to overwrite existing split records; choose a new output directory')
    files = []
    for name, data in splits.items():
        path = output / f'{name}.json'
        path.write_text(json.dumps(data, indent=2) + '\n')
        files.append({'file': path.name, 'sha256': sha256(path), 'scenario': data['scenario']})
    manifest = {'schema_version': 1, 'files': files,
                'input_metadata_sha256': metadata_hashes,
                'total_files': len(files)}
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('skin', 'cardiac', 'melanoma', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    paths = {name: getattr(args, name) for name in ('skin', 'cardiac', 'melanoma')}
    samples = collect_samples(*(read_rows(paths[name]) for name in ('skin', 'cardiac', 'melanoma')))
    manifest = write_splits(create_splits(samples), paths, args.output)
    print(f'Wrote {manifest["total_files"]} exploratory split files.')


if __name__ == '__main__':
    main()
