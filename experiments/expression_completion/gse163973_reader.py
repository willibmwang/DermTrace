"""Read the deposited GSE163973 processed 10x matrices one sample at a time.

No archive member is extracted to disk. Source feature IDs and symbols remain
separate; model alignment and per-cell eligibility belong to the caller.
"""
from __future__ import annotations

from contextlib import nullcontext
import gzip
from pathlib import Path
import re
import tarfile

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import issparse

SAMPLE_ALIASES = {'KL1': 'KF1', 'KL2': 'KF2', 'KL3': 'KF3',
                  'NS1': 'NF1', 'NS2': 'NF2', 'NS3': 'NF3'}


def iter_matrices(archive, samples=None, *, strict_members=False, thread_limit=None):
    """Yield sample, original feature table, barcode index and sparse counts.

    ``samples`` selects source sample IDs. Fail on missing or duplicated sample
    matrices, missing/duplicated nested components, invalid axes or counts.
    Input hash verification is the caller's responsibility.
    """
    selected = set(SAMPLE_ALIASES.values()) if samples is None else set(samples)
    if not selected or not selected <= set(SAMPLE_ALIASES.values()):
        raise ValueError('Unknown or empty source sample selection')
    with tarfile.open(Path(archive), 'r:') as outer:
        members = {}
        for member in outer.getmembers():
            if not member.isfile():
                continue
            match = re.search(r'_(KL[123]|NS[123])_matrix\.tar\.gz$', member.name)
            if not match:
                if strict_members:
                    raise ValueError('Unexpected file in count archive: ' + member.name)
                continue
            sample = SAMPLE_ALIASES[match.group(1)]
            if sample in members:
                raise ValueError('Duplicate sample matrix in archive')
            members[sample] = member
        if not selected <= set(members):
            raise ValueError('Requested source sample missing from archive')
        for sample, member in members.items():
            if sample not in selected:
                continue
            nested = outer.extractfile(member)
            if nested is None:
                raise ValueError('Unable to read selected sample matrix')
            with tarfile.open(fileobj=nested, mode='r:gz') as inner:
                components = {}
                for component in inner.getmembers():
                    if not component.isfile():
                        continue
                    name = Path(component.name).name
                    if name in components:
                        raise ValueError('Duplicate component in sample matrix')
                    components[name] = component
                needed = {'features.tsv.gz', 'barcodes.tsv.gz', 'matrix.mtx.gz'}
                if not needed <= components.keys():
                    raise ValueError('Missing component in sample matrix')
                with gzip.open(inner.extractfile(components['features.tsv.gz']), 'rt') as f:
                    features = pd.read_csv(f, sep='\t', header=None, dtype=str)
                with gzip.open(inner.extractfile(components['barcodes.tsv.gz']), 'rt') as f:
                    barcodes = pd.Index([line.strip() for line in f if line.strip()])
                if thread_limit is None:
                    thread_context = nullcontext()
                else:
                    from threadpoolctl import threadpool_limits
                    thread_context = threadpool_limits(limits=thread_limit)
                with gzip.open(inner.extractfile(components['matrix.mtx.gz']), 'rb') as f, thread_context:
                    raw_matrix = mmread(f)
                if not issparse(raw_matrix):
                    raise ValueError('Expected a sparse coordinate count matrix')
                if (not np.isfinite(raw_matrix.data).all() or (raw_matrix.data < 0).any()
                        or (raw_matrix.data != np.floor(raw_matrix.data)).any()):
                    raise ValueError('Invalid cell counts')
                coo = raw_matrix.tocoo(copy=False)
                if coo.shape[0] * coo.shape[1] > np.iinfo(np.int64).max:
                    raise ValueError('Count matrix dimensions exceed the coordinate check')
                coordinates = coo.row.astype(np.int64) * coo.shape[1] + coo.col
                coordinates.sort()
                if np.any(coordinates[1:] == coordinates[:-1]):
                    raise ValueError('Repeated count coordinates in sample matrix')
                del coordinates
                matrix = coo.tocsc()
                del coo, raw_matrix
            if (features.shape[1] < 2 or features.empty or features.iloc[:, 0].isna().any()
                    or features.iloc[:, 0].duplicated().any() or not barcodes.is_unique
                    or len(barcodes) == 0 or matrix.shape != (len(features), len(barcodes))):
                raise ValueError('Cell matrix feature or barcode axes are invalid')
            # Explicit stored zeros are not detected features.
            matrix.eliminate_zeros()
            yield sample, features, barcodes, matrix
            # Release the prior sparse sample before reading the next archive.
            del matrix, features, barcodes
