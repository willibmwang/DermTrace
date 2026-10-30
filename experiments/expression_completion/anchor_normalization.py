"""Candidate completion normalizers. Validated mechanics, not a chosen harmonizer.

Inputs are gene-indexed linear nonnegative measurements from one declared
measurement family. References are fitted on source training groups only.
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class FrozenAnchorNormalizer:
    anchors: tuple[str, ...]
    reference: tuple[float, ...]
    scale: float = 10000.0

    @classmethod
    def fit(cls, training, anchors, groups, scale=10000.0):
        anchors = tuple(anchors)
        if not anchors or len(set(anchors)) != len(anchors):
            raise ValueError("Supply a nonempty unique visible anchor panel.")
        if not training.index.is_unique or not training.columns.is_unique:
            raise ValueError("Training gene and sample axes must be unique.")
        a = training.loc[list(anchors)].to_numpy(float)
        if not np.isfinite(a).all() or not (a > 0).all():
            raise ValueError("This candidate requires anchors positive in every training sample.")
        groups = groups.reindex(training.columns)
        if groups.isna().any():
            raise ValueError("Every source training sample needs a biological group.")
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError("Use a positive fixed output scale.")
        weights = groups.map(1 / groups.value_counts()).to_numpy(float)
        weights /= weights.sum()
        totals = a.sum(axis=0)
        if not np.isfinite(totals).all() or not (totals > 0).all():
            raise ValueError("Training anchor totals must be finite and positive.")
        relative = a / totals
        reference = np.exp(np.sum(np.log(relative) * weights[None, :], axis=1))
        reference /= reference.sum()
        if not np.isfinite(reference).all() or not (reference > 0).all():
            raise ValueError("Training anchors cannot define a positive finite reference.")
        return cls(anchors, tuple(reference.tolist()), float(scale))

    def factor(self, observed):
        if not observed.index.is_unique:
            raise ValueError("Observed gene IDs must be unique.")
        a = observed.loc[list(self.anchors)].to_numpy(float)
        if not np.isfinite(a).all() or (a < 0).any():
            raise ValueError("Observed anchors must be finite nonnegative linear values.")
        factor = float(np.median(a / np.asarray(self.reference)))
        if not np.isfinite(factor) or factor <= 0:
            raise ValueError("Visible anchor panel cannot define a valid positive factor.")
        return factor

    def transform(self, values, factor):
        v = values.to_numpy(float)
        if not np.isfinite(v).all() or (v < 0).any():
            raise ValueError("This target requires finite nonnegative linear measurements.")
        if not np.isfinite(factor) or factor <= 0:
            raise ValueError("Use the fixed visible-anchor factor.")
        with np.errstate(divide="ignore"):
            transformed = np.logaddexp(0, np.log(v) + np.log(self.scale) - np.log(factor)) / np.log(2)
        if not np.isfinite(transformed).all():
            raise ValueError("Normalized output must remain finite.")
        return pd.Series(transformed, index=values.index)


def visible_total_factor(observed, anchors):
    a = observed.loc[list(anchors)].to_numpy(float)
    with np.errstate(over="ignore"):
        total = float(a.sum())
    if not np.isfinite(a).all() or (a < 0).any() or not np.isfinite(total) or total <= 0:
        raise ValueError("Visible anchor panel cannot define a valid positive total.")
    return total


def percentile_relative_to_anchors(values, observed_anchors):
    a = np.sort(observed_anchors.to_numpy(float))
    v = values.to_numpy(float)
    if len(a) == 0 or not np.isfinite(a).all() or not np.isfinite(v).all():
        raise ValueError("Relative assay-signal ranks need finite values and nonempty anchors.")
    # Signed values are permitted here because logged assays can be negative.
    low = np.searchsorted(a, v, side="left")
    high = np.searchsorted(a, v, side="right")
    return pd.Series((low + high) / (2 * len(a)), index=values.index)
