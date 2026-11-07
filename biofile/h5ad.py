"""Single-cell AnnData (h5ad) summaries."""
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

# Number of cells read from X when deciding whether it holds raw counts
_SAMPLE_CELLS = 100


def h5ad_summary(path: str, column: str | None = None) -> dict:
    """Summarize an AnnData ``.h5ad`` file without loading the expression matrix.

    The file is opened in backed (read-only) mode, so only metadata and a small
    sample of ``X`` are read from disk. This keeps the call fast and memory use
    low for files of any size.

    Parameters
    ----------
    path : str
        Path to an ``.h5ad`` file.
    column : str, optional
        Name of a column in ``obs`` (for example ``"cell_type"``, ``"leiden"``
        or ``"sample"``). When given, the number of cells in each category is
        returned as well.

    Returns
    -------
    dict
        ``file``
            The path that was read.
        ``n_cells``, ``n_genes``
            Shape of the expression matrix.
        ``cell_columns``, ``gene_columns``
            Column names in ``obs`` and ``var``.
        ``embeddings``
            Keys in ``obsm``, e.g. ``X_pca`` or ``X_umap``.
        ``neighbor_graphs``
            Keys in ``obsp``, e.g. ``connectivities``.
        ``X_looks_like``
            ``"raw counts"`` if the sampled values are non-negative integers,
            otherwise ``"normalized or scaled"``. Methods such as DESeq2 or
            scVI expect raw counts.
        ``cells_per_<column>``
            Present only when ``column`` is given: cell count per category,
            largest first.

    Raises
    ------
    ValueError
        If ``column`` is not a column of ``obs``. The message lists the
        available columns.
    OSError
        If the file is missing or is not a valid HDF5 file.
    """
    # Files written by older anndata versions trigger format-migration
    # warnings on load; they don't affect the result.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        adata = ad.read_h5ad(path, backed="r")

    try:
        sample = adata.X[:_SAMPLE_CELLS]
        values = sample.data if sp.issparse(sample) else np.asarray(sample).ravel()
        raw_counts = bool(np.all(values == np.round(values)) and values.min() >= 0)

        summary = {
            "file": path,
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "cell_columns": list(adata.obs.columns),
            "gene_columns": list(adata.var.columns),
            "embeddings": list(adata.obsm.keys()),
            "neighbor_graphs": list(adata.obsp.keys()),
            "X_looks_like": "raw counts" if raw_counts else "normalized or scaled",
        }

        if column:
            if column not in adata.obs.columns:
                raise ValueError(
                    f"Column '{column}' not found. Available: {list(adata.obs.columns)}"
                )
            counts = adata.obs[column].value_counts()
            summary[f"cells_per_{column}"] = {str(k): int(v) for k, v in counts.items()}

        return summary
    finally:
        adata.file.close()


def h5ad_column_stats(path: str, column: str, bins: int = 40) -> dict:
    """Distribution of one cell-metadata (``obs``) column, such as a QC metric.

    For numeric columns (for example ``n_counts``, ``n_genes`` or
    ``percent_mito``) returns summary statistics and a histogram; for
    categorical or text columns returns cell counts per category.

    Parameters
    ----------
    path : str
        Path to an ``.h5ad`` file.
    column : str
        Name of a column in ``obs``.
    bins : int, default 40
        Number of histogram bins for numeric columns.

    Returns
    -------
    dict
        ``file``, ``column``, ``kind`` (``"numeric"`` or ``"categorical"``)
        and ``n_cells``. Numeric columns add ``min``, ``q1``, ``median``,
        ``mean``, ``q3``, ``max``, ``missing`` and a ``histogram`` with
        ``start``, ``end`` and ``count`` lists. Categorical columns add
        ``counts``: cells per category, largest first.

    Raises
    ------
    ValueError
        If ``column`` is not a column of ``obs``.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        adata = ad.read_h5ad(path, backed="r")
    try:
        if column not in adata.obs.columns:
            raise ValueError(f"Column '{column}' not found. Available: {list(adata.obs.columns)}")
        values = adata.obs[column]
        result = {"file": path, "column": column, "n_cells": int(adata.n_obs)}

        if pd.api.types.is_numeric_dtype(values) and not pd.api.types.is_bool_dtype(values):
            numbers = values.to_numpy(dtype=float)
            present = numbers[~np.isnan(numbers)]
            if present.size == 0:
                return {**result, "kind": "numeric", "missing": int(numbers.size)}
            counts, edges = np.histogram(present, bins=bins)
            q1, median, q3 = np.percentile(present, [25, 50, 75])
            return {
                **result,
                "kind": "numeric",
                "min": round(float(present.min()), 4),
                "q1": round(float(q1), 4),
                "median": round(float(median), 4),
                "mean": round(float(present.mean()), 4),
                "q3": round(float(q3), 4),
                "max": round(float(present.max()), 4),
                "missing": int(numbers.size - present.size),
                "histogram": {
                    "start": np.round(edges[:-1], 4).tolist(),
                    "end": np.round(edges[1:], 4).tolist(),
                    "count": counts.tolist(),
                },
            }

        counts = values.astype(str).value_counts()
        return {**result, "kind": "categorical",
                "counts": {str(k): int(v) for k, v in counts.items()}}
    finally:
        adata.file.close()


def h5ad_embedding(path: str, basis: str | None = None, color: str | None = None,
                   max_points: int = 20_000, seed: int = 0) -> dict:
    """Coordinates of a 2-D embedding (UMAP, t-SNE, PCA) for plotting.

    Intended for drawing, not for an AI client: the result can hold tens of
    thousands of points. Large datasets are randomly subsampled.

    Parameters
    ----------
    path : str
        Path to an ``.h5ad`` file.
    basis : str, optional
        Key in ``obsm``, e.g. ``"X_umap"``. Defaults to UMAP, then t-SNE,
        then the first available embedding.
    color : str, optional
        ``obs`` column whose values are returned alongside the coordinates.
    max_points : int, default 20000
        Most cells returned; larger datasets are sampled without replacement.
    seed : int, default 0
        Random seed for the sample, so the same cells are drawn each time.

    Returns
    -------
    dict
        ``file``, ``basis``, ``n_cells``, ``sampled``, ``x``, ``y`` and, if
        ``color`` is given, ``color`` (the column name) and ``values``.

    Raises
    ------
    ValueError
        If the file has no 2-D embedding, ``basis`` is not in ``obsm``, or
        ``color`` is not a column of ``obs``.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        adata = ad.read_h5ad(path, backed="r")
    try:
        keys = list(adata.obsm.keys())
        if basis is None:
            preferred = [k for k in ("X_umap", "X_tsne") if k in keys]
            basis = (preferred or keys or [None])[0]
        if basis is None:
            raise ValueError("This file has no embeddings in obsm (e.g. X_umap).")
        if basis not in keys:
            raise ValueError(f"Embedding '{basis}' not found. Available: {keys}")
        if color is not None and color not in adata.obs.columns:
            raise ValueError(f"Column '{color}' not found. Available: {list(adata.obs.columns)}")

        coords = np.asarray(adata.obsm[basis])[:, :2]
        n = coords.shape[0]
        keep = np.arange(n)
        if n > max_points:
            keep = np.sort(np.random.default_rng(seed).choice(n, max_points, replace=False))

        result = {
            "file": path,
            "basis": basis,
            "n_cells": int(n),
            "sampled": bool(n > max_points),
            "x": np.round(coords[keep, 0], 4).tolist(),
            "y": np.round(coords[keep, 1], 4).tolist(),
        }
        if color is not None:
            values = adata.obs[color].iloc[keep]
            numeric = pd.api.types.is_numeric_dtype(values) and not pd.api.types.is_bool_dtype(values)
            result["color"] = color
            result["values"] = (np.round(values.to_numpy(dtype=float), 4).tolist() if numeric
                                else values.astype(str).tolist())
        return result
    finally:
        adata.file.close()
