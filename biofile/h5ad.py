"""Single-cell AnnData (h5ad) summaries."""
import warnings

import anndata as ad
import numpy as np
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
