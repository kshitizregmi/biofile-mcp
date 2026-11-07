"""Single-cell (AnnData .h5ad) tools."""
import warnings

import anndata as ad
import numpy as np
import scipy.sparse as sp


def h5ad_summary(path: str, column: str | None = None) -> dict:
    """Summarize a single-cell .h5ad file without loading the full matrix:
    number of cells and genes, metadata columns, embeddings, neighbor graphs,
    and whether X looks like raw counts. If column is given (e.g. 'louvain' or
    'cell_type'), also return the number of cells in each category."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        adata = ad.read_h5ad(path, backed="r")

    sample = adata.X[:100]
    if sp.issparse(sample):
        values = sample.data
    else:
        values = np.asarray(sample).ravel()
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
            adata.file.close()
            raise ValueError(f"Column '{column}' not found. Available: {list(adata.obs.columns)}")
        counts = adata.obs[column].value_counts()
        summary["cells_per_" + column] = {str(k): int(v) for k, v in counts.items()}

    adata.file.close()
    return summary
