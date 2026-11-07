import anndata as ad
import numpy as np
import pytest

from biofile import h5ad_summary


def test_shape_and_metadata(data_dir):
    r = h5ad_summary(str(data_dir / "small.h5ad"))
    assert (r["n_cells"], r["n_genes"]) == (200, 100)
    assert r["cell_columns"] == ["cell_type"]
    assert r["embeddings"] == ["X_umap"]


def test_detects_raw_counts(data_dir):
    r = h5ad_summary(str(data_dir / "small.h5ad"))
    assert r["X_looks_like"] == "raw counts"


def test_detects_normalized_data(tmp_path):
    rng = np.random.default_rng(0)
    adata = ad.AnnData(X=rng.normal(size=(50, 20)).astype(np.float32))  # scaled values, negatives
    path = tmp_path / "scaled.h5ad"
    adata.write_h5ad(path)

    r = h5ad_summary(str(path))
    assert r["X_looks_like"] == "normalized or scaled"


def test_cells_per_category(data_dir):
    r = h5ad_summary(str(data_dir / "small.h5ad"), column="cell_type")
    assert r["cells_per_cell_type"] == {"T cell": 100, "B cell": 60, "NK cell": 40}


def test_unknown_column_lists_available_ones(data_dir):
    with pytest.raises(ValueError, match="Available: \\['cell_type'\\]"):
        h5ad_summary(str(data_dir / "small.h5ad"), column="nope")
