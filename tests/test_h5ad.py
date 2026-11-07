"""Tests for :func:`biofile.h5ad_summary`."""
import anndata as ad
import numpy as np
import pandas as pd
import pytest

from biofile import h5ad_column_stats, h5ad_embedding, h5ad_summary


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
    adata = ad.AnnData(X=rng.normal(size=(50, 20)).astype(np.float32))  # scaled data includes negatives
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


def qc_h5ad(tmp_path):
    rng = np.random.default_rng(1)
    obs = pd.DataFrame({
        "n_counts": rng.integers(500, 5000, 300).astype(float),
        "percent_mito": rng.random(300) * 0.1,
        "leiden": pd.Categorical(rng.choice(["0", "1", "2"], 300)),
    }, index=[f"c{i}" for i in range(300)])
    adata = ad.AnnData(X=rng.poisson(1, (300, 20)).astype(np.float32), obs=obs)
    adata.obsm["X_umap"] = rng.normal(size=(300, 2))
    path = tmp_path / "qc.h5ad"
    adata.write_h5ad(path)
    return str(path), obs


def test_numeric_column_stats(tmp_path):
    path, obs = qc_h5ad(tmp_path)
    r = h5ad_column_stats(path, "n_counts")
    assert r["kind"] == "numeric" and r["missing"] == 0
    assert r["median"] == pytest.approx(obs["n_counts"].median())
    assert (r["min"], r["max"]) == (obs["n_counts"].min(), obs["n_counts"].max())
    assert sum(r["histogram"]["count"]) == 300 and len(r["histogram"]["count"]) == 40


def test_categorical_column_stats(tmp_path):
    path, obs = qc_h5ad(tmp_path)
    r = h5ad_column_stats(path, "leiden")
    assert r["kind"] == "categorical"
    assert r["counts"] == obs["leiden"].value_counts().to_dict()
    with pytest.raises(ValueError, match="Available"):
        h5ad_column_stats(path, "nope")


def test_embedding(tmp_path, data_dir):
    path, _ = qc_h5ad(tmp_path)
    e = h5ad_embedding(path, color="leiden")
    assert (e["basis"], e["n_cells"], e["sampled"]) == ("X_umap", 300, False)
    assert len(e["x"]) == len(e["y"]) == len(e["values"]) == 300
    sample = h5ad_embedding(path, color="n_counts", max_points=50)
    assert sample["sampled"] and len(sample["values"]) == 50
    assert isinstance(sample["values"][0], float)
    assert h5ad_embedding(str(data_dir / "small.h5ad"))["basis"] == "X_umap"
    with pytest.raises(ValueError, match="not found"):
        h5ad_embedding(path, basis="X_tsne")
