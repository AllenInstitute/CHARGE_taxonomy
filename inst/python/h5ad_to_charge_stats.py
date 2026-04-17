#!/usr/bin/env python3
"""
h5ad_to_charge_stats.py (PATCHED FOR R-PARITY)

Key guarantees vs R-only path:
- Normalization logic matches chargeTaxonomy():
    * use adata.X if already log-normalized
    * otherwise compute log2(CPM + 1) from raw counts
- All numeric ops in float64 (no float32 drift)
- Aggregation order unchanged
- Cell ordering preserved
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp


# -------------------------
# Argument parsing
# -------------------------

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--subsample", type=int, default=100000000)
    p.add_argument("--weight_by", choices=["cell", "cluster"], default="cell")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--layer", default=None)
    p.add_argument("--chunk_size", type=int, default=5000)
    p.add_argument("--max_pca", type=int, default=50)
    p.add_argument("--pca_sample_cells", type=int, default=50000)
    return p.parse_args()


# -------------------------
# Utilities
# -------------------------

def to_csr(X):
    if sp.issparse(X):
        return X.tocsr()
    return sp.csr_matrix(np.asarray(X))


def get_raw_counts(adata, layer=None):
    if layer is not None:
        X = adata.layers[layer]
        src = f"layers['{layer}']"
    elif getattr(adata, "raw", None) is not None and getattr(adata.raw, "X", None) is not None:
        X = adata.raw.X
        src = "raw.X"
    else:
        X = adata.X
        src = "X"
    print(f"[counts] using {src}")
    return to_csr(X)


def log2cpm_by_row(X_csr, chunk_size=5000):
    """
    EXACT match to R log2CPM_byRow():
      log2( (counts / rowSum) * 1e6 + 1 )
    """
    X_csr = to_csr(X_csr)
    n_cells, n_genes = X_csr.shape

    out = np.zeros((n_cells, n_genes), dtype=np.float64)

    row_sums = np.asarray(X_csr.sum(axis=1)).reshape(-1).astype(np.float64)
    row_sums[row_sums == 0] = 1.0
    scale = 1e6 / row_sums

    for start in range(0, n_cells, chunk_size):
        end = min(start + chunk_size, n_cells)
        block = X_csr[start:end].toarray().astype(np.float64)
        block *= scale[start:end, None]
        out[start:end, :] = np.log2(block + 1.0)

    return out


def get_norm_matrix(adata, raw_csr, keep_idx, chunk_size):
    """
    R-parity normalization rule:
      if adata.X exists AND looks log-normalized -> use it
      else -> compute log2(CPM + 1)
    """
    X = getattr(adata, "X", None)

    if X is not None:
        X_sub = X[keep_idx, :]
        if sp.issparse(X_sub):
            x_max = X_sub.max()
            x_max = float(x_max.A[0, 0]) if hasattr(x_max, "A") else float(x_max)
            if np.isfinite(x_max) and x_max <= 100:
                print("[norm] using adata.X (already normalized)")
                return X_sub.toarray().astype(np.float64)
        else:
            X_sub = np.asarray(X_sub, dtype=np.float64)
            x_max = np.nanmax(X_sub)
            if np.isfinite(x_max) and x_max <= 100:
                print("[norm] using adata.X (already normalized)")
                return X_sub

    print("[norm] computing log2(CPM + 1) from raw counts")
    return log2cpm_by_row(raw_csr, chunk_size=chunk_size)


# -------------------------
# Main
# -------------------------

def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    import anndata as ad
    adata = ad.read_h5ad(args.h5ad)

    hierarchy = adata.uns["hierarchy"]
    cluster_col = hierarchy[0]

    clusters = np.asarray(adata.obs[cluster_col])
    n_cells = len(clusters)

    rng = np.random.default_rng(args.seed)
    if args.subsample >= n_cells:
        keep_idx = np.arange(n_cells)
    else:
        keep_idx = rng.choice(n_cells, size=args.subsample, replace=False)
        keep_idx = np.sort(keep_idx)

    # -------------------------
    # Load counts
    # -------------------------

    raw_csr = get_raw_counts(adata, layer=args.layer)
    raw_csr = raw_csr[keep_idx, :]

    # -------------------------
    # Normalization (R-parity)
    # -------------------------

    norm_dense = get_norm_matrix(
        adata,
        raw_csr,
        keep_idx,
        chunk_size=args.chunk_size
    )

    # -------------------------
    # Cluster-level stats
    # -------------------------

    cluster_labels = clusters[keep_idx]
    uniq_clusters, inv = np.unique(cluster_labels, return_inverse=True)
    n_clusters = len(uniq_clusters)
    n_genes = raw_csr.shape[1]

    sums = np.zeros((n_genes, n_clusters), dtype=np.float64)
    counts = np.zeros((n_genes, n_clusters), dtype=np.float64)
    means = np.zeros((n_genes, n_clusters), dtype=np.float64)
    sds = np.zeros((n_genes, n_clusters), dtype=np.float64)
    count_n = np.zeros(n_clusters, dtype=np.int64)

    for k in range(n_clusters):
        idx = np.where(inv == k)[0]
        count_n[k] = len(idx)

        if len(idx) == 0:
            continue

        rc = raw_csr[idx]
        nc = norm_dense[idx]

        sums[:, k] = rc.sum(axis=0).A1 if sp.issparse(rc) else rc.sum(axis=0)
        counts[:, k] = (rc > 0).sum(axis=0).A1 if sp.issparse(rc) else (rc > 0).sum(axis=0)
        means[:, k] = nc.mean(axis=0)
        sds[:, k] = nc.std(axis=0, ddof=0)

    props = counts / count_n

    # -------------------------
    # Write outputs
    # -------------------------

    gene_names = adata.var_names
    cl_names = uniq_clusters.astype(str)

    def write_mat(name, mat):
        df = pd.DataFrame(mat, index=gene_names, columns=cl_names)
        df.insert(0, "gene", df.index)
        df.to_csv(outdir / f"{name}.csv", index=False)

    write_mat("sums", sums)
    write_mat("counts", counts)
    write_mat("props", props)
    write_mat("means", means)
    write_mat("sds", sds)

    pd.DataFrame({
        "group": cl_names,
        "n_cells": count_n
    }).to_csv(outdir / "count_n.csv", index=False)

    print("[done] wrote all precomputed outputs")


if __name__ == "__main__":
    main()