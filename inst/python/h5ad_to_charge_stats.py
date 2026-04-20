
#!/usr/bin/env python3
"""
h5ad_to_charge_stats.py

Hybrid precompute script for CHARGE_taxonomy.

This version supports exact subsample parity with the R-only path by accepting
a CSV of cell IDs chosen by R (via subsampleCells()).

Required for exact parity:
  --subsample_ids stats.dir/subsample_ids.csv

The CSV must contain one column:
  cell_id

If --subsample_ids is provided, Python uses those exact cells and ignores its
own internal subsampling logic.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp


def parse_args():
    p = argparse.ArgumentParser(description="Precompute CHARGE stats from an h5ad file")
    p.add_argument("--h5ad", required=True, help="Path to input .h5ad file (local)")
    p.add_argument("--outdir", default="./charge_stats", help="Output directory")
    p.add_argument("--subsample", type=int, default=100000000, help="Max cells to keep")
    p.add_argument("--subsample_ids", default=None,
                   help="Optional CSV of exact cell IDs chosen by R; if provided, overrides internal subsampling")
    p.add_argument("--weight_by", default="cell", choices=["cell", "cluster"],
                   help="Aggregate higher levels by cells or by cluster")
    p.add_argument("--seed", type=int, default=42, help="Random seed")
    p.add_argument("--cluster_col", default=None,
                   help="Optional override for highest-resolution cluster column")
    p.add_argument("--layer", default=None,
                   help="Optional layer name for raw counts (else raw.X then X)")
    p.add_argument("--chunk_size", type=int, default=5000, help="Chunk size for normalization")
    p.add_argument("--max_pca", type=int, default=50, help="Maximum PCs to retain after thresholding")
    p.add_argument("--pca_th", type=float, default=0.5,
                   help="Z-score threshold on explained-variance ratio")
    p.add_argument("--pca_sample_cells", type=int, default=50000,
                   help="How many cells to use to estimate PCA loadings")
    return p.parse_args()


def load_anndata(path):
    import anndata as ad
    print(f"[load] reading {path}")
    adata = ad.read_h5ad(path)
    print(f"[load] shape = {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    return adata


def get_hierarchy(adata):
    hier = adata.uns.get("hierarchy", None)
    if hier is None:
        raise ValueError("uns['hierarchy'] is required")
    if isinstance(hier, dict):
        keys = list(hier.keys())
        vals = np.array([float(hier[k]) for k in keys])
        order = np.argsort(-vals)
        return [keys[i] for i in order]
    if isinstance(hier, (list, tuple, np.ndarray)):
        return list(hier)
    raise TypeError(f"Unsupported hierarchy type: {type(hier)}")


def to_csr(X):
    if sp.issparse(X):
        return X.tocsr()
    return sp.csr_matrix(np.asarray(X))


def get_raw_counts(adata, layer=None):
    if layer is not None:
        if layer not in adata.layers:
            raise KeyError(f"Layer {layer!r} not found in adata.layers")
        X = adata.layers[layer]
        src = f"layers[{layer!r}]"
    elif getattr(adata, "raw", None) is not None and getattr(adata.raw, "X", None) is not None:
        X = adata.raw.X
        src = "raw.X"
    else:
        X = adata.X
        src = "X"
    print(f"[counts] using {src}")
    return to_csr(X)


def subsample_cells(cluster_labels, subsample, seed=42):
    """
    Fallback only. This is NOT guaranteed to match R exactly.
    Exact matching requires using --subsample_ids generated in R.
    """
    n = len(cluster_labels)
    if subsample is None or subsample >= n:
        return np.arange(n)
    rng = np.random.default_rng(seed)
    cluster_labels = np.asarray(cluster_labels)
    clusters, _ = np.unique(cluster_labels, return_counts=True)
    per_cluster_target = max(1, subsample // len(clusters))
    picked = []
    for cl in clusters:
        idx = np.where(cluster_labels == cl)[0]
        k = min(len(idx), per_cluster_target)
        if k > 0:
            picked.append(rng.choice(idx, size=k, replace=False))
    picked = np.sort(np.concatenate(picked)) if picked else np.array([], dtype=int)
    if len(picked) < subsample:
        remaining = np.setdiff1d(np.arange(n), picked, assume_unique=False)
        extra = min(len(remaining), subsample - len(picked))
        if extra > 0:
            picked = np.sort(np.concatenate([picked, rng.choice(remaining, size=extra, replace=False)]))
    elif len(picked) > subsample:
        picked = np.sort(rng.choice(picked, size=subsample, replace=False))
    return picked


def get_keep_idx(adata, cluster_col, subsample, seed=42, subsample_ids=None):
    """
    If subsample_ids is provided, use it exactly.
    Otherwise use Python fallback subsampling.
    """
    if subsample_ids is not None:
        ids_df = pd.read_csv(subsample_ids)
        if "cell_id" not in ids_df.columns:
            raise ValueError(f"{subsample_ids} must contain a column named 'cell_id'")
        wanted = ids_df["cell_id"].astype(str).tolist()

        obs_names = pd.Index(adata.obs_names.astype(str))
        missing = [x for x in wanted if x not in obs_names]
        if len(missing) > 0:
            raise ValueError(
                f"{len(missing)} cell_id values in {subsample_ids} were not found in the h5ad obs_names. "
                f"First few missing: {missing[:5]}"
            )

        keep_idx = obs_names.get_indexer(wanted)
        print(f"[subsample] using exact R-generated subsample_ids ({len(keep_idx)} cells)")
        return keep_idx

    cluster_labels_full = np.asarray(adata.obs[cluster_col]).astype(str)
    keep_idx = subsample_cells(cluster_labels_full, subsample, seed=seed)
    print(f"[subsample] using Python fallback subsampling ({len(keep_idx)} cells)")
    return keep_idx


def log2cpm_by_row(X_csr, chunk_size=5000):
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


def get_norm_matrix(adata, raw_csr_sub, keep_idx, chunk_size=5000):
    """
    Match R-only normalization logic:
      - use adata.X directly if it appears already normalized (max <= 100)
      - otherwise compute log2(CPM + 1) from raw counts
    """
    Xn = getattr(adata, "X", None)
    if Xn is not None:
        Xn_sub = Xn[keep_idx, :]
        if sp.issparse(Xn_sub):
            x_max = Xn_sub.max()
            if hasattr(x_max, "A"):
                x_max = float(x_max.A[0, 0])
            else:
                x_max = float(x_max)
            if np.isfinite(x_max) and x_max <= 100:
                print("[norm] using adata.X (already log-normalized)")
                return Xn_sub.toarray().astype(np.float64)
        else:
            Xn_sub = np.asarray(Xn_sub, dtype=np.float64)
            x_max = np.nanmax(Xn_sub)
            if np.isfinite(x_max) and x_max <= 100:
                print("[norm] using adata.X (already log-normalized)")
                return Xn_sub

    print("[norm] computing log2(CPM + 1) from raw counts")
    return log2cpm_by_row(raw_csr_sub, chunk_size=chunk_size)


def beta_score_from_props(props_mat):
    p = np.asarray(props_mat, dtype=np.float64)
    return np.mean(4.0 * p * (1.0 - p), axis=1)


def select_variable_genes(props_df, all_gene_names, n_max=1200):
    props_arr = props_df.to_numpy(dtype=float)
    max_prop = props_arr.max(axis=1)
    mask = max_prop > 0.5
    if mask.sum() == 0:
        mask = np.ones_like(max_prop, dtype=bool)
    scores = beta_score_from_props(props_arr[mask, :])
    genes = np.asarray(all_gene_names)[mask]
    order = np.argsort(scores)
    return list(genes[order[: min(n_max, len(order))]])


def compute_group_stats(raw_csr, norm_dense, groups, ordered_groups):
    """
    Exact intended CHARGE statistics:
      sums   = raw count sums
      counts = number of cells with raw count > 0
      props  = counts / n_cells
      means  = mean(normalized expression)
      sds    = sqrt(E[x^2] - mean^2)   [population variance]
    """
    raw_csr = to_csr(raw_csr)
    groups = np.asarray(groups)
    ordered_groups = list(ordered_groups)

    n_cells, n_genes = raw_csr.shape
    n_groups = len(ordered_groups)

    sums = np.zeros((n_genes, n_groups), dtype=np.float64)
    counts = np.zeros((n_genes, n_groups), dtype=np.float64)
    means = np.zeros((n_genes, n_groups), dtype=np.float64)
    sds = np.zeros((n_genes, n_groups), dtype=np.float64)
    count_n = np.zeros(n_groups, dtype=np.int64)

    for gi, g in enumerate(ordered_groups):
        idx = np.where(groups == g)[0]
        count_n[gi] = len(idx)
        if len(idx) == 0:
            continue

        sub_raw = raw_csr[idx, :]
        sub_norm = norm_dense[idx, :]

        sums[:, gi] = np.asarray(sub_raw.sum(axis=0)).reshape(-1)
        counts[:, gi] = np.asarray((sub_raw > 0).sum(axis=0)).reshape(-1)
        means[:, gi] = sub_norm.mean(axis=0)

        sqr_means = np.mean(np.square(sub_norm, dtype=np.float64), axis=0)
        vars_ = sqr_means - np.square(means[:, gi])
        vars_[vars_ < 0] = 0.0
        sds[:, gi] = np.sqrt(vars_)

    props = counts / np.maximum(count_n.astype(np.float64), 1.0)[None, :]

    return {
        "sums": sums,
        "counts": counts,
        "props": props,
        "means": means,
        "sds": sds,
        "count_n": dict(zip(ordered_groups, [int(x) for x in count_n])),
    }


def run_rd_pca_no_sklearn(norm_dense, gene_names, sample_names, variable_genes,
                          max_pca=50, th=0.5, seed=42, sample_cells=50000):
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    use_idx = [gene_to_idx[g] for g in variable_genes if g in gene_to_idx]
    if len(use_idx) < 2:
        raise ValueError("Too few variable genes for PCA")

    X_all = norm_dense[:, use_idx].astype(np.float64)
    n_cells = X_all.shape[0]
    rng = np.random.default_rng(seed)

    if n_cells > sample_cells:
        samp_idx = np.sort(rng.choice(np.arange(n_cells), size=sample_cells, replace=False))
    else:
        samp_idx = np.arange(n_cells)

    X_s = X_all[samp_idx, :]
    mu = X_s.mean(axis=0, keepdims=True)
    X_s = X_s - mu

    U, S, Vt = np.linalg.svd(X_s, full_matrices=False)

    eig = S ** 2
    total = eig.sum()
    v = eig / total if total > 0 else np.ones_like(eig) / len(eig)

    if np.std(v) == 0:
        selected = np.arange(min(len(v), max_pca))
    else:
        z = (v - np.mean(v)) / np.std(v)
        selected = np.where(z > th)[0]
        selected = selected[:max_pca]
        if len(selected) == 0:
            raise RuntimeError("No PCs selected by rd_PCA threshold logic")

    V = Vt.T[:, selected]
    Xc_all = X_all - mu
    scores = Xc_all @ V
    cols = [f"PC{i+1}" for i in selected]
    return pd.DataFrame(scores, index=sample_names, columns=cols), selected


def extract_umap(adata, keep_idx, rd_df, seed=42):
    default_embedding = None
    if "default_embedding" in adata.uns:
        de = adata.uns["default_embedding"]
        if isinstance(de, (list, tuple, np.ndarray)) and len(de) > 0:
            default_embedding = str(de[0])
        elif isinstance(de, str):
            default_embedding = de

    if default_embedding is not None and default_embedding in adata.obsm:
        emb = np.asarray(adata.obsm[default_embedding])[keep_idx, :2]
        return pd.DataFrame(emb, index=rd_df.index, columns=["UMAP_1", "UMAP_2"]), default_embedding, "obsm"

    if len(adata.obsm.keys()) > 0:
        key = list(adata.obsm.keys())[0]
        emb = np.asarray(adata.obsm[key])[keep_idx, :2]
        return pd.DataFrame(emb, index=rd_df.index, columns=["UMAP_1", "UMAP_2"]), key, "obsm"

    try:
        import umap
    except Exception:
        print("[umap] umap-learn not available; skipping UMAP computation", file=sys.stderr)
        return None, None, "missing"

    reducer = umap.UMAP(n_components=2, random_state=seed)
    emb = reducer.fit_transform(rd_df.to_numpy())
    return pd.DataFrame(emb, index=rd_df.index, columns=["UMAP_1", "UMAP_2"]), None, "computed"


def save_matrix(path, mat, row_names, col_names):
    df = pd.DataFrame(mat, index=row_names, columns=col_names)
    df.index.name = "gene"
    df.to_csv(path)


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    adata = load_anndata(args.h5ad)
    hierarchy = get_hierarchy(adata)
    if args.cluster_col is not None:
        if args.cluster_col not in adata.obs.columns:
            raise KeyError(f"cluster_col {args.cluster_col!r} not found in obs")
        hierarchy = [args.cluster_col] + [h for h in hierarchy if h != args.cluster_col]

    cluster_col = hierarchy[0]

    # Exact R-generated subsample if provided
    keep_idx = get_keep_idx(
        adata,
        cluster_col=cluster_col,
        subsample=args.subsample,
        seed=args.seed,
        subsample_ids=args.subsample_ids,
    )

    obs_sub = adata.obs.iloc[keep_idx].copy()
    sample_names = list(obs_sub.index.astype(str))

    raw_csr = get_raw_counts(adata, layer=args.layer)[keep_idx, :]
    norm_dense = get_norm_matrix(adata, raw_csr, keep_idx, chunk_size=args.chunk_size)
    gene_names = list(adata.var_names.astype(str))

    cluster_labels = np.asarray(obs_sub[cluster_col]).astype(str)
    all_clusters = list(pd.unique(cluster_labels))
    print(f"[stats] cluster level: {len(all_clusters)} groups")
    st0 = compute_group_stats(raw_csr, norm_dense, cluster_labels, all_clusters)

    col_names = list(all_clusters)
    sums = st0["sums"]
    counts = st0["counts"]
    props = st0["props"]
    means = st0["means"]
    sds = st0["sds"]
    count_n = dict(st0["count_n"])

    if args.weight_by == "cell" and len(hierarchy) > 1:
        for level_col in hierarchy[1:]:
            level_labels = np.asarray(obs_sub[level_col]).astype(str)
            groups = list(pd.unique(level_labels))
            print(f"[stats] level {level_col}: {len(groups)} groups")
            st = compute_group_stats(raw_csr, norm_dense, level_labels, groups)
            sums = np.concatenate([sums, st["sums"]], axis=1)
            counts = np.concatenate([counts, st["counts"]], axis=1)
            props = np.concatenate([props, st["props"]], axis=1)
            means = np.concatenate([means, st["means"]], axis=1)
            sds = np.concatenate([sds, st["sds"]], axis=1)
            count_n.update(st["count_n"])
            col_names.extend(groups)

    save_matrix(outdir / "sums.csv", sums, gene_names, col_names)
    save_matrix(outdir / "counts.csv", counts, gene_names, col_names)
    save_matrix(outdir / "props.csv", props, gene_names, col_names)
    save_matrix(outdir / "means.csv", means, gene_names, col_names)
    save_matrix(outdir / "sds.csv", sds, gene_names, col_names)

    pd.DataFrame({
        "group": list(count_n.keys()),
        "n_cells": list(count_n.values())
    }).to_csv(outdir / "count_n.csv", index=False)

    obs_sub.to_csv(outdir / "obs_info.csv")
    adata.var.to_csv(outdir / "var_info.csv")

    print(f"[save] sums.csv: {sums.shape}")
    print(f"[save] counts.csv: {counts.shape}")
    print(f"[save] props.csv: {props.shape}")
    print(f"[save] means.csv: {means.shape}")
    print(f"[save] sds.csv: {sds.shape}")

    if "highly_variable_genes_standard" in adata.var.columns:
        mask = np.asarray(adata.var["highly_variable_genes_standard"]).astype(bool)
        variable_genes = [g for g, m in zip(gene_names, mask) if m]
        if len(variable_genes) == 0:
            variable_genes = select_variable_genes(
                pd.DataFrame(props[:, :len(all_clusters)], index=gene_names, columns=all_clusters),
                gene_names
            )
    else:
        variable_genes = select_variable_genes(
            pd.DataFrame(props[:, :len(all_clusters)], index=gene_names, columns=all_clusters),
            gene_names
        )

    variable_genes = [g for g in variable_genes if g in set(gene_names)]
    if len(variable_genes) < 50:
        raise RuntimeError("<50 valid variable genes after selection")

    rd_df, selected = run_rd_pca_no_sklearn(
        norm_dense, gene_names, sample_names, variable_genes,
        max_pca=args.max_pca, th=args.pca_th,
        seed=args.seed, sample_cells=args.pca_sample_cells
    )
    rd_df.to_csv(outdir / "rd_dat.csv")
    print(f"[save] rd_dat.csv: {rd_df.shape}")

    umap_df, embedding_key, embedding_source = extract_umap(adata, keep_idx, rd_df, seed=args.seed)
    if umap_df is not None:
        umap_df.to_csv(outdir / "umap.csv")
        print("[save] umap.csv")

    uns_data = {
        "hierarchy": hierarchy,
        "cluster_col": cluster_col,
        "default_embedding": embedding_key,
        "embedding_source": embedding_source,
        "weight_by": args.weight_by,
        "selected_pcs": [int(x) for x in selected],
        "variable_genes": variable_genes,
        "group_columns": col_names,
    }
    with open(outdir / "uns_data.json", "w") as f:
        json.dump(uns_data, f, indent=2)

    manifest = {
        "h5ad": args.h5ad,
        "outdir": str(outdir),
        "subsample": args.subsample,
        "subsample_ids": args.subsample_ids,
        "weight_by": args.weight_by,
        "seed": args.seed,
        "cluster_col": args.cluster_col,
        "layer": args.layer,
        "chunk_size": args.chunk_size,
        "pca_sample_cells": args.pca_sample_cells,
        "n_cells_kept": int(len(keep_idx)),
        "n_genes": int(len(gene_names)),
        "n_groups_written": int(len(col_names)),
    }
    with open(outdir / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)

    print("[done] wrote all precomputed outputs")


if __name__ == "__main__":
    main()
