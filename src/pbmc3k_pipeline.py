#!/usr/bin/env python
"""Reusable scverse best-practices scRNA-seq pipeline (10x PBMC 3k).

Steps: load 10x -> QC metrics -> MAD-based cell filtering -> normalise + log1p
-> HVG selection -> PCA -> kNN graph -> Leiden -> UMAP -> marker-based
cell-type annotation.  Every figure is written to ``--figdir``.

Usage
-----
    python src/pbmc3k_pipeline.py \
        --input data/filtered_gene_bc_matrices/hg19 \
        --figdir figures --output results/pbmc3k_annotated.h5ad

The same functions are imported by ``notebooks/pbmc3k_walkthrough.ipynb``.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation

SEED = 0

# Canonical PBMC markers. Order matters only for display.
MARKERS: dict[str, list[str]] = {
    "Naive CD4 T": ["IL7R", "CCR7", "LEF1", "TCF7"],
    "Memory CD4 T": ["IL7R", "S100A4", "CD2"],
    "CD8 T": ["CD8A", "CD8B", "GZMK", "CCL5"],
    "NK": ["GNLY", "NKG7", "GZMB", "PRF1"],
    "B": ["MS4A1", "CD79A", "CD79B", "CD74"],
    "CD14+ Mono": ["CD14", "LYZ", "S100A8", "S100A9"],
    "FCGR3A+ Mono": ["FCGR3A", "MS4A7", "COTL1", "LST1"],
    "Dendritic": ["FCER1A", "CLEC10A", "HLA-DPB1"],
    "Megakaryocyte": ["PPBP", "PF4", "GNG11"],
}


# --------------------------------------------------------------------------- #
# QC
# --------------------------------------------------------------------------- #
def is_outlier(adata: ad.AnnData, metric: str, nmads: float) -> pd.Series:
    """True where ``metric`` deviates more than ``nmads`` MADs from the median."""
    m = adata.obs[metric]
    med, mad = np.median(m), median_abs_deviation(m)
    return (m < med - nmads * mad) | (m > med + nmads * mad)


def add_qc_metrics(adata: ad.AnnData) -> None:
    """Annotate mito / ribosomal / haemoglobin genes and compute QC metrics."""
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains(r"^HB[^P]")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )


def flag_low_quality(adata: ad.AnnData, nmads: float = 5, mt_nmads: float = 5,
                     mt_floor: float = 8.0) -> None:
    """Add boolean ``outlier`` and ``mt_outlier`` columns to ``adata.obs``.

    Follows the sc-best-practices recipe: 5 MADs on log1p(total_counts),
    log1p(n_genes) and % counts in top-20 genes; the same 5 MADs on % mito
    (3 MADs would drop ~8% of PBMC 3k because the mito MAD is tiny), plus an
    absolute ceiling (``mt_floor`` %).
    """
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", nmads)
        | is_outlier(adata, "log1p_n_genes_by_counts", nmads)
        | is_outlier(adata, "pct_counts_in_top_20_genes", nmads)
    )
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", mt_nmads) | (
        adata.obs["pct_counts_mt"] > mt_floor
    )


# --------------------------------------------------------------------------- #
# Annotation
# --------------------------------------------------------------------------- #
def annotate_clusters(adata: ad.AnnData, groupby: str = "leiden",
                      markers: dict[str, list[str]] = MARKERS) -> pd.DataFrame:
    """Assign each cluster the cell type with the highest mean marker score.

    A per-cell score (``sc.tl.score_genes``) is computed for each marker set and
    averaged per cluster; the argmax is the label.  Returns the cluster x
    cell-type score table so the decision can be inspected.
    """
    score_cols = []
    for ct, genes in markers.items():
        genes = [g for g in genes if g in adata.var_names]
        col = f"score_{ct}"
        sc.tl.score_genes(adata, genes, score_name=col, random_state=SEED)
        score_cols.append(col)
    table = adata.obs.groupby(groupby, observed=True)[score_cols].mean()
    table.columns = [c.removeprefix("score_") for c in table.columns]
    return table


# --------------------------------------------------------------------------- #
# Pipeline
# --------------------------------------------------------------------------- #
def savefig(fig, figdir: Path, name: str) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figdir / f"{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def run(input_dir: Path, figdir: Path, output: Path | None, resolution: float = 1.0,
        n_pcs: int = 30, n_hvg: int = 2000) -> ad.AnnData:
    sc.settings.verbosity = 2
    sc.set_figure_params(dpi=100, frameon=False, figsize=(5, 5))

    # 1. Load -------------------------------------------------------------- #
    adata = sc.read_10x_mtx(input_dir, var_names="gene_symbols", cache=False)
    adata.var_names_make_unique()
    print(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")

    # 2. QC ---------------------------------------------------------------- #
    add_qc_metrics(adata)
    flag_low_quality(adata)
    n0 = adata.n_obs
    keep = ~(adata.obs["outlier"] | adata.obs["mt_outlier"])
    print(f"QC: removing {n0 - keep.sum()} / {n0} cells "
          f"({adata.obs['outlier'].sum()} MAD outliers, "
          f"{adata.obs['mt_outlier'].sum()} mito outliers)")

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    sc.pl.violin(adata, "n_genes_by_counts", groupby=None, ax=axes[0], show=False)
    sc.pl.violin(adata, "total_counts", ax=axes[1], show=False)
    sc.pl.violin(adata, "pct_counts_mt", ax=axes[2], show=False)
    savefig(fig, figdir, "01_qc_violins_prefilter")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt",
                  ax=axes[0], show=False)
    sc.pl.scatter(adata, "total_counts", "pct_counts_mt", color="outlier",
                  ax=axes[1], show=False)
    savefig(fig, figdir, "02_qc_scatter_prefilter")

    adata = adata[keep].copy()
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"After filtering: {adata.n_obs} cells x {adata.n_vars} genes")

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    sc.pl.violin(adata, "n_genes_by_counts", ax=axes[0], show=False)
    sc.pl.violin(adata, "total_counts", ax=axes[1], show=False)
    sc.pl.violin(adata, "pct_counts_mt", ax=axes[2], show=False)
    savefig(fig, figdir, "03_qc_violins_postfilter")

    # 3. Normalise --------------------------------------------------------- #
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # keep full log-normalised matrix for DE / plotting

    # 4. HVG --------------------------------------------------------------- #
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat")
    print(f"Highly variable genes: {int(adata.var['highly_variable'].sum())}")
    fig = plt.figure(figsize=(6, 4.5))
    sc.pl.highly_variable_genes(adata, show=False)
    savefig(plt.gcf(), figdir, "04_highly_variable_genes")

    # 5. PCA --------------------------------------------------------------- #
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=50, svd_solver="arpack", random_state=SEED)
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.uns["pca"] = adata_hvg.uns["pca"]
    adata.varm["PCs"] = np.zeros((adata.n_vars, 50))
    adata.varm["PCs"][adata.var["highly_variable"].to_numpy()] = adata_hvg.varm["PCs"]
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
    savefig(plt.gcf(), figdir, "05_pca_variance_ratio")

    # 6. Graph / cluster / embed ------------------------------------------ #
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs, random_state=SEED)
    sc.tl.leiden(adata, resolution=resolution, flavor="leidenalg", n_iterations=-1,
                 directed=False, random_state=SEED, key_added="leiden")
    sc.tl.umap(adata, random_state=SEED)
    print("Clusters:", adata.obs["leiden"].nunique())

    ax = sc.pl.umap(adata, color="leiden", legend_loc="on data", show=False,
                    title=f"Leiden (resolution={resolution})")
    savefig(ax.figure, figdir, "06_umap_leiden")

    axs = sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "pct_counts_mt"],
                     ncols=3, show=False)
    savefig(axs[0].figure, figdir, "07_umap_qc_covariates")

    # 7. Markers ----------------------------------------------------------- #
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=True)
    sc.pl.rank_genes_groups(adata, n_genes=8, sharey=False, show=False)
    savefig(plt.gcf(), figdir, "08_rank_genes_groups")
    top = pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(5)
    print("Top-5 DE genes per cluster:\n", top.T.to_string())

    # 8. Annotation -------------------------------------------------------- #
    table = annotate_clusters(adata)
    mapping = table.idxmax(axis=1).to_dict()
    adata.obs["cell_type"] = adata.obs["leiden"].map(mapping).astype("category")
    print("Cluster -> cell type:\n", pd.Series(mapping).to_string())

    fig, ax = plt.subplots(figsize=(8, 5))
    im = ax.imshow(table.values, aspect="auto", cmap="viridis")
    ax.set_xticks(range(table.shape[1]), table.columns, rotation=45, ha="right")
    ax.set_yticks(range(table.shape[0]), table.index)
    ax.set_ylabel("Leiden cluster")
    fig.colorbar(im, label="mean marker score")
    savefig(fig, figdir, "09_cluster_marker_scores")

    flat = list(dict.fromkeys(g for v in MARKERS.values() for g in v if g in adata.raw.var_names))
    ax = sc.pl.dotplot(adata, flat, groupby="leiden", standard_scale="var",
                       show=False, return_fig=True)
    ax.savefig(figdir / "10_marker_dotplot_by_cluster.png", dpi=150, bbox_inches="tight")
    plt.close("all")

    ax = sc.pl.umap(adata, color="cell_type", show=False, title="PBMC 3k - annotated")
    savefig(ax.figure, figdir, "11_umap_cell_types")

    canon = ["IL7R", "CCR7", "CD8A", "NKG7", "GNLY", "MS4A1", "CD14", "LYZ",
             "FCGR3A", "FCER1A", "PPBP"]
    axs = sc.pl.umap(adata, color=canon, ncols=4, show=False, use_raw=True)
    savefig(axs[0].figure, figdir, "12_umap_canonical_markers")

    dp = sc.pl.dotplot(adata, flat, groupby="cell_type", standard_scale="var",
                       show=False, return_fig=True)
    dp.savefig(figdir / "13_marker_dotplot_by_celltype.png", dpi=150, bbox_inches="tight")
    plt.close("all")

    counts = adata.obs["cell_type"].value_counts()
    fig, ax = plt.subplots(figsize=(7, 4))
    (counts / counts.sum() * 100).sort_values().plot.barh(ax=ax, color="steelblue")
    ax.set_xlabel("% of cells")
    savefig(fig, figdir, "14_cell_type_proportions")
    print("Cell-type counts:\n", counts.to_string())

    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(output, compression="gzip")
        print(f"Saved {output}")
    return adata


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", type=Path,
                   default=Path("data/filtered_gene_bc_matrices/hg19"))
    p.add_argument("--figdir", type=Path, default=Path("figures"))
    p.add_argument("--output", type=Path, default=Path("results/pbmc3k_annotated.h5ad"))
    p.add_argument("--resolution", type=float, default=1.0)
    p.add_argument("--n-pcs", type=int, default=30)
    p.add_argument("--n-hvg", type=int, default=2000)
    a = p.parse_args()
    run(a.input, a.figdir, a.output, a.resolution, a.n_pcs, a.n_hvg)


if __name__ == "__main__":
    main()
