"""
scripts/stage8_heatmaps.py
────────────────────────────────────────────────────────────────────────────
STAGE 8 — DOMAIN-FUNCTION / DOMAIN-BINDING HEATMAPS
Builds association matrices and publication-style clustered heatmaps.
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
import sys as _sys
import os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from scripts.utils import (get_input_path, load_config, load_table,
                            make_output_dir, rename_to_canonical,
                            save_tsv, setup_logging, split_cell, stage_banner)


def _build_matrix(df: pd.DataFrame,
                  row_col: str, col_col: str,
                  row_delim: str = ";", col_delim: str = ";",
                  max_rows: int = 40, max_cols: int = 30) -> pd.DataFrame:
    """Build a co-occurrence count matrix from two multi-value columns."""
    counts: dict[str, dict[str, int]] = {}
    for _, row in df.iterrows():
        row_vals = split_cell(row.get(row_col, ""), row_delim)
        col_vals = split_cell(row.get(col_col, ""), col_delim)
        for rv in row_vals:
            for cv in col_vals:
                counts.setdefault(rv, {}).setdefault(cv, 0)
                counts[rv][cv] += 1

    mat = pd.DataFrame(counts).T.fillna(0).astype(int)
    # Trim to top entries by total count
    if len(mat) > max_rows:
        mat = mat.loc[mat.sum(axis=1).nlargest(max_rows).index]
    if len(mat.columns) > max_cols:
        mat = mat[mat.sum(axis=0).nlargest(max_cols).index]
    return mat


def _cluster_matrix(mat: pd.DataFrame,
                    method: str = "ward") -> pd.DataFrame:
    """Reorder matrix rows by hierarchical clustering."""
    if len(mat) < 2:
        return mat
    try:
        dist = pdist(mat.values, metric="euclidean")
        link = linkage(dist, method=method)
        order = dendrogram(link, no_plot=True)["leaves"]
        return mat.iloc[order]
    except Exception:
        return mat


def _heatmap(mat: pd.DataFrame, title: str, out_path: Path,
             cmap: str = "YlOrRd", dpi: int = 150,
             font_size: int = 8) -> None:
    """Render a publication-style heatmap using matplotlib only."""
    if mat.empty:
        return

    row_h = max(4, len(mat) * 0.25)
    col_w = max(5, len(mat.columns) * 0.3)
    fig, ax = plt.subplots(figsize=(col_w, row_h), dpi=dpi)

    im = ax.imshow(mat.values, aspect="auto", cmap=cmap,
                   interpolation="nearest")
    plt.colorbar(im, ax=ax, shrink=0.6, label="Count")

    ax.set_xticks(range(len(mat.columns)))
    ax.set_xticklabels(mat.columns, rotation=45, ha="right",
                       fontsize=font_size)
    ax.set_yticks(range(len(mat)))
    ax.set_yticklabels(mat.index, fontsize=font_size)
    ax.set_title(title, fontsize=font_size + 2, fontweight="bold", pad=10)

    # Annotate cells if matrix is small enough
    if len(mat) * len(mat.columns) <= 400:
        vmax = mat.values.max()
        for i in range(len(mat)):
            for j in range(len(mat.columns)):
                val = mat.iloc[i, j]
                if val > 0:
                    color = "white" if val > vmax * 0.6 else "black"
                    ax.text(j, i, str(val), ha="center", va="center",
                            fontsize=max(5, font_size - 2), color=color)

    plt.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _diversity_plots(occ: pd.DataFrame, domain_col: str, function_col: str,
                     binding_col: str, out_dir: Path, dpi: int,
                     fs: int, cmap: str) -> None:
    """Extra summary charts that make the sparse heatmaps more convincing."""
    from collections import Counter

    # ── 1. Binding diversity per domain (bar chart) ───────────────────────
    # Collect globally unique binding targets per domain type.
    # Using a set per domain ensures the same target is counted only once
    # even when it appears in multiple occurrences of the same domain type.
    _div_targets: dict[str, set] = {}
    for _, row in occ.iterrows():
        dt  = str(row.get(domain_col, "")).strip()
        bps = [b.strip() for b in str(row.get(binding_col, "")).split(";") if b.strip()]
        if dt and bps:
            if dt not in _div_targets:
                _div_targets[dt] = set()
            _div_targets[dt].update(bps)
    diversity: dict[str, int] = {dt: len(targets) for dt, targets in _div_targets.items()}

    if diversity:
        top_div = sorted(diversity.items(), key=lambda x: -x[1])[:30]
        labels = [x.replace("INTERPRO:", "").replace("PFAM:", "").replace("CDD:", "") for x, _ in top_div]
        vals   = [v for _, v in top_div]
        fig, ax = plt.subplots(figsize=(12, max(5, len(labels) * 0.35)), dpi=dpi)
        bars = ax.barh(labels[::-1], vals[::-1], color="#4878D0", edgecolor="white", linewidth=0.5)
        ax.set_xlabel("Number of distinct binding targets", fontsize=fs + 1)
        ax.set_title("Domain Binding Diversity\n(distinct binding targets per domain type)",
                     fontsize=fs + 3, fontweight="bold")
        ax.bar_label(bars, padding=3, fontsize=max(6, fs - 1))
        ax.spines[["top", "right"]].set_visible(False)
        plt.tight_layout()
        fig.savefig(out_dir / "domain_binding_diversity.png", dpi=dpi, bbox_inches="tight")
        plt.close(fig)

    # ── 2. Function tag frequency (bar chart) ─────────────────────────────
    func_counts: Counter = Counter()
    for _, row in occ.iterrows():
        for f in str(row.get(function_col, "")).split(";"):
            f = f.strip()
            if f:
                func_counts[f] += 1

    if func_counts:
        top_f = func_counts.most_common(30)
        labels_f = [x for x, _ in top_f]
        vals_f   = [v for _, v in top_f]
        fig, ax = plt.subplots(figsize=(12, max(5, len(labels_f) * 0.35)), dpi=dpi)
        bars_f = ax.barh(labels_f[::-1], vals_f[::-1], color="#6ACC65", edgecolor="white", linewidth=0.5)
        ax.set_xlabel("Number of domain occurrences", fontsize=fs + 1)
        ax.set_title("Functional Tag Frequency\n(occurrences per function tag across all domain types)",
                     fontsize=fs + 3, fontweight="bold")
        ax.bar_label(bars_f, padding=3, fontsize=max(6, fs - 1))
        ax.spines[["top", "right"]].set_visible(False)
        plt.tight_layout()
        fig.savefig(out_dir / "function_tag_frequency.png", dpi=dpi, bbox_inches="tight")
        plt.close(fig)

    # ── 3. Binding target coverage (how many domains link to each target) ──
    target_domain_count: Counter = Counter()
    for _, row in occ.iterrows():
        dt  = str(row.get(domain_col, "")).strip()
        bps = [b.strip() for b in str(row.get(binding_col, "")).split(";") if b.strip()]
        for bp in set(bps):
            target_domain_count[bp] += 1

    if target_domain_count:
        top_t = target_domain_count.most_common(30)
        labels_t = [x for x, _ in top_t]
        vals_t   = [v for _, v in top_t]
        fig, ax = plt.subplots(figsize=(12, max(5, len(labels_t) * 0.35)), dpi=dpi)
        bars_t = ax.barh(labels_t[::-1], vals_t[::-1], color="#D65F5F", edgecolor="white", linewidth=0.5)
        ax.set_xlabel("Number of domain types associating with this target", fontsize=fs + 1)
        ax.set_title("Binding Target Connectivity\n(domain types per binding target)",
                     fontsize=fs + 3, fontweight="bold")
        ax.bar_label(bars_t, padding=3, fontsize=max(6, fs - 1))
        ax.spines[["top", "right"]].set_visible(False)
        plt.tight_layout()
        fig.savefig(out_dir / "binding_target_connectivity.png", dpi=dpi, bbox_inches="tight")
        plt.close(fig)


def run(cfg: dict) -> None:
    logger = setup_logging(cfg.get("log_dir","logs"), "stage8_heatmaps")
    stage_banner(logger, "STAGE 8", "Domain-Function / Domain-Binding Heatmaps")
    out_dir = make_output_dir(cfg, "stage8_heatmaps")

    plot_cfg = cfg.get("plotting", {})
    dpi     = plot_cfg.get("figure_dpi", 150)
    cmap    = plot_cfg.get("colormap_heatmap", "YlOrRd")
    fs      = plot_cfg.get("font_size", 8)
    max_r   = plot_cfg.get("max_heatmap_rows", 40)
    max_c   = plot_cfg.get("max_heatmap_cols", 30)
    cluster_method = cfg.get("clustering", {}).get("linkage_method", "ward")

    occ_path = get_input_path(cfg, "domain_occurrence_table")
    if not occ_path:
        logger.warning("domain_occurrence_table not found — skipping Stage 8")
        return

    occ = load_table(occ_path)
    aliases = cfg.get("column_aliases", {})
    occ = rename_to_canonical(occ, aliases)

    domain_col   = "domain_type_id"
    function_col = "function_tags"
    # After rename_to_canonical the typo 'binding_parnter' is corrected to
    # the canonical name 'binding_partner' — check both to be safe.
    binding_col = next(
        (c for c in ["binding_partner", "binding_parnter"] if c in occ.columns),
        "binding_partner",
    )

    for col in [domain_col, function_col, binding_col]:
        if col not in occ.columns:
            logger.warning(f"Column '{col}' not found — some heatmaps skipped")

    # Domain × Function matrix
    if domain_col in occ.columns and function_col in occ.columns:
        logger.info("Building domain × function matrix …")
        mat_func = _build_matrix(occ, domain_col, function_col,
                                 max_rows=max_r, max_cols=max_c)
        mat_func = _cluster_matrix(mat_func, cluster_method)
        save_tsv(mat_func.reset_index().rename(columns={"index": domain_col}),
                 out_dir / "domain_function_matrix.tsv", logger)
        _heatmap(mat_func, "Domain Type × Function Tag", 
                 out_dir / "domain_function_heatmap.png", cmap, dpi, fs)
        logger.info(f"Domain×function matrix: {mat_func.shape}")

    # Domain × Binding matrix
    if domain_col in occ.columns and binding_col in occ.columns:
        logger.info("Building domain × binding matrix …")
        mat_bind = _build_matrix(occ, domain_col, binding_col,
                                 max_rows=max_r, max_cols=max_c)
        mat_bind = _cluster_matrix(mat_bind, cluster_method)
        save_tsv(mat_bind.reset_index().rename(columns={"index": domain_col}),
                 out_dir / "domain_binding_matrix.tsv", logger)
        _heatmap(mat_bind, "Domain Type × Binding Target",
                 out_dir / "domain_binding_heatmap.png", cmap, dpi, fs)
        logger.info(f"Domain×binding matrix: {mat_bind.shape}")

    # ── Supplementary plots: diversity summaries ──────────────────────────
    if domain_col in occ.columns and function_col in occ.columns and binding_col in occ.columns:
        _diversity_plots(occ, domain_col, function_col, binding_col,
                         out_dir, dpi, fs, cmap)

    logger.info("Stage 8 complete.")


if __name__ == "__main__":
    import sys, os
    # Always run from the ontology_pipeline project root so that
    # relative paths (input/, output/, config/) resolve correctly
    # regardless of where PyCharm / the terminal launched the script from.
    _project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(_project_root)
    _cfg_path = sys.argv[1] if len(sys.argv) > 1 else "config/pipeline_config.yaml"
    run(load_config(_cfg_path))
