"""
scripts/stage13_protein_modules.py
────────────────────────────────────────────────────────────────────────────
STAGE 13 — PROTEIN MODULE SEPARATION & CLUSTERING
Converts ontology-derived protein annotations into a feature matrix,
performs PCA and PCoA, hierarchical clustering, and benchmarks cluster
quality against known biological groups and MMseqs2 sequence clusters.
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations

import warnings
warnings.filterwarnings("ignore", message="Mean of empty slice")
warnings.filterwarnings("ignore", message="invalid value encountered in scalar divide")
warnings.filterwarnings("ignore", category=RuntimeWarning, module="numpy")
warnings.filterwarnings("ignore", category=FutureWarning, module="sklearn")
warnings.filterwarnings("ignore", category=UserWarning, message=".*unfilled marker.*")

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    silhouette_score,
    adjusted_rand_score,
    normalized_mutual_info_score,
)

import sys as _sys
import os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))

from scripts.utils import (
    get_input_path, load_config, load_table,
    make_output_dir, rename_to_canonical,
    save_tsv, setup_logging, split_cell, stage_banner
)


# ── Build binary feature matrix ───────────────────────────────────────────────

def build_feature_matrix(
    occ: pd.DataFrame,
    feature_cols: list[str],
    prot_col: str = "protein_accession",
    binary: bool = True,
) -> pd.DataFrame:
    """
    Returns a DataFrame: proteins × features (binary or count).
    Multi-value cells are split by semicolon.
    """
    counts: dict[str, dict[str, int]] = {}

    for _, row in occ.iterrows():
        prot = str(row.get(prot_col, "")).strip()
        if not prot or prot == "nan":
            continue
        counts.setdefault(prot, {})
        for col in feature_cols:
            for tok in split_cell(row.get(col, "")):
                feat_key = f"{col}::{tok}"
                counts[prot][feat_key] = counts[prot].get(feat_key, 0) + 1

    df = pd.DataFrame(counts).T.fillna(0)
    if df.empty:
        return df
    if binary:
        df = (df > 0).astype(int)
    return df


# ── Similarity matrix ─────────────────────────────────────────────────────────

def compute_similarity(feat_mat: pd.DataFrame, metric: str = "jaccard") -> pd.DataFrame:
    """
    Returns a protein × protein similarity DataFrame.
    Uses 1 - distance for distance metrics.
    """
    if feat_mat.empty:
        return pd.DataFrame()
    dist = squareform(pdist(feat_mat.values, metric=metric))
    sim = 1 - dist
    np.fill_diagonal(sim, 1.0)
    return pd.DataFrame(sim, index=feat_mat.index, columns=feat_mat.index)


# ── PCA ───────────────────────────────────────────────────────────────────────

def run_pca(feat_mat: pd.DataFrame, n_components: int = 10) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns (pca_coords DataFrame, explained_variance DataFrame).
    """
    if feat_mat.empty or len(feat_mat) < 2:
        return pd.DataFrame(), pd.DataFrame()

    n_comp = min(n_components, len(feat_mat) - 1, feat_mat.shape[1])
    pca = PCA(n_components=n_comp, random_state=42)
    scaled = StandardScaler().fit_transform(feat_mat.values.astype(float))
    coords = pca.fit_transform(scaled)

    coord_df = pd.DataFrame(
        coords,
        index=feat_mat.index,
        columns=[f"PC{i+1}" for i in range(n_comp)],
    )
    var_df = pd.DataFrame({
        "component": [f"PC{i+1}" for i in range(n_comp)],
        "explained_variance_ratio": pca.explained_variance_ratio_.round(4),
        "cumulative_variance": np.cumsum(pca.explained_variance_ratio_).round(4),
    })
    return coord_df, var_df


# ── True classical PCoA ──────────────────────────────────────────────────────

def run_pcoa(sim_mat: pd.DataFrame, n_components: int = 2) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Classical PCoA from a similarity matrix.
    Input: sim_mat = protein × protein similarity matrix (e.g. Jaccard similarity)

    Returns:
      coord_df = coordinates (PCoA1, PCoA2, ...)
      var_df   = explained variance / goodness-of-fit table
    """
    if sim_mat.empty or len(sim_mat) < 2:
        return pd.DataFrame(), pd.DataFrame()

    # Convert similarity to distance
    D = 1.0 - sim_mat.values.astype(float)
    D = np.clip(D, 0.0, 1.0)
    np.fill_diagonal(D, 0.0)

    n = D.shape[0]

    # Double-center squared distance matrix
    D2 = D ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J

    # Eigen decomposition
    eigvals, eigvecs = np.linalg.eigh(B)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    # Keep positive eigenvalues only
    pos = eigvals > 1e-10
    eigvals_pos = eigvals[pos]
    eigvecs_pos = eigvecs[:, pos]

    if len(eigvals_pos) == 0:
        return pd.DataFrame(), pd.DataFrame()

    n_comp = min(n_components, len(eigvals_pos))
    coords = eigvecs_pos[:, :n_comp] * np.sqrt(eigvals_pos[:n_comp])

    coord_df = pd.DataFrame(
        coords,
        index=sim_mat.index,
        columns=[f"PCoA{i+1}" for i in range(n_comp)],
    )

    explained = eigvals_pos / eigvals_pos.sum()
    var_df = pd.DataFrame({
        "component": [f"PCoA{i+1}" for i in range(len(explained))],
        "explained_variance_ratio": explained.round(4),
        "cumulative_variance": np.cumsum(explained).round(4),
    })

    coord_df.attrs["gof_2d"] = float(explained[:n_comp].sum())
    return coord_df, var_df


# ── Hierarchical clustering ───────────────────────────────────────────────────

def run_clustering(feat_mat: pd.DataFrame, n_clusters: int = 6, method: str = "ward") -> pd.Series:
    """
    Returns a Series: protein_accession → cluster_id (1-indexed).
    """
    if len(feat_mat) < 2:
        return pd.Series(dtype=int)

    # Ward linkage expects Euclidean distance
    dist = pdist(feat_mat.values, metric="euclidean")
    link = linkage(dist, method=method)
    labels = fcluster(link, t=n_clusters, criterion="maxclust")
    return pd.Series(labels, index=feat_mat.index, name="cluster_id")


# ── Group label assignment ────────────────────────────────────────────────────

def _assign_group(
    prot: str,
    group_map: dict[str, str],
    grp_df: pd.DataFrame | None,
    prot_col: str | None,
    grp_col: str | None,
) -> str:
    """Look up known biological group from protein_group_labels.tsv or fallback map."""
    if grp_df is not None and prot_col and grp_col and prot_col in grp_df.columns and grp_col in grp_df.columns:
        rows = grp_df[grp_df[prot_col].astype(str) == str(prot)]
        if not rows.empty:
            return str(rows.iloc[0][grp_col])

    for group, keywords in group_map.items():
        for kw in keywords:
            if kw.lower() in prot.lower():
                return group
    return "other"


# ── Within vs between group similarity ───────────────────────────────────────

def _group_similarity_stats(sim_mat: pd.DataFrame, group_series: pd.Series) -> pd.DataFrame:
    """Compute mean within-group and between-group similarity."""
    common = sim_mat.index.intersection(group_series.index)
    if common.empty:
        return pd.DataFrame()

    sim = sim_mat.loc[common, common]
    grp = group_series.loc[common]
    rows = []

    for g in grp.unique():
        members = grp[grp == g].index.tolist()
        others = grp[grp != g].index.tolist()
        if len(members) < 2:
            continue

        within_vals = [sim.loc[a, b] for i, a in enumerate(members) for b in members[i+1:]]
        between_vals = [sim.loc[a, b] for a in members for b in others]

        rows.append({
            "group": g,
            "n_proteins": len(members),
            "within_mean_sim": round(np.mean(within_vals), 4) if within_vals else 0,
            "between_mean_sim": round(np.mean(between_vals), 4) if between_vals else 0,
            "separation_ratio": round(
                np.mean(within_vals) / (np.mean(between_vals) + 1e-9), 4
            ) if within_vals else 0,
        })

    return pd.DataFrame(rows).sort_values("separation_ratio", ascending=False)


# ── Plot styling ──────────────────────────────────────────────────────────────

_GROUP_COLOURS: dict[str, str] = {
    "GPCR_class_A": "#1565C0",
    "GPCR_class_B": "#42A5F5",
    "GPCR_class_C": "#90CAF9",
    "AGC_kinase_PKA_like": "#B71C1C",
    "AGC_kinase_other": "#E53935",
    "AKT_kinase": "#EF9A9A",
    "CaMK": "#E65100",
    "MAP2K": "#1B5E20",
    "MAPK": "#388E3C",
    "MAPK_related_kinase_or_review": "#66BB6A",
    "MAPK_scaffold": "#A5D6A7",
    "RAF_kinase": "#81C784",
    "PAK_kinase": "#6D4C41",
    "ROCK_kinase": "#A1887F",
    "PI3K_catalytic": "#006064",
    "PI3K_regulatory": "#4DD0E1",
    "ABC_channel": "#4A148C",
    "ABC_transporter": "#7B1FA2",
    "HCN_channel": "#CE93D8",
    "RapGEF_EPAC": "#F57F17",
    "Rho_GEF": "#FBC02D",
    "adenylyl_cyclase": "#AD1457",
    "cAMP_effector_POPDC": "#546E7A",
    "phosphodiesterase": "#827717",
    "heterotrimeric_G_alpha": "#00838F",
    "small_GTPase_Rac_family": "#00695C",
    "small_GTPase_other": "#4DB6AC",
    "calcium_pump": "#1A237E",
    "calmodulin": "#3949AB",
    "voltage_gated_calcium_channel": "#5C6BC0",
    "cyclic_nucleotide_gated_channel": "#9FA8DA",
    "ryanodine_receptor": "#283593",
    "phospholipase_C": "#4E342E",
    "phospholipase_D": "#8D6E63",
    "protein_phosphatase_PP1": "#607D8B",
    "peptide_ligand": "#BF360C",
    "secreted_binding_protein": "#FF7043",
    "secreted_growth_factor": "#FFAB91",
    "transcription_factor": "#311B92",
    "transcriptional_coactivator": "#7E57C2",
    "other_or_review": "#BDBDBD",
    "other": "#9E9E9E",
}

_GROUP_MARKERS: dict[str, str] = {
    "GPCR_class_A": "o",
    "GPCR_class_B": "o",
    "GPCR_class_C": "o",
    "AGC_kinase_PKA_like": "s",
    "AGC_kinase_other": "s",
    "AKT_kinase": "s",
    "CaMK": "s",
    "MAP2K": "^",
    "MAPK": "^",
    "MAPK_related_kinase_or_review": "^",
    "MAPK_scaffold": "^",
    "RAF_kinase": "^",
    "PAK_kinase": "P",
    "ROCK_kinase": "P",
    "PI3K_catalytic": "D",
    "PI3K_regulatory": "D",
    "ABC_channel": "v",
    "ABC_transporter": "v",
    "HCN_channel": "v",
    "RapGEF_EPAC": "*",
    "Rho_GEF": "*",
    "adenylyl_cyclase": "h",
    "cAMP_effector_POPDC": "X",
    "phosphodiesterase": "8",
    "heterotrimeric_G_alpha": "^",
    "small_GTPase_Rac_family": ">",
    "small_GTPase_other": ">",
    "calcium_pump": "v",
    "calmodulin": "v",
    "voltage_gated_calcium_channel": "v",
    "cyclic_nucleotide_gated_channel": "v",
    "ryanodine_receptor": "v",
    "phospholipase_C": "p",
    "phospholipase_D": "p",
    "protein_phosphatase_PP1": "H",
    "peptide_ligand": "d",
    "secreted_binding_protein": "d",
    "secreted_growth_factor": "d",
    "transcription_factor": "+",
    "transcriptional_coactivator": "1",
    "other_or_review": "o",
    "other": "o",
}

_SUPERFAMILY_ORDER: list[tuple[str, list[str]]] = [
    ("GPCRs", ["GPCR_class_A", "GPCR_class_B", "GPCR_class_C"]),
    ("AGC / PKA kinases", ["AGC_kinase_PKA_like", "AGC_kinase_other", "AKT_kinase", "CaMK"]),
    ("MAPK cascade", ["MAP2K", "MAPK", "MAPK_related_kinase_or_review", "MAPK_scaffold", "RAF_kinase"]),
    ("Other kinases", ["PAK_kinase", "ROCK_kinase"]),
    ("PI3K", ["PI3K_catalytic", "PI3K_regulatory"]),
    ("Channels/transporters", ["ABC_channel", "ABC_transporter", "HCN_channel"]),
    ("GEF / RAS", ["RapGEF_EPAC", "Rho_GEF"]),
    ("Adenylyl cyclase", ["adenylyl_cyclase"]),
    ("cAMP effectors", ["cAMP_effector_POPDC"]),
    ("Phosphodiesterases", ["phosphodiesterase"]),
    ("G proteins", ["heterotrimeric_G_alpha", "small_GTPase_Rac_family", "small_GTPase_other"]),
    ("Calcium signalling", ["calcium_pump", "calmodulin", "voltage_gated_calcium_channel", "cyclic_nucleotide_gated_channel", "ryanodine_receptor"]),
    ("Phospholipases", ["phospholipase_C", "phospholipase_D"]),
    ("Phosphatases", ["protein_phosphatase_PP1"]),
    ("Peptide ligands", ["peptide_ligand", "secreted_binding_protein", "secreted_growth_factor"]),
    ("Transcription", ["transcription_factor", "transcriptional_coactivator"]),
    ("Other", ["other_or_review", "other"]),
]


def _get_group_colour(group: str) -> str:
    if group in _GROUP_COLOURS:
        return _GROUP_COLOURS[group]
    import hashlib
    idx = int(hashlib.md5(group.encode()).hexdigest(), 16) % 20
    return plt.cm.tab20(idx / 20)


def _mmseqs2_border_colours(n_clusters: int) -> list[str]:
    base = [
        "#000000", "#FF6B00", "#00A896", "#9B2335", "#4B0082",
        "#228B22", "#8B4513", "#2F4F4F", "#DC143C", "#006400",
        "#00008B", "#8B008B", "#FF1493", "#FF8C00", "#1E90FF",
        "#32CD32", "#FFD700", "#FF4500", "#9400D3", "#00CED1",
    ]
    return [base[i % len(base)] for i in range(n_clusters)]


# ── Generic ordination plot (PCA / PCoA) ─────────────────────────────────────

def _ordination_plot(
    coords_df: pd.DataFrame,
    group_series: pd.Series,
    out_path: Path,
    mmseqs2_series: pd.Series | None = None,
    dpi: int = 150,
    font_size: int = 9,
    ordination_name: str = "PCA",
) -> None:
    """
    Generic ordination plot for PCA or PCoA.
    First two columns of coords_df are used as x/y axes.
    """
    if coords_df.empty or coords_df.shape[1] < 2:
        return

    x_col = coords_df.columns[0]
    y_col = coords_df.columns[1]

    has_mmseqs2 = mmseqs2_series is not None and len(mmseqs2_series) > 0
    present_groups = sorted(group_series.unique())

    if has_mmseqs2:
        fig, (ax, ax2) = plt.subplots(1, 2, figsize=(22, 9), dpi=dpi)
    else:
        fig, ax = plt.subplots(figsize=(12, 9), dpi=dpi)
        ax2 = None

    if has_mmseqs2:
        unique_clusters = sorted(
            mmseqs2_series.unique(),
            key=lambda x: int(x) if str(x).isdigit() else str(x)
        )
        border_palette = _mmseqs2_border_colours(len(unique_clusters))
        cluster_to_border = {c: border_palette[i] for i, c in enumerate(unique_clusters)}
    else:
        unique_clusters = []
        cluster_to_border = {}

    # panel 1
    for grp in present_groups:
        idx = group_series[group_series == grp].index
        pts = coords_df.loc[coords_df.index.intersection(idx)]
        col = _get_group_colour(grp)
        mrk = _GROUP_MARKERS.get(grp, "o")

        if has_mmseqs2:
            for prot in pts.index:
                sc = mmseqs2_series.get(prot, None)
                border_col = cluster_to_border.get(str(sc), "#cccccc") if sc is not None else "#cccccc"
                border_width = 2.0 if sc is not None else 0.5
                ax.scatter(
                    coords_df.loc[prot, x_col],
                    coords_df.loc[prot, y_col],
                    color=col,
                    marker=mrk,
                    s=95,
                    alpha=0.90,
                    edgecolors=border_col,
                    linewidths=border_width,
                    zorder=3,
                )
        else:
            ax.scatter(
                pts[x_col], pts[y_col],
                color=col, marker=mrk,
                s=75, alpha=0.88,
                edgecolors="white", linewidths=0.5, zorder=3
            )

    if len(coords_df) <= 80:
        for prot in coords_df.index:
            ax.annotate(
                prot,
                xy=(coords_df.loc[prot, x_col], coords_df.loc[prot, y_col]),
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=max(3, font_size - 4),
                alpha=0.50,
                color="#333333",
                zorder=2,
            )

    ax.set_xlabel(x_col, fontsize=font_size)
    ax.set_ylabel(y_col, fontsize=font_size)
    title1 = (
        f"{ordination_name} of Ontology Space (fill=ontology group, border=MMseqs2 cluster)"
        if has_mmseqs2 else
        f"{ordination_name} of Protein Ontology Space"
    )
    ax.set_title(title1, fontsize=font_size + 1, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=font_size - 1)

    leg_handles, leg_labels = [], []
    for superfamily, members in _SUPERFAMILY_ORDER:
        present_m = [g for g in members if g in present_groups]
        if not present_m:
            continue
        leg_handles.append(plt.Line2D([0], [0], linestyle="none", marker="none", alpha=0))
        leg_labels.append(f"── {superfamily} ──")
        for grp in present_m:
            col = _get_group_colour(grp)
            mrk = _GROUP_MARKERS.get(grp, "o")
            leg_handles.append(
                plt.Line2D([0], [0], linestyle="none", marker=mrk,
                           markerfacecolor=col, markeredgecolor="white",
                           markeredgewidth=0.5, markersize=7)
            )
            leg_labels.append(grp)

    covered = {g for _, members in _SUPERFAMILY_ORDER for g in members}
    extras = [g for g in present_groups if g not in covered]
    if extras:
        leg_handles.append(plt.Line2D([0], [0], linestyle="none", marker="none", alpha=0))
        leg_labels.append("── Unclassified ──")
        for grp in extras:
            col = _get_group_colour(grp)
            leg_handles.append(
                plt.Line2D([0], [0], linestyle="none", marker="o",
                           markerfacecolor=col, markeredgecolor="white", markersize=7)
            )
            leg_labels.append(grp)

    if has_mmseqs2:
        leg_handles.append(plt.Line2D([0], [0], linestyle="none", marker="none", alpha=0))
        leg_labels.append("── MMseqs2 clusters (border) ──")
        for clust, bcol in cluster_to_border.items():
            leg_handles.append(
                plt.Line2D([0], [0], linestyle="none", marker="o",
                           markerfacecolor="#aaaaaa",
                           markeredgecolor=bcol,
                           markeredgewidth=2.0, markersize=7)
            )
            leg_labels.append(f"Seq cluster {clust}")

    ax.legend(
        leg_handles, leg_labels,
        fontsize=max(5, font_size - 2),
        bbox_to_anchor=(1.01, 1), loc="upper left",
        frameon=True, framealpha=0.9, edgecolor="#cccccc",
        handlelength=1.2, handletextpad=0.5,
        borderpad=0.6, labelspacing=0.3
    )

    # panel 2: MMseqs2 only
    if has_mmseqs2 and ax2 is not None:
        for clust in unique_clusters:
            m_idx = mmseqs2_series[mmseqs2_series == clust].index
            pts = coords_df.loc[coords_df.index.intersection(m_idx)]
            if pts.empty:
                continue
            col = cluster_to_border[clust]
            ax2.scatter(
                pts[x_col], pts[y_col],
                color=col, marker="o", s=75,
                alpha=0.85, edgecolors="white", linewidths=0.4,
                label=f"Seq cluster {clust}", zorder=3
            )

        if len(coords_df) <= 80:
            for prot in coords_df.index:
                ax2.annotate(
                    prot,
                    xy=(coords_df.loc[prot, x_col], coords_df.loc[prot, y_col]),
                    xytext=(3, 3),
                    textcoords="offset points",
                    fontsize=max(3, font_size - 4),
                    alpha=0.50,
                    color="#333333",
                    zorder=2,
                )

        ax2.set_xlabel(x_col, fontsize=font_size)
        ax2.set_ylabel(y_col, fontsize=font_size)
        ax2.set_title(
            f"{ordination_name}: MMseqs2 Sequence Clusters",
            fontsize=font_size + 1, fontweight="bold"
        )
        ax2.spines[["top", "right"]].set_visible(False)
        ax2.tick_params(labelsize=font_size - 1)
        ax2.legend(
            fontsize=max(5, font_size - 2),
            bbox_to_anchor=(1.01, 1), loc="upper left",
            frameon=True, framealpha=0.9, edgecolor="#cccccc",
            handlelength=1.0, labelspacing=0.3
        )

    plt.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# Enhanced clustering visualisation suite
# ══════════════════════════════════════════════════════════════════════════════

def _top_pair_table(
    sim_mat: pd.DataFrame,
    cluster_series: pd.Series,
    group_series: pd.Series,
    feat_mat: pd.DataFrame,
    out_path: Path,
    top_n: int = 30,
) -> pd.DataFrame:
    proteins = sim_mat.index.tolist()
    rows = []

    for i in range(len(proteins)):
        for j in range(i + 1, len(proteins)):
            a, b = proteins[i], proteins[j]
            sim = float(sim_mat.loc[a, b])
            if sim < 0.3:
                continue

            same_clust = (str(cluster_series.get(a, "")) == str(cluster_series.get(b, "")))

            if not feat_mat.empty and a in feat_mat.index and b in feat_mat.index:
                shared = feat_mat.columns[
                    (feat_mat.loc[a] == 1) & (feat_mat.loc[b] == 1)
                ].tolist()
                n_shared = len(shared)
                short = [f.split("::")[-1][:25] for f in shared[:5]]
                shared_str = "; ".join(short)
            else:
                n_shared, shared_str = 0, ""

            rows.append({
                "protein_A": a,
                "protein_B": b,
                "jaccard_similarity": round(sim, 4),
                "same_cluster": same_clust,
                "group_A": str(group_series.get(a, "other")),
                "group_B": str(group_series.get(b, "other")),
                "n_shared_features": n_shared,
                "top_shared_features": shared_str,
            })

    df = pd.DataFrame(rows).sort_values("jaccard_similarity", ascending=False).reset_index(drop=True)
    if top_n and len(df) > top_n:
        df = df.head(top_n)
    df.to_csv(out_path, sep="\t", index=False)
    return df


def _within_between_boxplot(
    sim_mat: pd.DataFrame,
    cluster_series: pd.Series,
    mmseqs2_series: pd.Series | None,
    out_path: Path,
    dpi: int = 150,
    font_size: int = 9,
) -> None:
    """
    Boxplot (+ jitter scatter overlay) of within-cluster vs between-cluster
    Jaccard similarity.

    Root cause of the previous bug
    ───────────────────────────────
    The original _draw_panel helper built per-cluster similarity lists as a
    Python list-of-lists (one sub-list per cluster, variable length) and then
    passed that object directly to ax.scatter as the y argument:

        jitter_y = [cluster_sims[ci] for ci in clusters]   # (k,) + ragged
        ax.scatter(jitter_x, jitter_y, ...)                 # ValueError

    NumPy cannot convert a ragged list-of-lists into a 2-D array, so it raises:
        "setting an array element with a sequence. The requested array has an
         inhomogeneous shape after 1 dimensions. Detected shape: (12,) + ..."

    Fix: always flatten per-cluster lists into 1-D arrays BEFORE passing to
    scatter.  np.concatenate([list_of_lists]) or an explicit loop that calls
    .extend() both produce a flat, homogeneous 1-D sequence that scatter
    accepts.
    """
    from scipy.stats import mannwhitneyu
    rng = np.random.default_rng(42)

    proteins = [p for p in sim_mat.index if p in cluster_series.index]

    # ── collect flat within/between lists for the overall boxplot ─────────
    within_ont, between_ont = [], []
    for i in range(len(proteins)):
        for j in range(i + 1, len(proteins)):
            a, b = proteins[i], proteins[j]
            s = float(sim_mat.loc[a, b])
            if cluster_series[a] == cluster_series[b]:
                within_ont.append(s)
            else:
                between_ont.append(s)

    # ── per-cluster within-similarity for scatter jitter (flat lists) ─────
    clusters_sorted = sorted(cluster_series.unique())
    # jitter_x and jitter_y must be 1-D flat arrays — NO list-of-lists
    jitter_x_ont: list[float] = []
    jitter_y_ont: list[float] = []
    for pos, ci in enumerate(clusters_sorted, start=1):
        members = [p for p in proteins if cluster_series[p] == ci]
        for ii in range(len(members)):
            for jj in range(ii + 1, len(members)):
                a, b = members[ii], members[jj]
                jitter_x_ont.append(
                    1.0 + rng.uniform(-0.18, 0.18)  # jitter around box 1
                )
                jitter_y_ont.append(float(sim_mat.loc[a, b]))

    has_seq = mmseqs2_series is not None and len(mmseqs2_series) > 0
    within_seq, between_seq = [], []
    jitter_x_seq: list[float] = []
    jitter_y_seq: list[float] = []

    if has_seq:
        seq_prots = [p for p in proteins if p in mmseqs2_series.index]
        seq_clusters_sorted = sorted(mmseqs2_series.unique(),
                                     key=lambda x: int(x) if str(x).isdigit() else str(x))
        for i in range(len(seq_prots)):
            for j in range(i + 1, len(seq_prots)):
                a, b = seq_prots[i], seq_prots[j]
                s = float(sim_mat.loc[a, b])
                if str(mmseqs2_series[a]) == str(mmseqs2_series[b]):
                    within_seq.append(s)
                    jitter_x_seq.append(1.0 + rng.uniform(-0.18, 0.18))
                    jitter_y_seq.append(s)
                else:
                    between_seq.append(s)

    n_panels = 2 if has_seq else 1
    fig, axes = plt.subplots(1, n_panels, figsize=(6 * n_panels, 6), dpi=dpi)
    axes_main = [axes] if n_panels == 1 else list(axes)

    def _draw_panel(ax, within, between, jx, jy, title):
        """
        Draw boxplot for overall within/between distributions plus a jitter
        scatter for the within-cluster pairs.

        Parameters
        ----------
        within, between : list[float]   — flat 1-D lists of scalar similarities
        jx, jy          : list[float]   — flat 1-D jitter coordinates (same length)

        The critical invariant: jx and jy must be 1-D flat numeric sequences.
        Passing a list-of-lists (ragged) to ax.scatter raises ValueError.
        """
        data   = [within, between]
        labels = [f"Within\n(n={len(within)})", f"Between\n(n={len(between)})"]
        bp = ax.boxplot(
            data, labels=labels, patch_artist=True,
            medianprops=dict(color="black", linewidth=2),
            whiskerprops=dict(linewidth=1.2),
            flierprops=dict(marker=".", markersize=2, alpha=0.3),
            zorder=2,
        )
        bp["boxes"][0].set_facecolor("#90CAF9")
        bp["boxes"][1].set_facecolor("#EF9A9A")

        # scatter jitter for within-cluster pairs
        # jx and jy are guaranteed flat 1-D — no inhomogeneous-array risk
        if len(jx) > 0:
            jx_arr = np.asarray(jx, dtype=float)   # explicit 1-D float array
            jy_arr = np.asarray(jy, dtype=float)   # same
            ax.scatter(
                jx_arr, jy_arr,
                s=8, alpha=0.35, color="#1565C0",
                edgecolors="none", zorder=3,
            )

        if len(within) > 1 and len(between) > 1:
            _, p = mannwhitneyu(within, between, alternative="greater")
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
            y_max = max(max(within), max(between)) + 0.05
            ax.plot([1, 2], [y_max, y_max], "k-", linewidth=1)
            ax.text(1.5, y_max + 0.01, f"{sig}\np={p:.4f}",
                    ha="center", fontsize=font_size - 1)

        ax.set_ylabel("Jaccard Similarity", fontsize=font_size)
        ax.set_title(title, fontsize=font_size + 1, fontweight="bold")
        ax.set_ylim(-0.05, 1.15)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=font_size - 1)

    _draw_panel(
        axes_main[0],
        within_ont, between_ont,
        jitter_x_ont, jitter_y_ont,
        "Ontology Clusters\nWithin vs Between Jaccard",
    )
    if has_seq:
        _draw_panel(
            axes_main[1],
            within_seq, between_seq,
            jitter_x_seq, jitter_y_seq,
            "MMseqs2 Sequence Clusters\nWithin vs Between Jaccard",
        )

    fig.suptitle(
        "Within-cluster similarity should be higher than between-cluster\n"
        "if ontology annotations reflect biologically meaningful groupings",
        fontsize=font_size, style="italic",
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _cluster_summary_heatmap(
    sim_mat: pd.DataFrame,
    cluster_series: pd.Series,
    group_series: pd.Series,
    out_path: Path,
    dpi: int = 150,
    font_size: int = 9,
) -> None:
    clusters = sorted(cluster_series.unique())
    n = len(clusters)
    mat = np.zeros((n, n))
    labels = []

    for i, ci in enumerate(clusters):
        pi = cluster_series[cluster_series == ci].index.tolist()
        pi = [p for p in pi if p in sim_mat.index]

        grp_counts = group_series.loc[[p for p in pi if p in group_series.index]].value_counts()
        dominant = grp_counts.index[0] if not grp_counts.empty else f"C{ci}"
        labels.append(f"C{ci}\n{dominant[:18]}")

        for j, cj in enumerate(clusters):
            pj = cluster_series[cluster_series == cj].index.tolist()
            pj = [p for p in pj if p in sim_mat.index]
            vals = [float(sim_mat.loc[a, b]) for a in pi for b in pj if a != b or ci != cj]
            mat[i, j] = np.mean(vals) if vals else 0.0

    fig, ax = plt.subplots(figsize=(max(8, n * 0.7), max(7, n * 0.65)), dpi=dpi)
    im = ax.imshow(mat, cmap="YlOrRd", vmin=0, vmax=1, aspect="auto", interpolation="nearest")
    plt.colorbar(im, ax=ax, shrink=0.7, label="Mean Jaccard Similarity")

    for i in range(n):
        for j in range(n):
            val = mat[i, j]
            colour = "white" if val > 0.6 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=max(5, font_size - 2), color=colour)

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=font_size - 1)
    ax.set_yticklabels(labels, fontsize=font_size - 1)
    ax.set_title(
        "Cluster-Level Mean Jaccard Similarity\n"
        "(diagonal = within-cluster coherence, off-diagonal = between-cluster overlap)",
        fontsize=font_size + 1, fontweight="bold"
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _per_cluster_heatmaps(
    sim_mat: pd.DataFrame,
    cluster_series: pd.Series,
    group_series: pd.Series,
    out_dir: Path,
    dpi: int = 150,
    font_size: int = 8,
) -> None:
    clusters = sorted(cluster_series.unique())
    for ci in clusters:
        members = cluster_series[cluster_series == ci].index.tolist()
        members = [p for p in members if p in sim_mat.index]
        if len(members) < 2:
            continue

        mat = sim_mat.loc[members, members].values
        labels = [f"{p}\n({str(group_series.get(p,''))[:16]})" for p in members]
        n = len(members)
        size = max(4, n * 0.55)

        fig, ax = plt.subplots(figsize=(size + 1.5, size), dpi=dpi)
        im = ax.imshow(mat, cmap="YlOrRd", vmin=0, vmax=1, aspect="auto", interpolation="nearest")
        plt.colorbar(im, ax=ax, shrink=0.8, label="Jaccard")

        for i in range(n):
            for j in range(n):
                val = mat[i, j]
                col = "white" if val > 0.65 else "black"
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=max(5, font_size - 1), color=col)

        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=max(5, font_size - 1))
        ax.set_yticklabels(labels, fontsize=max(5, font_size - 1))

        grp_counts = group_series.loc[[p for p in members if p in group_series.index]].value_counts()
        dom_grp = grp_counts.index[0] if not grp_counts.empty else ""
        ax.set_title(f"Cluster {ci}  —  {dom_grp}  ({n} proteins)",
                     fontsize=font_size + 1, fontweight="bold")

        plt.tight_layout()
        fig.savefig(out_dir / f"cluster_{ci}_heatmap.png", dpi=dpi, bbox_inches="tight")
        plt.close(fig)


def _clustering_heatmap(
    sim_mat: pd.DataFrame,
    cluster_series: pd.Series,
    group_series: pd.Series,
    out_path: Path,
    dpi: int = 150,
    font_size: int = 7,
) -> None:
    if sim_mat.empty:
        return

    dist_vec = squareform(1 - sim_mat.values.clip(0, 1), checks=False)
    link = linkage(dist_vec, method="ward")
    dend = dendrogram(link, no_plot=True)
    order = [sim_mat.index[i] for i in dend["leaves"]]
    mat = sim_mat.loc[order, order].values.astype(float)
    n = len(order)

    mask = np.tril(np.ones((n, n), dtype=bool), k=-1)
    mat_masked = np.where(mask, np.nan, mat)

    fig, ax = plt.subplots(figsize=(12, 10), dpi=dpi)
    cmap = plt.cm.YlOrRd.copy()
    cmap.set_bad(color="#f5f5f5")
    im = ax.imshow(mat_masked, cmap=cmap, vmin=0, vmax=1, aspect="auto", interpolation="nearest")
    plt.colorbar(im, ax=ax, shrink=0.55, label="Jaccard Similarity", pad=0.02)

    if n <= 50:
        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels(order, rotation=90, fontsize=max(4, font_size - 1))
        ax.set_yticklabels(order, fontsize=max(4, font_size - 1))
    else:
        ax.set_xticks([])
        ax.set_yticks([])

    clust_order = [cluster_series.get(p, 0) for p in order]
    boundaries = []
    prev = None
    for i, c in enumerate(clust_order):
        if c != prev and prev is not None:
            boundaries.append(i - 0.5)
        prev = c

    for b in boundaries:
        ax.axhline(b, color="#222222", linewidth=1.2, zorder=5)
        ax.axvline(b, color="#222222", linewidth=1.2, zorder=5)

    block_starts = [0] + [int(b + 0.5) for b in boundaries]
    block_ends = [int(b + 0.5) for b in boundaries] + [n]
    for start, end in zip(block_starts, block_ends):
        mid = (start + end) / 2
        block_prots = order[start:end]
        grp_counts = group_series.loc[[p for p in block_prots if p in group_series.index]].value_counts()
        label = grp_counts.index[0][:20] if not grp_counts.empty else ""
        ax.text(-0.8, mid, label, ha="right", va="center",
                fontsize=max(4, font_size - 1), color="#333333",
                transform=ax.get_yaxis_transform())

    if n <= 40:
        for i in range(n):
            for j in range(i, n):
                val = mat[i, j]
                if val > 0.70 and i != j:
                    ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                            fontsize=max(4, font_size - 2),
                            color="white" if val > 0.8 else "#333333")

    ax.set_title(
        "Protein Ontology Similarity Matrix\n"
        "(hierarchical clustering order, upper triangle, block labels = dominant biological group)",
        fontsize=font_size + 1, fontweight="bold"
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# ── Stage runner ──────────────────────────────────────────────────────────────

def run(cfg: dict) -> None:
    logger = setup_logging(cfg.get("log_dir", "logs"), "stage13_protein_modules")
    stage_banner(logger, "STAGE 13", "Protein Module Separation & Clustering")
    out_dir = make_output_dir(cfg, "stage13_modules")

    occ_path = get_input_path(cfg, "domain_occurrence_table")
    if not occ_path:
        logger.warning("domain_occurrence_table not found — skipping Stage 13")
        return

    occ = load_table(occ_path)
    occ = rename_to_canonical(occ, cfg.get("column_aliases", {}))

    cluster_cfg = cfg.get("clustering", {})
    plot_cfg = cfg.get("plotting", {})
    dpi = plot_cfg.get("figure_dpi", 150)
    fs = plot_cfg.get("font_size", 9)
    link_method = cluster_cfg.get("linkage_method", "ward")

    # auto-select k
    k_hint = cluster_cfg.get("n_clusters", 6)
    max_k = cluster_cfg.get("max_clusters", 12)
    auto_k = cluster_cfg.get("auto_select_k", True)

    feature_types = cluster_cfg.get(
        "feature_types",
        ["domain_types", "function_tags", "topology_context", "binding_partner"]
    )

    feat_col_map = {
        "domain_types": "domain_type_id",
        "function_tags": "function_tags",
        "topology_context": "topology_context",
        "binding_partner": next(
            (c for c in ["binding_partner", "binding_parnter"] if c in occ.columns),
            "binding_partner",
        ),
        "proximity_categories": "proximity_categories",
        "process": "process",
    }

    feat_cols = [feat_col_map[k] for k in feature_types if feat_col_map.get(k) in occ.columns]
    if not feat_cols:
        logger.warning("No feature columns found — check clustering.feature_types in config")
        return

    prot_col = "protein_accession"
    if prot_col not in occ.columns:
        logger.error(f"'{prot_col}' column not found")
        return

    logger.info(f"Feature columns: {feat_cols}")
    feat_mat = build_feature_matrix(occ, feat_cols, prot_col, binary=True)
    logger.info(f"Feature matrix: {feat_mat.shape} (proteins × features)")

    if feat_mat.empty:
        logger.warning("Feature matrix is empty — skipping Stage 13")
        return

    save_tsv(feat_mat.reset_index().rename(columns={"index": prot_col}),
             out_dir / "protein_feature_matrix.tsv", logger)

    # Similarity
    sim_mat = compute_similarity(feat_mat, metric="jaccard")
    save_tsv(sim_mat.reset_index().rename(columns={"index": prot_col}),
             out_dir / "protein_similarity_matrix.tsv", logger)

    # PCA
    pca_coords, pca_var = run_pca(feat_mat, n_components=min(10, len(feat_mat) - 1))
    if not pca_coords.empty:
        save_tsv(pca_coords.reset_index().rename(columns={"index": prot_col}),
                 out_dir / "pca_coordinates.tsv", logger)
        save_tsv(pca_var, out_dir / "pca_variance_explained.tsv", logger)

    # PCoA
    logger.info("Running PCoA on Jaccard distance matrix …")
    pcoa_coords = pd.DataFrame()
    pcoa_var = pd.DataFrame()
    try:
        pcoa_coords, pcoa_var = run_pcoa(sim_mat, n_components=2)
        if not pcoa_coords.empty:
            save_tsv(pcoa_coords.reset_index().rename(columns={"index": prot_col}),
                     out_dir / "pcoa_coordinates.tsv", logger)
            save_tsv(pcoa_var, out_dir / "pcoa_variance_explained.tsv", logger)
            gof = pcoa_coords.attrs.get("gof_2d", None)
            if gof is not None:
                logger.info(f"  PCoA 2D goodness-of-fit={gof:.4f}")
    except Exception as e:
        logger.warning(f"PCoA failed: {e}")

    # Choose k by silhouette
    if auto_k and len(feat_mat) >= 6:
        upper = min(max_k + 1, max(4, len(feat_mat) // 2))
        k_range = range(3, upper)
        k_scores: list[tuple[int, float]] = []

        for k in k_range:
            try:
                lbls = run_clustering(feat_mat, k, link_method)
                if lbls.nunique() < 2:
                    continue
                sil_k = silhouette_score(feat_mat.values.astype(float), lbls.values)
                k_scores.append((k, round(float(sil_k), 4)))
            except Exception:
                pass

        if k_scores:
            best_k, best_sil = max(k_scores, key=lambda x: x[1])
            n_clusters = best_k
            logger.info(f"Optimal k selection: tested k={[k for k, _ in k_scores]}")
            logger.info(f"  Silhouette scores: {k_scores}")
            logger.info(f"  Best k={best_k}  (silhouette={best_sil:.4f})")
        else:
            n_clusters = k_hint
            logger.warning(f"k selection failed — using hint k={n_clusters}")

        save_tsv(pd.DataFrame(k_scores, columns=["k", "silhouette_score"]),
                 out_dir / "k_selection_silhouette.tsv", logger)
    else:
        n_clusters = k_hint
        logger.info(f"Using configured n_clusters={n_clusters} (auto_select_k=False)")

    # Clustering
    cluster_series = run_clustering(feat_mat, n_clusters, link_method)
    cluster_df = cluster_series.reset_index()
    cluster_df.columns = [prot_col, "cluster_id"]
    save_tsv(cluster_df, out_dir / "clustering_assignments.tsv", logger)

    # Silhouette for chosen k
    sil = None
    if len(feat_mat) > n_clusters:
        try:
            sil = silhouette_score(feat_mat.values.astype(float), cluster_series.values)
            logger.info(f"Silhouette score (k={n_clusters}): {sil:.4f}")
        except Exception as e:
            logger.warning(f"Silhouette score failed: {e}")

    # Group label loading
    grp_df = None
    grp_path = get_input_path(cfg, "protein_group_labels")
    if grp_path and Path(grp_path).exists():
        grp_df = load_table(grp_path)
        grp_col = next((c for c in ["protein_group", "group", "protein_family"] if c in grp_df.columns), None)
        p_col2 = next((c for c in ["protein_accession", "accession"] if c in grp_df.columns), None)
        if grp_col and p_col2:
            logger.info(f"Protein group labels loaded: {grp_df[grp_col].nunique()} groups for {len(grp_df)} proteins from {grp_path}")
        else:
            logger.warning("protein_group_labels.tsv found but missing accession/group columns — falling back to biological_groups keyword matching")
            grp_df = None
    elif grp_path:
        logger.warning(f"protein_group_labels path configured but file not found: {grp_path}")
        grp_df = None
    else:
        logger.info("No protein_group_labels configured — using biological_groups keyword matching from config")

    grp_col = next((c for c in ["protein_group", "group", "protein_family"] if grp_df is not None and c in grp_df.columns), None)
    p_col2 = next((c for c in ["protein_accession", "accession"] if grp_df is not None and c in grp_df.columns), None)
    bio_groups = cfg.get("biological_groups", {})

    group_series = pd.Series(
        {p: _assign_group(p, bio_groups, grp_df, p_col2, grp_col) for p in feat_mat.index},
        name="biological_group",
    )
    n_assigned = (group_series != "other").sum()
    logger.info(f"Group assignments: {n_assigned}/{len(group_series)} proteins assigned to named groups ({group_series.nunique()} unique groups)")

    group_df = group_series.reset_index()
    group_df.columns = [prot_col, "biological_group"]
    cluster_df = cluster_df.merge(group_df, on=prot_col, how="left")

    # Within / between known group similarity
    sim_stats = _group_similarity_stats(sim_mat, group_series)
    save_tsv(sim_stats, out_dir / "group_similarity_stats.tsv", logger)

    # MMseqs2 loading and comparison
    mmseqs2_series: pd.Series | None = None
    ari_val = None
    nmi_val = None

    mmseqs2_path = get_input_path(cfg, "mmseqs2_clusters")
    if mmseqs2_path and Path(mmseqs2_path).exists():
        try:
            mdf = load_table(mmseqs2_path)
            acc_c = next((c for c in ["accession", "protein_accession"] if c in mdf.columns), None)
            cid_c = next((c for c in ["cluster_id", "cluster"] if c in mdf.columns), None)

            if acc_c and cid_c:
                mmseqs2_series = pd.Series(
                    mdf.set_index(acc_c)[cid_c].astype(str),
                    name="mmseqs2_cluster",
                )
                n_overlap = len(set(mmseqs2_series.index) & set(feat_mat.index))
                logger.info(f"MMseqs2 clusters loaded: {mmseqs2_series.nunique()} clusters, {len(mmseqs2_series)} proteins, {n_overlap} overlap with feature matrix")

                if n_overlap == 0:
                    logger.warning("MMseqs2 accessions do not overlap with the occurrence table — check accession format")
                    mmseqs2_series = None
                else:
                    common = sorted(set(feat_mat.index) & set(mmseqs2_series.index))
                    if len(common) >= 4:
                        ont = cluster_series.loc[common].astype(str).values
                        seq = mmseqs2_series.loc[common].astype(str).values
                        ari_val = round(float(adjusted_rand_score(ont, seq)), 4)
                        nmi_val = round(float(normalized_mutual_info_score(ont, seq, average_method="arithmetic")), 4)
                        logger.info(f"ARI vs MMseqs2 clusters: {ari_val:.4f}  NMI: {nmi_val:.4f}")
            else:
                logger.warning(
                    "MMseqs2 file missing required columns.\n"
                    f"  Found columns: {mdf.columns.tolist()}\n"
                    "  Required: accession/protein_accession and cluster_id/cluster"
                )
        except Exception as e:
            logger.warning(f"Could not load MMseqs2 clusters: {e}")
    elif mmseqs2_path:
        logger.info(f"MMseqs2 cluster file not found: {mmseqs2_path}")
    else:
        logger.info("No mmseqs2_clusters configured — sequence comparison skipped")

    # Summary (after MMseqs2 metrics are known)
    summary_rows = [
        {"metric": "n_proteins", "value": len(feat_mat)},
        {"metric": "n_features", "value": feat_mat.shape[1]},
        {"metric": "n_clusters_chosen", "value": n_clusters},
        {"metric": "silhouette_score", "value": round(sil, 4) if sil is not None else "N/A"},
        {"metric": "ARI_vs_mmseqs2", "value": ari_val if ari_val is not None else "N/A"},
        {"metric": "NMI_vs_mmseqs2", "value": nmi_val if nmi_val is not None else "N/A"},
        {"metric": "pcoa_gof_2d", "value": round(pcoa_coords.attrs.get("gof_2d", 0), 4) if not pcoa_coords.empty else "N/A"},
    ]
    save_tsv(pd.DataFrame(summary_rows), out_dir / "module_separation_summary.tsv", logger)

    # figures / tables
    logger.info("Building top protein-pair similarity table …")
    _top_pair_table(sim_mat, cluster_series, group_series, feat_mat,
                    out_dir / "top_similar_pairs.tsv", top_n=50)

    logger.info("Generating global clustering heatmap …")
    _clustering_heatmap(sim_mat, cluster_series, group_series,
                        out_dir / "clustering_heatmap.png", dpi, fs)

    logger.info("Generating cluster-level summary heatmap …")
    _cluster_summary_heatmap(sim_mat, cluster_series, group_series,
                             out_dir / "cluster_summary_heatmap.png", dpi, fs)

    logger.info("Generating per-cluster sub-heatmaps …")
    sub_dir = out_dir / "per_cluster_heatmaps"
    sub_dir.mkdir(exist_ok=True)
    _per_cluster_heatmaps(sim_mat, cluster_series, group_series,
                          sub_dir, dpi, fs)

    logger.info("Generating within vs between cluster boxplot …")
    _within_between_boxplot(sim_mat, cluster_series, mmseqs2_series,
                            out_dir / "within_between_similarity.png", dpi, fs)

    if not pca_coords.empty:
        logger.info("Generating PCA plot …")
        _ordination_plot(
            pca_coords,
            group_series,
            out_path=out_dir / "pca_plot.png",
            mmseqs2_series=mmseqs2_series,
            dpi=dpi,
            font_size=fs,
            ordination_name="PCA",
        )

    if not pcoa_coords.empty:
        logger.info("Generating PCoA plot (primary ordination on Jaccard distances) …")
        _ordination_plot(
            pcoa_coords,
            group_series,
            out_path=out_dir / "pcoa_plot.png",
            mmseqs2_series=mmseqs2_series,
            dpi=dpi,
            font_size=fs,
            ordination_name="PCoA",
        )

    k_tsv = out_dir / "k_selection_silhouette.tsv"
    if k_tsv.exists():
        try:
            k_df = pd.read_csv(k_tsv, sep="\t")
            if not k_df.empty:
                fig_k, ax_k = plt.subplots(figsize=(7, 4), dpi=dpi)
                ax_k.plot(k_df["k"], k_df["silhouette_score"],
                          "o-", color="#1565C0", linewidth=2, markersize=7)
                best_row = k_df.loc[k_df["silhouette_score"].idxmax()]
                ax_k.axvline(best_row["k"], color="#B71C1C", linestyle="--",
                             linewidth=1.5,
                             label=f"Best k={int(best_row['k'])} (sil={best_row['silhouette_score']:.3f})")
                ax_k.set_xlabel("Number of clusters (k)", fontsize=fs)
                ax_k.set_ylabel("Silhouette Score", fontsize=fs)
                ax_k.set_title("Optimal k via Silhouette Score", fontsize=fs + 1, fontweight="bold")
                ax_k.legend(fontsize=fs - 1)
                ax_k.spines[["top", "right"]].set_visible(False)
                plt.tight_layout()
                fig_k.savefig(out_dir / "k_selection_plot.png", dpi=dpi, bbox_inches="tight")
                plt.close(fig_k)
        except Exception as e:
            logger.warning(f"k-selection plot failed: {e}")

    logger.info("Stage 13 complete.")
    logger.info("")
    logger.info("Output files:")
    logger.info("  clustering_heatmap.png         — global heatmap (dendrogram order)")
    logger.info("  cluster_summary_heatmap.png    — cluster × cluster mean Jaccard")
    logger.info("  per_cluster_heatmaps/          — one readable heatmap per cluster")
    logger.info("  within_between_similarity.png  — statistical boxplot (Mann–Whitney)")
    logger.info("  top_similar_pairs.tsv          — ranked table of similar pairs")
    logger.info("  pcoa_plot.png                  — PCoA on Jaccard distance (primary)")
    logger.info("  pca_plot.png                   — PCA on feature variance")
    logger.info("  k_selection_plot.png           — silhouette vs k")
    logger.info("  k_selection_silhouette.tsv     — silhouette scores per k")


if __name__ == "__main__":
    import sys
    import os
    _project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(_project_root)
    _cfg_path = sys.argv[1] if len(sys.argv) > 1 else "config/pipeline_config.yaml"
    run(load_config(_cfg_path))