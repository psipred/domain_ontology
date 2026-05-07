"""
scripts/stage16_sequence_ontology_validation.py
────────────────────────────────────────────────────────────────────────────
STAGE 16 — SEQUENCE vs ONTOLOGY CONGRUENCE ANALYSIS

Scientific question
───────────────────
Do the domain/function/topology annotations in this ontology reflect
the underlying evolutionary relationships between proteins?

Evidence hierarchy (strongest → weakest)
─────────────────────────────────────────
1. Mantel test          — correlation between two independent distance
                          matrices (sequence vs ontology). Permutation-based
                          p-value. This is the primary statistical test.

2. ARI / NMI            — agreement between hard cluster assignments.
                          Adjusted Rand Index corrects for chance.
                          Normalised Mutual Information measures shared
                          information between partitions.

3. Enrichment test      — for each MMseqs2 sequence cluster, are specific
                          ontology terms (domain types, function tags)
                          over-represented relative to background?
                          Fisher's exact test + Benjamini-Hochberg FDR.

4. PCA / PCoA           — visualisation only. Positions proteins in 2D
                          using ontology distances (PCoA) or feature
                          variance (PCA). Colour by BOTH ontology group
                          and MMseqs2 cluster. Congruence is visually
                          confirmed but NOT statistically inferred from
                          this plot.

Inputs required
───────────────
  domain_occurrence_table   — occurrence TSV with domain/function annotations
  mmseqs2_clusters          — TSV from mmseqs2_cluster.sh:
                              columns: accession | cluster_id | representative
  protein_group_labels      — optional TSV: accession | protein_group

Outputs
───────
  mantel_test_result.tsv        — r statistic, p-value, n_permutations
  ari_nmi_scores.tsv            — ARI, NMI, per threshold
  enrichment_results.tsv        — per-cluster per-term enrichment
  distance_matrix_ontology.tsv  — protein × protein ontology distances
  distance_matrix_sequence.tsv  — protein × protein seq identity distances
  pcoa_coordinates.tsv          — MDS/PCoA coords for plotting
  validation_summary.tsv        — one-row summary of all key metrics
  congruence_plot.png           — 4-panel figure (distance scatter,
                                  ARI bar, enrichment heatmap, PCoA)
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.spatial.distance import pdist, squareform
from scipy.stats import fisher_exact, pearsonr, spearmanr
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.manifold import MDS

import sys as _sys
import os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from scripts.utils import (
    build_cv_lookup, get_input_path, load_config, load_table,
    make_output_dir, rename_to_canonical, save_tsv,
    setup_logging, split_cell, stage_banner,
)

warnings.filterwarnings("ignore", category=FutureWarning)


# ══════════════════════════════════════════════════════════════════════════════
# 1 — Feature matrix
# ══════════════════════════════════════════════════════════════════════════════

def build_feature_matrix(occ: pd.DataFrame,
                          feature_cols: list[str],
                          prot_col: str = "protein_accession") -> pd.DataFrame:
    """
    Binary protein × ontology-feature matrix.
    Each cell = 1 if protein has that feature in ≥1 occurrence, else 0.
    """
    from collections import defaultdict
    counts: dict[str, dict[str, int]] = defaultdict(dict)
    for _, row in occ.iterrows():
        prot = str(row.get(prot_col, "")).strip()
        if not prot or prot == "nan":
            continue
        for col in feature_cols:
            for tok in split_cell(row.get(col, "")):
                key = f"{col}::{tok}"
                counts[prot][key] = counts[prot].get(key, 0) + 1
    df = pd.DataFrame(counts).T.fillna(0)
    return (df > 0).astype(int)


# ══════════════════════════════════════════════════════════════════════════════
# 2 — Distance matrices
# ══════════════════════════════════════════════════════════════════════════════

def ontology_distance_matrix(feat_mat: pd.DataFrame,
                              metric: str = "jaccard") -> pd.DataFrame:
    """
    Protein × protein distance matrix from the ontology feature matrix.
    Jaccard distance = 1 − Jaccard similarity.
    """
    if feat_mat.empty:
        return pd.DataFrame()
    dist = squareform(pdist(feat_mat.values.astype(float), metric=metric))
    return pd.DataFrame(dist, index=feat_mat.index, columns=feat_mat.index)


def sequence_distance_matrix(mmseqs_df: pd.DataFrame,
                              proteins: list[str],
                              acc_col: str = "accession",
                              cid_col: str = "cluster_id") -> pd.DataFrame:
    """
    Build a proxy sequence distance matrix from MMseqs2 cluster assignments.

    MMseqs2 easy-cluster does not directly output pairwise sequence identities
    unless you run easy-search instead.  We use cluster membership as a proxy:
        distance(A, B) = 0.0   if A and B are in the same cluster
        distance(A, B) = 1.0   if A and B are in different clusters

    This is a conservative approximation — within-cluster pairs may still
    differ in sequence, but it cleanly separates same-cluster from
    cross-cluster comparisons for the Mantel test.

    If you ran MMseqs2 easy-search with --format-output
    'query,target,pident' you can provide a richer identity matrix instead
    by calling sequence_distance_matrix_from_identity().
    """
    clust = pd.Series(
        mmseqs_df.set_index(acc_col)[cid_col].astype(str)
    )
    common = [p for p in proteins if p in clust.index]
    n = len(common)
    mat = np.ones((n, n), dtype=float)
    np.fill_diagonal(mat, 0.0)
    for i in range(n):
        for j in range(i + 1, n):
            if clust[common[i]] == clust[common[j]]:
                mat[i, j] = mat[j, i] = 0.0
    return pd.DataFrame(mat, index=common, columns=common)


def sequence_distance_matrix_from_identity(
        search_tsv: str | Path,
        proteins: list[str],
        query_col: str = "query",
        target_col: str = "target",
        pident_col: str = "pident") -> pd.DataFrame:
    """
    Build a pairwise distance matrix from MMseqs2 easy-search output.
    distance = 1 - (pident / 100)

    To generate the input file on the HPC:
        mmseqs easy-search camp_proteins.fasta camp_proteins.fasta \\
            search_results.tsv tmp \\
            --format-output 'query,target,pident' \\
            --min-seq-id 0.0 -c 0.0
    """
    df = pd.read_csv(search_tsv, sep="\t",
                     names=[query_col, target_col, pident_col])
    # Extract bare accessions from full FASTA headers if needed
    for col in [query_col, target_col]:
        df[col] = df[col].apply(
            lambda x: x.split("|")[1] if "|" in str(x) else str(x).split()[0]
        )
    common = [p for p in proteins
              if p in df[query_col].values or p in df[target_col].values]
    idx = {p: i for i, p in enumerate(common)}
    mat = np.ones((len(common), len(common)), dtype=float)
    np.fill_diagonal(mat, 0.0)
    for _, row in df.iterrows():
        q, t, pid = row[query_col], row[target_col], float(row[pident_col])
        if q in idx and t in idx:
            d = 1.0 - pid / 100.0
            mat[idx[q], idx[t]] = d
            mat[idx[t], idx[q]] = d
    return pd.DataFrame(mat, index=common, columns=common)


# ══════════════════════════════════════════════════════════════════════════════
# 3 — Mantel test
# ══════════════════════════════════════════════════════════════════════════════

def mantel_test(dist_a: pd.DataFrame,
                dist_b: pd.DataFrame,
                n_permutations: int = 9999,
                method: str = "pearson",
                logger=None) -> dict:
    """
    Mantel test: correlation between two distance matrices.

    Returns
    -------
    dict with keys:
        r          — observed correlation coefficient
        p_value    — permutation p-value (two-tailed)
        n_perms    — number of permutations performed
        method     — 'pearson' or 'spearman'
        n_proteins — number of proteins in the analysis
        interpretation — plain-English summary

    The Mantel test is appropriate here because:
      - Both matrices are derived from the SAME set of proteins
      - The observations are not independent (each protein appears in
        n-1 pairs), so standard correlation tests are invalid
      - Permutation of one matrix's row/column order preserves the
        internal structure while breaking the inter-matrix relationship
    """
    # Align to common proteins
    common = sorted(set(dist_a.index) & set(dist_b.index))
    if len(common) < 4:
        return {"r": None, "p_value": None,
                "n_proteins": len(common),
                "error": "Fewer than 4 common proteins"}

    a = dist_a.loc[common, common].values.astype(float)
    b = dist_b.loc[common, common].values.astype(float)
    n = len(common)

    # Extract upper triangle (exclude diagonal)
    tri_idx = np.triu_indices(n, k=1)
    a_flat = a[tri_idx]
    b_flat = b[tri_idx]

    # Observed correlation
    if method == "spearman":
        obs_r, _ = spearmanr(a_flat, b_flat)
    else:
        obs_r, _ = pearsonr(a_flat, b_flat)

    if logger:
        logger.info(f"  Mantel test: n={n} proteins, "
                    f"{len(a_flat)} pairwise distances, "
                    f"observed r={obs_r:.4f}")

    # Permutation distribution
    rng = np.random.default_rng(42)
    perm_rs = np.empty(n_permutations)
    perm_order = np.arange(n)

    for i in range(n_permutations):
        rng.shuffle(perm_order)
        b_perm  = b[np.ix_(perm_order, perm_order)]
        b_flat_ = b_perm[tri_idx]
        if method == "spearman":
            r_, _ = spearmanr(a_flat, b_flat_)
        else:
            r_, _ = pearsonr(a_flat, b_flat_)
        perm_rs[i] = r_

    # Two-tailed p-value
    p_value = (np.sum(np.abs(perm_rs) >= np.abs(obs_r)) + 1) / (n_permutations + 1)

    # Interpretation
    if p_value < 0.001:
        sig = "highly significant (p < 0.001)"
    elif p_value < 0.01:
        sig = "significant (p < 0.01)"
    elif p_value < 0.05:
        sig = "significant (p < 0.05)"
    else:
        sig = "not significant (p >= 0.05)"

    if obs_r > 0.7:
        strength = "strong positive"
    elif obs_r > 0.4:
        strength = "moderate positive"
    elif obs_r > 0.1:
        strength = "weak positive"
    elif obs_r > -0.1:
        strength = "negligible"
    else:
        strength = "negative"

    interp = (
        f"{strength} correlation (r={obs_r:.3f}) between ontology and "
        f"sequence distances — {sig}. "
    )
    if p_value < 0.05 and obs_r > 0.3:
        interp += (
            "The ontology annotations reflect evolutionary relationships. "
            "Proteins that are sequence-similar tend to share ontology features."
        )
    elif p_value < 0.05 and obs_r < 0:
        interp += (
            "Ontology distance is inversely correlated with sequence distance — "
            "possible annotation inconsistency or strong functional convergence."
        )
    elif p_value >= 0.05:
        interp += (
            "No significant matrix correlation. The ontology may capture "
            "functional properties that are independent of sequence similarity, "
            "which can indicate functional convergence or annotation gaps."
        )

    return {
        "r":              round(float(obs_r), 6),
        "p_value":        round(float(p_value), 6),
        "n_permutations": n_permutations,
        "method":         method,
        "n_proteins":     n,
        "n_pairs":        len(a_flat),
        "perm_r_mean":    round(float(perm_rs.mean()), 6),
        "perm_r_std":     round(float(perm_rs.std()), 6),
        "interpretation": interp,
        "_perm_rs":       perm_rs,   # kept for plotting, not saved to TSV
    }


# ══════════════════════════════════════════════════════════════════════════════
# 4 — ARI / NMI
# ══════════════════════════════════════════════════════════════════════════════

def cluster_agreement(ontology_labels: pd.Series,
                      mmseqs2_labels: pd.Series) -> dict:
    """
    Compute Adjusted Rand Index and Normalised Mutual Information between
    two sets of cluster assignments.

    ARI = 1.0  → perfect agreement
    ARI = 0.0  → agreement no better than random chance
    ARI < 0    → less agreement than expected by chance

    NMI = 1.0  → labels share all information
    NMI = 0.0  → labels are independent

    Both metrics require no assumed correspondence between label values,
    i.e. ontology cluster 1 does not need to equal MMseqs2 cluster 1.
    """
    common = sorted(set(ontology_labels.index) & set(mmseqs2_labels.index))
    if len(common) < 2:
        return {"ARI": None, "NMI": None, "n_proteins": len(common)}

    a = ontology_labels.loc[common].astype(str).values
    b = mmseqs2_labels.loc[common].astype(str).values

    ari = adjusted_rand_score(a, b)
    nmi = normalized_mutual_info_score(a, b, average_method="arithmetic")

    # Interpretation
    if ari >= 0.8:
        ari_interp = "near-perfect agreement"
    elif ari >= 0.6:
        ari_interp = "strong agreement"
    elif ari >= 0.4:
        ari_interp = "moderate agreement"
    elif ari >= 0.2:
        ari_interp = "weak agreement"
    elif ari >= 0.0:
        ari_interp = "near-random agreement"
    else:
        ari_interp = "worse than random (possible systematic disagreement)"

    return {
        "ARI":            round(float(ari), 6),
        "NMI":            round(float(nmi), 6),
        "n_proteins":     len(common),
        "n_ontology_clusters": len(set(a)),
        "n_mmseqs2_clusters":  len(set(b)),
        "ARI_interpretation":  ari_interp,
    }


# ══════════════════════════════════════════════════════════════════════════════
# 5 — Enrichment test
# ══════════════════════════════════════════════════════════════════════════════

def enrichment_analysis(feat_mat: pd.DataFrame,
                         mmseqs2_labels: pd.Series,
                         min_count: int = 2,
                         fdr_threshold: float = 0.05) -> pd.DataFrame:
    """
    For each MMseqs2 sequence cluster, test which ontology features
    (domain types, function tags, topology terms) are significantly
    over-represented using Fisher's exact test + BH-FDR correction.

    This answers: "given that these proteins are sequence-similar,
    do they also share specific ontology annotations more than expected?"

    Parameters
    ----------
    feat_mat      : binary protein × feature matrix (from build_feature_matrix)
    mmseqs2_labels: Series mapping protein_accession → cluster_id
    min_count     : minimum number of proteins with a feature to test it
    fdr_threshold : BH-FDR threshold for significance

    Returns
    -------
    DataFrame with columns:
        cluster_id, feature, n_in_cluster, n_total, odds_ratio,
        p_value, fdr, significant
    """
    common = sorted(set(feat_mat.index) & set(mmseqs2_labels.index))
    if not common:
        return pd.DataFrame()

    fm = feat_mat.loc[common]
    labels = mmseqs2_labels.loc[common]
    n_total = len(common)
    clusters = sorted(labels.unique(), key=lambda x: int(x) if str(x).isdigit() else 0)
    features = fm.columns.tolist()

    rows = []
    for clust in clusters:
        in_clust   = labels[labels == clust].index
        out_clust  = labels[labels != clust].index
        n_in = len(in_clust)

        for feat in features:
            a = int(fm.loc[in_clust, feat].sum())   # in cluster, has feature
            b = int(fm.loc[out_clust, feat].sum())  # out of cluster, has feature
            c = n_in - a                             # in cluster, no feature
            d = n_total - n_in - b                   # out of cluster, no feature

            total_with_feat = a + b
            if total_with_feat < min_count:
                continue

            _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
            or_val = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))

            rows.append({
                "cluster_id":      str(clust),
                "feature":         feat,
                "n_cluster":       n_in,
                "n_with_feature_in_cluster":  a,
                "n_with_feature_total":       total_with_feat,
                "odds_ratio":      round(or_val, 4),
                "p_value":         round(p, 6),
            })

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)

    # BH-FDR correction across all tests
    p_vals = df["p_value"].values
    n = len(p_vals)
    order = np.argsort(p_vals)
    fdr = np.empty(n)
    fdr[order] = p_vals[order] * n / (np.arange(n) + 1)
    # Enforce monotonicity
    for i in range(n - 2, -1, -1):
        fdr[order[i]] = min(fdr[order[i]], fdr[order[i + 1]])
    fdr = np.minimum(fdr, 1.0)

    df["fdr"] = np.round(fdr, 6)
    df["significant"] = df["fdr"] < fdr_threshold

    return df.sort_values(["cluster_id", "fdr"]).reset_index(drop=True)


# ══════════════════════════════════════════════════════════════════════════════
# 6 — PCoA (metric MDS on ontology distance matrix)
# ══════════════════════════════════════════════════════════════════════════════

def run_pcoa(dist_mat: pd.DataFrame, n_components: int = 3) -> pd.DataFrame:
    """
    Principal Coordinates Analysis (metric MDS) on the ontology distance
    matrix.  This preserves pairwise distances better than PCA on the raw
    feature matrix and is the standard ordination for distance data.

    Unlike PCA, PCoA does not require linearity assumptions and works
    directly with any distance metric.
    """
    if dist_mat.empty or len(dist_mat) < 3:
        return pd.DataFrame()
    n = min(n_components, len(dist_mat) - 1)
    mds = MDS(n_components=n, dissimilarity="precomputed",
              metric=True, random_state=42, n_init=4)
    coords = mds.fit_transform(dist_mat.values.astype(float))
    df = pd.DataFrame(
        coords,
        index=dist_mat.index,
        columns=[f"PCoA{i+1}" for i in range(n)],
    )
    df.attrs["stress"] = float(mds.stress_)
    return df


# ══════════════════════════════════════════════════════════════════════════════
# 7 — Group colour palettes (reused from stage 13)
# ══════════════════════════════════════════════════════════════════════════════

_GROUP_COLOURS: dict[str, str] = {
    "GPCR_class_A": "#1565C0", "GPCR_class_B": "#42A5F5",
    "GPCR_class_C": "#90CAF9",
    "AGC_kinase_PKA_like": "#B71C1C", "AGC_kinase_other": "#E53935",
    "AKT_kinase": "#EF9A9A", "CaMK": "#E65100",
    "MAP2K": "#1B5E20", "MAPK": "#388E3C",
    "MAPK_related_kinase_or_review": "#66BB6A",
    "MAPK_scaffold": "#A5D6A7", "RAF_kinase": "#81C784",
    "PAK_kinase": "#6D4C41", "ROCK_kinase": "#A1887F",
    "PI3K_catalytic": "#006064", "PI3K_regulatory": "#4DD0E1",
    "ABC_channel": "#4A148C", "ABC_transporter": "#7B1FA2",
    "HCN_channel": "#CE93D8",
    "RapGEF_EPAC": "#F57F17", "Rho_GEF": "#FBC02D",
    "adenylyl_cyclase": "#AD1457",
    "cAMP_effector_POPDC": "#546E7A",
    "other": "#9E9E9E",
}

_MMSEQS2_PALETTE = [
    "#000000","#E6194B","#3CB44B","#4363D8","#F58231",
    "#911EB4","#42D4F4","#F032E6","#BFEF45","#FABED4",
    "#469990","#DCBEFF","#9A6324","#FFFAC8","#800000",
    "#AAFFC3","#808000","#FFD8B1","#000075","#A9A9A9",
]


def _seq_colour(cluster_id: str, unique_clusters: list) -> str:
    idx = unique_clusters.index(cluster_id) if cluster_id in unique_clusters else -1
    return _MMSEQS2_PALETTE[idx % len(_MMSEQS2_PALETTE)] if idx >= 0 else "#cccccc"


def _ont_colour(group: str) -> str:
    if group in _GROUP_COLOURS:
        return _GROUP_COLOURS[group]
    import hashlib
    idx = int(hashlib.md5(group.encode()).hexdigest(), 16) % 20
    return plt.cm.tab20(idx / 20)


# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
# 8 — 4-panel congruence figure  (publication quality)
# ══════════════════════════════════════════════════════════════════════════════

# ── Shared style constants ────────────────────────────────────────────────────
_PANEL_LABEL_KW = dict(fontsize=13, fontweight="bold",
                        transform_rotates_text=False)
_SPINE_COLS     = ["top", "right"]

# Colour palette: colourblind-friendly (Paul Tol muted)
_COL = dict(
    obs_line   = "#CC3311",   # red — observed statistic
    null_hist  = "#6699CC",   # muted blue — permutation distribution
    null_shade = "#DDAAAA",   # light red — critical region shade
    ari        = "#0077BB",   # blue
    nmi        = "#009988",   # teal
    sig_dot    = "#CC3311",   # red dot on heatmap
    insig      = "#F8F0E3",   # near-white — non-significant cell
)


def _label_panel(ax, letter: str, fs: int) -> None:
    """Add bold uppercase panel letter inside top-left of axes."""
    ax.text(-0.10, 1.06, letter, transform=ax.transAxes,
            fontsize=fs + 5, fontweight="bold", va="top", ha="left",
            color="black")


def _congruence_plot(mantel_result: dict,
                     ari_result: dict,
                     enrichment_df: pd.DataFrame,
                     pcoa_coords: pd.DataFrame,
                     group_series: pd.Series,
                     mmseqs2_labels: pd.Series,
                     out_path: Path,
                     dpi: int = 150,
                     font_size: int = 9) -> None:
    """
    Publication-quality 4-panel figure.

    A (top-left)    — Mantel test: one-tailed permutation histogram.
                      Observed r sits outside null distribution.
    B (top-right)   — ARI/NMI with interpretation zones.
    C (bottom-left) — PCoA coloured by ontology superfamily.
    D (bottom-right)— Enrichment heatmap with odds-ratio dot overlay.
    """
    fs = font_size
    fig = plt.figure(figsize=(16, 13), dpi=dpi)
    gs  = gridspec.GridSpec(2, 2, figure=fig,
                             hspace=0.46, wspace=0.40,
                             left=0.09, right=0.97,
                             top=0.91, bottom=0.09)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    _draw_mantel_panel(ax_a, mantel_result, fs)
    _draw_ari_panel(ax_b, ari_result, fs)
    _draw_pcoa_panel(ax_c, pcoa_coords, group_series, mmseqs2_labels, fs)
    _draw_enrichment_panel(ax_d, enrichment_df, fs)

    fig.suptitle(
        "Sequence–Ontology Congruence Validation",
        fontsize=fs + 5, fontweight="bold", y=0.975,
    )
    fig.text(0.5, 0.955,
             "Primary statistics: Mantel test (A) and ARI/NMI (B)  ·  "
             "Visualisation: PCoA (C)  ·  Biological detail: Enrichment (D)",
             ha="center", fontsize=fs, style="italic", color="#444444")

    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Panel A — Mantel test (one-tailed permutation histogram)
# ─────────────────────────────────────────────────────────────────────────────

def _draw_mantel_panel(ax, mantel_result: dict, fs: int) -> None:
    """
    One-tailed permutation histogram for the Mantel test.

    Design decisions vs the original:
    • One-tailed shading only (right tail) — we have a directional hypothesis
      (positive correlation expected), so shading the left tail is misleading.
    • X-axis extended to include observed r with 10 % breathing room so the
      red line is never clipped.
    • Kernel density estimate overlaid on histogram gives the distributional
      shape without discretisation artefacts.
    • Exact permutation fraction shown ("0 / 9999 permutations ≥ observed r")
      which is more informative than just the p-value.
    • Effect-size interpretation appended below the p annotation.
    """
    perm_rs = mantel_result.get("_perm_rs")
    obs_r   = mantel_result.get("r")
    p_val   = mantel_result.get("p_value")
    n_perm  = mantel_result.get("n_permutations", len(perm_rs) if perm_rs is not None else 0)

    if perm_rs is None or obs_r is None:
        ax.text(0.5, 0.5, "Mantel test not available",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=fs, color="#888888")
        ax.set_title("A", fontsize=fs + 1, fontweight="bold")
        _label_panel(ax, "A", fs)
        return

    perm_rs = np.asarray(perm_rs)

    # ── histogram ─────────────────────────────────────────────────────────
    x_lo = min(perm_rs.min(), -0.01)
    x_hi = max(obs_r * 1.15, perm_rs.max() * 1.05)

    bins = np.linspace(x_lo, x_hi, 55)
    ax.hist(perm_rs, bins=bins,
            color=_COL["null_hist"], edgecolor="white",
            linewidth=0.3, alpha=0.80,
            label=f"Null distribution\n({n_perm:,} permutations)")

    # ── KDE overlay ───────────────────────────────────────────────────────
    try:
        from scipy.stats import gaussian_kde
        kde   = gaussian_kde(perm_rs, bw_method=0.12)
        x_kde = np.linspace(x_lo, x_hi, 400)
        scale = len(perm_rs) * (bins[1] - bins[0])
        ax.plot(x_kde, kde(x_kde) * scale,
                color="#1A4E8C", linewidth=1.5, alpha=0.75)
    except Exception:
        pass

    # ── one-tailed critical region (right tail only) ───────────────────────
    crit = np.percentile(perm_rs, 95)      # one-tailed 95th percentile
    ax.axvspan(crit, x_hi,
               alpha=0.18, color=_COL["null_shade"],
               label=f"One-tailed p < 0.05\n(r > {crit:.3f})")

    # ── observed r line ────────────────────────────────────────────────────
    ax.axvline(obs_r, color=_COL["obs_line"], linewidth=2.2, zorder=5,
               label=f"Observed r = {obs_r:.3f}")

    # ── p-value / permutation count annotation ─────────────────────────────
    n_exceed = int((perm_rs >= obs_r).sum())
    if p_val is not None and p_val < 0.0001:
        p_str = f"p < 0.0001"
    elif p_val is not None:
        p_str = f"p = {p_val:.4f}"
    else:
        p_str = f"p = {n_exceed}/{n_perm}"

    effect_str = (
        "Very strong" if obs_r >= 0.30 else
        "Strong"      if obs_r >= 0.20 else
        "Moderate"    if obs_r >= 0.10 else
        "Weak"
    )

    ax.text(0.97, 0.97,
            f"{p_str} (one-tailed)\n"
            f"{n_exceed}/{n_perm:,} permutations ≥ r\n"
            f"Effect: {effect_str} (r = {obs_r:.3f})",
            ha="right", va="top", transform=ax.transAxes,
            fontsize=max(6, fs - 1), color=_COL["obs_line"],
            bbox=dict(boxstyle="round,pad=0.4",
                      fc="white", ec=_COL["obs_line"], alpha=0.90,
                      linewidth=1.2))

    ax.set_xlim(x_lo, x_hi)
    ax.set_xlabel("Permuted Mantel r", fontsize=fs, labelpad=4)
    ax.set_ylabel("Frequency", fontsize=fs, labelpad=4)
    ax.set_title(
        "Mantel Test\n"
        r"Correlation: ontology distance vs sequence identity",
        fontsize=fs + 1, fontweight="bold", pad=8)
    ax.legend(fontsize=max(6, fs - 2), loc="upper left",
              framealpha=0.90, edgecolor="#cccccc")
    ax.spines[_SPINE_COLS].set_visible(False)
    ax.tick_params(labelsize=fs - 1)
    _label_panel(ax, "A", fs)


# ─────────────────────────────────────────────────────────────────────────────
# Panel B — ARI / NMI with interpretation zones
# ─────────────────────────────────────────────────────────────────────────────

def _draw_ari_panel(ax, ari_result: dict, fs: int) -> None:
    """
    Horizontal bar chart for ARI and NMI with coloured interpretation zones.

    Zones (ARI / NMI, roughly):
      0.00–0.10  poor / chance
      0.10–0.30  fair
      0.30–0.60  moderate–good
      0.60–1.00  strong–perfect
    """
    if not ari_result or ari_result.get("ARI") is None:
        ax.text(0.5, 0.5, "ARI/NMI not available",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=fs, color="#888888")
        ax.set_title("B", fontsize=fs + 1, fontweight="bold")
        _label_panel(ax, "B", fs)
        return

    ari = float(ari_result["ARI"])
    nmi = float(ari_result["NMI"])

    # ── background interpretation zones ───────────────────────────────────
    zone_bounds = [(0.00, 0.10, "#F5F5F5", "chance"),
                   (0.10, 0.30, "#E8F4F8", "fair"),
                   (0.30, 0.60, "#D4EDDA", "moderate–good"),
                   (0.60, 1.00, "#C3E6CB", "strong")]
    y_lo, y_hi = -0.35, 1.35
    for x0, x1, col, lbl in zone_bounds:
        ax.axvspan(x0, x1, ymin=0, ymax=1,
                   color=col, alpha=0.55, zorder=0)
        ax.text((x0 + x1) / 2, y_hi - 0.02, lbl,
                ha="center", va="top",
                fontsize=max(5, fs - 3), color="#666666",
                style="italic")

    # ── bars ───────────────────────────────────────────────────────────────
    metrics = {"ARI\n(Adjusted Rand Index)": ari,
               "NMI\n(Normalised Mutual Info)": nmi}
    colors  = [_COL["ari"], _COL["nmi"]]
    ys      = [0, 1]
    bars    = ax.barh(ys, list(metrics.values()),
                      color=colors, height=0.42,
                      alpha=0.88, edgecolor="white",
                      linewidth=0, zorder=3)

    # ── value labels ───────────────────────────────────────────────────────
    for bar, val in zip(bars, metrics.values()):
        x_txt = val + 0.02 if val < 0.85 else val - 0.04
        ha    = "left"      if val < 0.85 else "right"
        ax.text(x_txt, bar.get_y() + bar.get_height() / 2,
                f"{val:.3f}", va="center", ha=ha,
                fontsize=fs, fontweight="bold", color="#111111", zorder=4)

    # ── reference lines ───────────────────────────────────────────────────
    ax.axvline(0.0, color="#999999", linewidth=1.0,
               linestyle="--", alpha=0.7, zorder=2)
    ax.axvline(1.0, color="#B71C1C", linewidth=1.0,
               linestyle="--", alpha=0.5, zorder=2)

    ax.set_yticks(ys)
    ax.set_yticklabels(list(metrics.keys()), fontsize=fs)
    ax.set_xlim(-0.05, 1.10)
    ax.set_ylim(y_lo, y_hi)
    ax.set_xlabel("Score", fontsize=fs, labelpad=4)
    n_ont = ari_result.get("n_ontology_clusters", "?")
    n_seq = ari_result.get("n_mmseqs2_clusters", "?")
    n_p   = ari_result.get("n_proteins", "?")
    ax.set_title(
        "Cluster Agreement\n"
        f"Ontology: {n_ont} clusters  ·  MMseqs2: {n_seq} clusters  ·  "
        f"n = {n_p} proteins",
        fontsize=fs + 1, fontweight="bold", pad=8)
    ax.spines[_SPINE_COLS].set_visible(False)
    ax.tick_params(labelsize=fs - 1, left=False)
    _label_panel(ax, "B", fs)


# ─────────────────────────────────────────────────────────────────────────────
# Panel C — PCoA (ontology distances, superfamily colours)
# ─────────────────────────────────────────────────────────────────────────────

def _draw_pcoa_panel(ax,
                     pcoa_coords: pd.DataFrame,
                     group_series: pd.Series,
                     mmseqs2_labels: pd.Series,
                     fs: int) -> None:
    """
    PCoA on ontology Jaccard distances.
    Colour = ontology superfamily (same palette as Stage 13).
    Marker shape = superfamily.
    Small grey seq-cluster numbers annotated on each point.
    """
    if pcoa_coords.empty or "PCoA1" not in pcoa_coords.columns:
        ax.text(0.5, 0.5, "PCoA not available",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=fs, color="#888888")
        ax.set_title("C", fontsize=fs + 1, fontweight="bold")
        _label_panel(ax, "C", fs)
        return

    has_seq = mmseqs2_labels is not None and len(mmseqs2_labels) > 0
    pc2 = pcoa_coords.get("PCoA2", pd.Series(0, index=pcoa_coords.index))

    # ── convex hulls per ontology group (≥3 proteins) ─────────────────────
    from scipy.spatial import ConvexHull
    from matplotlib.patches import Polygon as MplPolygon
    for grp in sorted(group_series.unique()):
        idx = group_series[group_series == grp].index
        pts_idx = pcoa_coords.index.intersection(idx)
        if len(pts_idx) < 3:
            continue
        x = pcoa_coords.loc[pts_idx, "PCoA1"].values
        y = pc2.loc[pts_idx].values
        col = _ont_colour(grp)
        try:
            hull  = ConvexHull(np.column_stack([x, y]))
            verts = np.column_stack([x, y])[hull.vertices]
            cx, cy = verts.mean(axis=0)
            vp = verts + 0.015 * (verts - np.array([cx, cy]))
            ax.add_patch(MplPolygon(vp, closed=True,
                                    facecolor=col, alpha=0.10,
                                    edgecolor=col, linewidth=1.2,
                                    linestyle="--", zorder=1))
        except Exception:
            pass

    # ── scatter points ─────────────────────────────────────────────────────
    present_groups = sorted(group_series.unique())
    legend_handles, legend_labels = [], []
    for grp in present_groups:
        idx  = group_series[group_series == grp].index
        pts  = pcoa_coords.loc[pcoa_coords.index.intersection(idx)]
        pc2_ = pc2.loc[pts.index]
        if pts.empty:
            continue
        col = _ont_colour(grp)
        sc  = ax.scatter(pts["PCoA1"], pc2_,
                         color=col, s=70, alpha=0.88,
                         edgecolors="white", linewidths=0.5,
                         zorder=3, label=grp)
        legend_handles.append(sc)
        legend_labels.append(grp)

    # ── seq cluster number annotations ─────────────────────────────────────
    if has_seq:
        for prot in pcoa_coords.index:
            sc_id = mmseqs2_labels.get(prot, None)
            if sc_id is not None:
                ax.annotate(str(sc_id),
                            xy=(pcoa_coords.loc[prot, "PCoA1"], pc2.loc[prot]),
                            xytext=(0, 4), textcoords="offset points",
                            fontsize=max(4, fs - 5), alpha=0.45,
                            color="#333333", ha="center", zorder=4)

    ax.set_xlabel("PCoA1", fontsize=fs, labelpad=4)
    ax.set_ylabel("PCoA2", fontsize=fs, labelpad=4)
    seq_note = "  ·  grey number = seq cluster" if has_seq else ""
    ax.set_title(
        f"PCoA on Ontology Jaccard Distances\n"
        f"colour = ontology group  ·  hull = group boundary{seq_note}",
        fontsize=fs + 1, fontweight="bold", pad=8)

    # Compact legend (≤16 entries, outside right)
    if len(legend_handles) <= 20:
        ax.legend(legend_handles, legend_labels,
                  fontsize=max(5, fs - 3),
                  bbox_to_anchor=(1.01, 1), loc="upper left",
                  frameon=True, framealpha=0.92, edgecolor="#cccccc",
                  handlelength=0.9, handletextpad=0.4,
                  borderpad=0.5, labelspacing=0.3,
                  markerscale=0.9)
    ax.spines[_SPINE_COLS].set_visible(False)
    ax.tick_params(labelsize=fs - 1)
    _label_panel(ax, "C", fs)


# ─────────────────────────────────────────────────────────────────────────────
# Panel D — Enrichment heatmap with odds-ratio dot overlay
# ─────────────────────────────────────────────────────────────────────────────

def _draw_enrichment_panel(ax, enrichment_df: pd.DataFrame, fs: int) -> None:
    """
    Heatmap: rows = MMseqs2 clusters with ≥1 significant hit,
             columns = top enriched features (sorted by max −log10 FDR).
    Cell background = −log10(FDR) on YlOrRd scale.
    White dot in significant cells scaled by odds ratio.
    Only clusters and features with at least one hit are shown.
    """
    if enrichment_df.empty:
        ax.text(0.5, 0.5, "No enrichment data",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=fs, color="#888888")
        ax.set_title("D", fontsize=fs + 1, fontweight="bold")
        _label_panel(ax, "D", fs)
        return

    sig = enrichment_df[enrichment_df["significant"]].copy()
    if sig.empty:
        ax.text(0.5, 0.5, "No enriched terms at FDR < 0.05",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=fs, color="#888888")
        ax.set_title("D", fontsize=fs + 1, fontweight="bold")
        _label_panel(ax, "D", fs)
        return

    sig["-log10_fdr"] = -np.log10(sig["fdr"].clip(lower=1e-300))

    # Shorten feature names
    def _shorten(c: str) -> str:
        parts = c.split("::")
        label = parts[-1] if len(parts) > 1 else c
        label = label.split(":")[-1] if ":" in label else label
        return label[:18]

    # Top features by max enrichment
    top_feats = (
        sig.groupby("feature")["-log10_fdr"].max()
           .sort_values(ascending=False)
           .head(25).index.tolist()
    )
    # Only clusters that have at least one hit in top_feats
    clusters = sorted(
        sig[sig["feature"].isin(top_feats)]["cluster_id"].unique(),
        key=lambda x: int(x) if str(x).isdigit() else 0,
    )

    # Build −log10 FDR matrix and odds ratio matrix
    short_feats = [_shorten(f) for f in top_feats]
    heat = pd.DataFrame(0.0, index=clusters, columns=top_feats)
    odds = pd.DataFrame(0.0, index=clusters, columns=top_feats)

    for _, row in sig[sig["feature"].isin(top_feats) &
                      sig["cluster_id"].isin(clusters)].iterrows():
        heat.loc[row["cluster_id"], row["feature"]] = row["-log10_fdr"]
        odds.loc[row["cluster_id"], row["feature"]] = float(
            row.get("odds_ratio", 0) or 0
        )

    heat.columns = short_feats
    odds.columns = short_feats

    # ── heatmap ───────────────────────────────────────────────────────────
    vmax = heat.values.max()
    im   = ax.imshow(heat.values, aspect="auto",
                     cmap="YlOrRd", vmin=0, vmax=vmax,
                     interpolation="nearest")

    cbar = plt.colorbar(im, ax=ax, shrink=0.65, pad=0.02)
    cbar.set_label("−log₁₀(FDR)", fontsize=max(6, fs - 1))
    cbar.ax.tick_params(labelsize=max(5, fs - 2))

    # ── odds-ratio dot overlay ─────────────────────────────────────────────
    or_max = max(odds.values.max(), 1.0)
    for i, cl in enumerate(clusters):
        for j, ft in enumerate(short_feats):
            v_fdr = heat.iloc[i, j]
            v_or  = odds.iloc[i, j]
            if v_fdr > 0:
                # dot size proportional to log(OR)
                dot_r = 100 * (np.log1p(v_or) / np.log1p(or_max + 1))
                ax.scatter(j, i, s=max(8, dot_r),
                           color="white", alpha=0.85,
                           edgecolors="#333333", linewidths=0.4,
                           zorder=4)

    # ── axes ─────────────────────────────────────────────────────────────
    ax.set_xticks(range(len(short_feats)))
    ax.set_xticklabels(short_feats, rotation=55, ha="right",
                       fontsize=max(5, fs - 2))
    ax.set_yticks(range(len(clusters)))
    ax.set_yticklabels([f"Seq {c}" for c in clusters],
                       fontsize=max(6, fs - 1))

    # Gridlines between cells
    for x in np.arange(-0.5, len(short_feats), 1):
        ax.axvline(x, color="white", linewidth=0.5, zorder=2)
    for y in np.arange(-0.5, len(clusters), 1):
        ax.axhline(y, color="white", linewidth=0.5, zorder=2)

    # Odds-ratio legend (dot size)
    from matplotlib.lines import Line2D
    dot_legend = [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor="white", markeredgecolor="#333333",
               markeredgewidth=0.5,
               markersize=np.sqrt(max(8, 100 * np.log1p(v) /
                                      np.log1p(or_max + 1))),
               label=f"OR = {v:.0f}")
        for v in [2, 5, 10, 50] if v <= or_max + 1
    ]
    if dot_legend:
        leg = ax.legend(handles=dot_legend,
                        title="Odds ratio", title_fontsize=max(5, fs - 2),
                        fontsize=max(5, fs - 2),
                        loc="lower right",
                        framealpha=0.88, edgecolor="#cccccc",
                        handlelength=0.5, borderpad=0.5,
                        labelspacing=0.4)

    ax.set_title(
        "Ontology Term Enrichment per Sequence Cluster\n"
        f"Fisher's exact, BH-FDR < 0.05  ·  top {len(top_feats)} features  ·  "
        "dot size = odds ratio",
        fontsize=fs + 1, fontweight="bold", pad=8)
    ax.tick_params(labelsize=fs - 1, bottom=False)
    _label_panel(ax, "D", fs)

# 9 — Stage runner
# ══════════════════════════════════════════════════════════════════════════════

def run(cfg: dict) -> None:
    logger = setup_logging(cfg.get("log_dir", "logs"),
                           "stage16_sequence_ontology_validation")
    stage_banner(logger, "STAGE 16",
                 "Sequence vs Ontology Congruence Analysis")
    out_dir  = make_output_dir(cfg, "stage16_congruence")
    plot_cfg = cfg.get("plotting", {})
    dpi      = plot_cfg.get("figure_dpi", 150)
    fs       = plot_cfg.get("font_size", 9)
    n_perms  = cfg.get("mantel_permutations", 9999)

    # ── Load occurrence table ────────────────────────────────────────────
    occ_path = get_input_path(cfg, "domain_occurrence_table")
    if not occ_path:
        logger.error("domain_occurrence_table not found — aborting Stage 16")
        return
    occ = load_table(occ_path)
    occ = rename_to_canonical(occ, cfg.get("column_aliases", {}))

    # ── Build feature matrix ─────────────────────────────────────────────
    cluster_cfg   = cfg.get("clustering", {})
    feature_types = cluster_cfg.get(
        "feature_types",
        ["domain_types", "function_tags", "topology_context", "binding_partner"]
    )
    feat_col_map = {
        "domain_types":        "domain_type_id",
        "function_tags":       "function_tags",
        "topology_context":    "topology_context",
        "binding_partner":     next(
            (c for c in ["binding_partner", "binding_parnter"] if c in occ.columns),
            "binding_partner",
        ),
        "proximity_categories": "proximity_categories",
    }
    feat_cols = [feat_col_map[k] for k in feature_types
                 if feat_col_map.get(k) in occ.columns]
    logger.info(f"Feature columns: {feat_cols}")

    feat_mat = build_feature_matrix(occ, feat_cols)
    logger.info(f"Feature matrix: {feat_mat.shape}  (proteins × features)")
    save_tsv(
        feat_mat.reset_index().rename(columns={"index": "protein_accession"}),
        out_dir / "feature_matrix.tsv", logger
    )

    # ── Ontology distance matrix ─────────────────────────────────────────
    ont_dist = ontology_distance_matrix(feat_mat)
    save_tsv(
        ont_dist.reset_index().rename(columns={"index": "protein_accession"}),
        out_dir / "distance_matrix_ontology.tsv", logger
    )

    # ── Load protein group labels ────────────────────────────────────────
    group_series = pd.Series("other", index=feat_mat.index, name="group")
    grp_path = get_input_path(cfg, "protein_group_labels_sequening")
    if grp_path:
        grp_df = load_table(grp_path)
        acc_c = next((c for c in ["protein_accession", "accession"]
                      if c in grp_df.columns), None)
        grp_c = next((c for c in ["protein_group", "group", "protein_family"]
                      if c in grp_df.columns), None)
        if acc_c and grp_c:
            for _, row in grp_df.iterrows():
                acc = str(row[acc_c]).strip()
                if acc in group_series.index:
                    group_series[acc] = str(row[grp_c]).strip()
            logger.info(f"Loaded {grp_df[grp_c].nunique()} biological groups")

    # ── PCoA ─────────────────────────────────────────────────────────────
    logger.info("Running PCoA on ontology distance matrix …")
    pcoa_coords = run_pcoa(ont_dist)
    if not pcoa_coords.empty:
        stress = pcoa_coords.attrs.get("stress", None)
        if stress is not None:
            logger.info(f"  PCoA stress: {stress:.4f}"
                        + ("  (good)" if stress < 0.1 else
                           "  (acceptable)" if stress < 0.2 else
                           "  (high — interpret with caution)"))
        save_tsv(
            pcoa_coords.reset_index().rename(
                columns={"index": "protein_accession"}),
            out_dir / "pcoa_coordinates.tsv", logger
        )

    # ── Load MMseqs2 clusters (optional) ────────────────────────────────
    mmseqs2_df     = None
    mmseqs2_labels = pd.Series(dtype=str)
    seq_dist       = pd.DataFrame()

    mmseqs2_path = get_input_path(cfg, "mmseqs2_clusters")
    search_path  = get_input_path(cfg, "mmseqs2_search_results")  # pairwise_identity.tsv

    if mmseqs2_path and Path(mmseqs2_path).exists():
        logger.info(f"Loading MMseqs2 clusters from {mmseqs2_path} …")
        mmseqs2_df = load_table(mmseqs2_path)
        acc_c = next((c for c in ["accession", "protein_accession"]
                      if c in mmseqs2_df.columns), None)
        cid_c = next((c for c in ["cluster_id", "cluster"]
                      if c in mmseqs2_df.columns), None)
        if acc_c and cid_c:
            mmseqs2_labels = pd.Series(
                mmseqs2_df.set_index(acc_c)[cid_c].astype(str)
            )
            logger.info(f"  {mmseqs2_labels.nunique()} MMseqs2 clusters, "
                        f"{len(mmseqs2_labels)} proteins")

            # Use continuous pairwise identity if available (richer Mantel test)
            # otherwise fall back to binary cluster-membership distance
            if search_path and Path(search_path).exists():
                logger.info(f"  Using pairwise identity matrix from {search_path}")
                logger.info("  (continuous distance — more sensitive Mantel test)")
                seq_dist = sequence_distance_matrix_from_identity(
                    search_path, feat_mat.index.tolist()
                )
            else:
                logger.info("  No pairwise_identity.tsv found — using binary "
                            "cluster-membership distance proxy")
                logger.info("  To upgrade: set mmseqs2_search_results: "
                            "input/pairwise_identity.tsv in config")
                seq_dist = sequence_distance_matrix(
                    mmseqs2_df, feat_mat.index.tolist(), acc_c, cid_c
                )

            save_tsv(
                seq_dist.reset_index().rename(
                    columns={"index": "protein_accession"}),
                out_dir / "distance_matrix_sequence.tsv", logger
            )
        else:
            logger.warning("MMseqs2 file missing 'accession' / 'cluster_id' columns")
    else:
        logger.info(
            "No MMseqs2 clusters file provided.\n"
            "  → Mantel test, ARI/NMI, and enrichment analysis will be skipped.\n"
            "  → To enable:\n"
            "      1. python fetch_sequences.py              (get FASTA)\n"
            "      2. sbatch mmseqs2_cluster.sh              (run on HPC)\n"
            "      3. Copy results/*.tsv to input/\n"
            "      4. Add to pipeline_config.yaml:\n"
            "           mmseqs2_clusters:       input/camp_family_cluster_ids.tsv\n"
            "           mmseqs2_search_results: input/pairwise_identity.tsv\n"
            "           protein_group_labels:   input/protein_group_labels_sequencing.tsv"
        )

    # ── Mantel test ──────────────────────────────────────────────────────
    mantel_result = {"r": None, "p_value": None, "_perm_rs": None}
    if not seq_dist.empty:
        logger.info(f"Running Mantel test ({n_perms} permutations) …")
        mantel_result = mantel_test(ont_dist, seq_dist, n_perms,
                                    method="spearman", logger=logger)
        logger.info(f"  r = {mantel_result['r']:.4f}, "
                    f"p = {mantel_result['p_value']:.4f}")
        logger.info(f"  {mantel_result['interpretation']}")

        mantel_save = {k: v for k, v in mantel_result.items()
                       if k != "_perm_rs"}
        save_tsv(pd.DataFrame([mantel_save]),
                 out_dir / "mantel_test_result.tsv", logger)

    # ── ARI / NMI ────────────────────────────────────────────────────────
    ari_result = {"ARI": None, "NMI": None}
    if not mmseqs2_labels.empty:
        # Use hierarchical clustering on ontology distances as ontology labels
        from scipy.cluster.hierarchy import linkage, fcluster
        from scipy.spatial.distance import squareform as sf
        common_prots = sorted(
            set(feat_mat.index) & set(mmseqs2_labels.index)
        )
        n_ont_clusters = mmseqs2_labels.loc[common_prots].nunique()
        dist_vec = sf(
            ont_dist.loc[common_prots, common_prots].values.astype(float),
            checks=False
        )
        link = linkage(dist_vec, method="ward")
        ont_clust_labels = pd.Series(
            fcluster(link, t=n_ont_clusters, criterion="maxclust"),
            index=common_prots,
            name="ontology_cluster",
        ).astype(str)

        logger.info(f"Computing ARI/NMI  "
                    f"(ontology: {n_ont_clusters} clusters, "
                    f"MMseqs2: {mmseqs2_labels.nunique()} clusters) …")
        ari_result = cluster_agreement(ont_clust_labels, mmseqs2_labels)
        logger.info(f"  ARI = {ari_result['ARI']:.4f}  "
                    f"NMI = {ari_result['NMI']:.4f}  "
                    f"({ari_result['ARI_interpretation']})")
        save_tsv(pd.DataFrame([ari_result]),
                 out_dir / "ari_nmi_scores.tsv", logger)

    # ── Enrichment analysis ──────────────────────────────────────────────
    enrichment_df = pd.DataFrame()
    if not mmseqs2_labels.empty:
        logger.info("Running enrichment analysis …")
        enrichment_df = enrichment_analysis(feat_mat, mmseqs2_labels)
        n_sig = enrichment_df["significant"].sum() if not enrichment_df.empty else 0
        logger.info(f"  {n_sig} significant enrichments  "
                    f"(FDR < 0.05) across "
                    f"{enrichment_df['cluster_id'].nunique() if not enrichment_df.empty else 0} clusters")
        save_tsv(enrichment_df, out_dir / "enrichment_results.tsv", logger)

    # ── Validation summary ────────────────────────────────────────────────
    summary = {
        "n_proteins":            len(feat_mat),
        "n_features":            feat_mat.shape[1],
        "n_ontology_clusters":   ari_result.get("n_ontology_clusters", "N/A"),
        "n_mmseqs2_clusters":    ari_result.get("n_mmseqs2_clusters", "N/A"),
        "mantel_r":              mantel_result.get("r", "N/A"),
        "mantel_p":              mantel_result.get("p_value", "N/A"),
        "mantel_n_perms":        mantel_result.get("n_permutations", "N/A"),
        "ARI":                   ari_result.get("ARI", "N/A"),
        "NMI":                   ari_result.get("NMI", "N/A"),
        "ARI_interpretation":    ari_result.get("ARI_interpretation", "N/A"),
        "n_enriched_terms":      int(enrichment_df["significant"].sum())
                                 if not enrichment_df.empty else "N/A",
        "pcoa_stress":           pcoa_coords.attrs.get("stress", "N/A")
                                 if not pcoa_coords.empty else "N/A",
        "mantel_interpretation": mantel_result.get("interpretation", "N/A"),
    }
    save_tsv(pd.DataFrame([summary]),
             out_dir / "validation_summary.tsv", logger)

    # ── 4-panel figure ────────────────────────────────────────────────────
    logger.info("Generating congruence plot …")
    _congruence_plot(
        mantel_result, ari_result, enrichment_df,
        pcoa_coords, group_series, mmseqs2_labels,
        out_dir / "congruence_plot.png", dpi=dpi, font_size=fs,
    )
    logger.info(f"  Saved → {out_dir}/congruence_plot.png")

    # ── Final interpretation log ──────────────────────────────────────────
    logger.info("")
    logger.info("=" * 66)
    logger.info("  STAGE 16 RESULTS SUMMARY")
    logger.info("=" * 66)
    if mantel_result.get("r") is not None:
        logger.info(f"  Mantel r    = {mantel_result['r']:.4f}  "
                    f"(p = {mantel_result['p_value']:.4f})")
    if ari_result.get("ARI") is not None:
        logger.info(f"  ARI         = {ari_result['ARI']:.4f}  "
                    f"NMI = {ari_result['NMI']:.4f}")
    if not enrichment_df.empty:
        n_sig = int(enrichment_df["significant"].sum())
        logger.info(f"  Enrichment  = {n_sig} significant term-cluster associations")
    logger.info("")
    logger.info("  Interpretation:")
    logger.info(f"  {mantel_result.get('interpretation', 'Run with MMseqs2 data for full analysis')}")
    logger.info("=" * 66)
    logger.info("Stage 16 complete.")


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import sys, os
    _project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(_project_root)
    _cfg_path = sys.argv[1] if len(sys.argv) > 1 else "config/pipeline_config.yaml"
    run(load_config(_cfg_path))
