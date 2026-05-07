"""
08_evaluate_transfer.py
────────────────────────────────────────────────────────────────────────────
Step 8 — Scientific evaluation of transfer predictions against gold labels.

Evaluates on the BLIND TEST SET (ortholog_pairs_test.tsv).
Calibration set was used only in steps 06–07.

Scientific improvements over original
──────────────────────────────────────
  1. Blind test set: evaluates only on held-out 20% pairs
  2. Bootstrap 95% CI: 1 000 bootstrap resamples per method × category
  3. Permutation test: 1 000 permutations to compute empirical p-value
     — establishes statistical significance vs random pair assignment
  4. Difficulty-stratified evaluation:
       Easy   gene_similarity == 1.0  (exact gene symbol match)
       Medium gene_similarity >= 0.7  (closely related genes)
       Hard   gene_similarity <  0.7  (low name similarity)
  5. Fixed is_transferable filter in random control (was invalid pandas .get)
  6. Four publication-quality figures (300 DPI, Nature style)

Outputs
───────
  evaluation/evaluation_per_protein.tsv      per-protein per-category scores
  evaluation/evaluation_summary.tsv          mean P/R/F1/J per method × cat
  evaluation/evaluation_bootstrap_ci.tsv     bootstrap 95% CI per method × cat
  evaluation/evaluation_difficulty.tsv       scores per difficulty tier
  evaluation/evaluation_permutation.tsv      permutation p-values
  evaluation/evaluation_controls.tsv         random-pair control results

  evaluation/fig1_main_results.png           F1 per method × category + CI
  evaluation/fig2_difficulty_stratified.png  F1 by difficulty tier
  evaluation/fig3_calibration.png            confidence tier vs actual precision
  evaluation/fig4_permutation.png            null distribution vs observed F1
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations
import argparse
import random
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd


# ── Global figure style ───────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":       "DejaVu Sans",
    "font.size":         8,
    "axes.titlesize":    9,
    "axes.labelsize":    8,
    "xtick.labelsize":   7,
    "ytick.labelsize":   7,
    "legend.fontsize":   7,
    "axes.linewidth":    0.7,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        300,
    "savefig.dpi":       300,
    "savefig.bbox":      "tight",
    "figure.facecolor":  "white",
    "axes.facecolor":    "white",
    "legend.frameon":    False,
    "pdf.fonttype":      42,
})

PAL = {
    "A":       "#0077BB",
    "B":       "#EE7733",
    "C":       "#009988",
    "random":  "#BBBBBB",
    "easy":    "#44AA99",
    "medium":  "#DDAA33",
    "hard":    "#BB5566",
}

METHOD_LABELS = {
    "A_direct_transfer":     "Method A\n(direct)",
    "B_domain_architecture": "Method B\n(domain arch.)",
    "C_nearest_sequence":    "Method C\n(nearest seq.)",
    "control_random":        "Random\ncontrol",
}

DIFFICULTY_ORDER = ["easy", "medium", "hard"]


# ══════════════════════════════════════════════════════════════════════════════
# Metrics
# ══════════════════════════════════════════════════════════════════════════════

def _prf(pred: set, gold: set) -> dict:
    if not pred and not gold:
        return {"precision":1.0,"recall":1.0,"f1":1.0,
                "jaccard":1.0,"n_pred":0,"n_gold":0,"n_tp":0}
    tp = len(pred & gold)
    p  = tp / len(pred) if pred else 0.0
    r  = tp / len(gold) if gold else 0.0
    f1 = 2*p*r/(p+r) if (p+r) else 0.0
    j  = len(pred & gold) / len(pred | gold) if (pred | gold) else 0.0
    return {"precision":round(p,4),"recall":round(r,4),"f1":round(f1,4),
            "jaccard":round(j,4),"n_pred":len(pred),"n_gold":len(gold),"n_tp":tp}


def _build_gold_sets(gold_df: pd.DataFrame) -> dict[str, dict[str, set]]:
    result: dict[str, dict[str, set]] = defaultdict(lambda: defaultdict(set))
    for _, row in gold_df.iterrows():
        acc = str(row.get("protein_accession","")).strip()
        cat = str(row.get("annotation_category","")).strip()
        val = str(row.get("term_value","")).strip()
        if acc and cat and val:
            result[acc][cat].add(val)
    return result


def _build_pred_sets(pred_df: pd.DataFrame,
                      method: str) -> dict[str, dict[str, set]]:
    sub    = pred_df[pred_df["method"] == method]
    result: dict[str, dict[str, set]] = defaultdict(lambda: defaultdict(set))
    for _, row in sub.iterrows():
        acc = str(row.get("mouse_accession","")).strip()
        cat = str(row.get("annotation_category","")).strip()
        val = str(row.get("term_value","")).strip()
        if acc and cat and val:
            result[acc][cat].add(val)
    return result


# ══════════════════════════════════════════════════════════════════════════════
# Bootstrap confidence intervals
# ══════════════════════════════════════════════════════════════════════════════

def bootstrap_ci(values: list[float],
                  n: int = 1000,
                  ci: float = 0.95,
                  rng: np.random.Generator | None = None) -> tuple[float,float,float]:
    """Return (mean, lower_ci, upper_ci) via bootstrap over protein-level scores."""
    if not values:
        return 0.0, 0.0, 0.0
    if rng is None:
        rng = np.random.default_rng(42)
    arr   = np.array(values, dtype=float)
    mean_ = float(np.mean(arr))
    boots = [float(np.mean(rng.choice(arr, size=len(arr), replace=True)))
             for _ in range(n)]
    alpha = (1 - ci) / 2
    lo    = float(np.percentile(boots, 100 * alpha))
    hi    = float(np.percentile(boots, 100 * (1 - alpha)))
    return mean_, lo, hi


# ══════════════════════════════════════════════════════════════════════════════
# Evaluation core
# ══════════════════════════════════════════════════════════════════════════════

def evaluate_method(method:     str,
                     pred_df:    pd.DataFrame,
                     gold_sets:  dict,
                     categories: list[str]) -> list[dict]:
    pred_sets = _build_pred_sets(pred_df, method)
    rows = []
    for m_acc, cat_gold in gold_sets.items():
        for cat in categories:
            gold_set = cat_gold.get(cat, set())
            pred_set = pred_sets.get(m_acc, {}).get(cat, set())
            metrics  = _prf(pred_set, gold_set)
            rows.append({
                "mouse_accession": m_acc,
                "method":          method,
                "category":        cat,
                **metrics,
            })
    return rows


def random_control(gold_sets:    dict,
                    human_ann_df: pd.DataFrame,
                    categories:   list[str],
                    n_random:     int = 100,
                    seed:         int = 42) -> list[dict]:
    """Random-pair control: predict from a random human protein (not the ortholog)."""
    if human_ann_df.empty:
        return []

    # Fix: correct pandas filtering for is_transferable column
    if "is_transferable" in human_ann_df.columns:
        mask     = human_ann_df["is_transferable"].astype(str).str.lower().isin(["true","yes","1"])
        h_ann_df = human_ann_df[mask]
    else:
        h_ann_df = human_ann_df

    human_accs = h_ann_df["protein_accession"].unique().tolist()
    if not human_accs:
        return []

    h_sets: dict[str, dict[str, set]] = defaultdict(lambda: defaultdict(set))
    for _, row in h_ann_df.iterrows():
        acc = str(row.get("protein_accession","")).strip()
        cat = str(row.get("annotation_category","")).strip()
        val = str(row.get("term_value","")).strip()
        if acc and cat and val:
            h_sets[acc][cat].add(val)

    rng  = random.Random(seed)
    rows = []
    for m_acc, cat_gold in gold_sets.items():
        for cat in categories:
            gold_set    = cat_gold.get(cat, set())
            jaccards    = []
            for _ in range(n_random):
                rand_acc = rng.choice(human_accs)
                pred_set = h_sets.get(rand_acc, {}).get(cat, set())
                jaccards.append(_prf(pred_set, gold_set)["jaccard"])
            rows.append({
                "mouse_accession":  m_acc,
                "method":           "control_random",
                "category":         cat,
                "jaccard":          round(float(np.mean(jaccards)), 4),
                "jaccard_std":      round(float(np.std(jaccards)), 4),
                "precision":        0.0,
                "recall":           0.0,
                "f1":               float(np.mean(jaccards)),  # use jaccard as F1 proxy
            })
    return rows


# ══════════════════════════════════════════════════════════════════════════════
# Permutation test
# ══════════════════════════════════════════════════════════════════════════════

def permutation_test(method:     str,
                      pairs:      pd.DataFrame,
                      human_ann:  dict[str, pd.DataFrame],
                      gold_sets:  dict,
                      categories: list[str],
                      n_perm:     int = 1000,
                      seed:       int = 42) -> dict[str, dict]:
    """
    Permutation test: shuffle mouse_accession in the pairs table.
    For each permutation, compute mean F1 per category.
    p-value = fraction of permuted F1 ≥ observed F1.
    """
    from importlib import import_module
    method_07 = import_module("07_generate_transfer_predictions")

    rng = np.random.default_rng(seed)
    m_accs = pairs["mouse_accession"].values.copy()

    # Observed F1
    real_preds = _build_pred_sets_from_method(method, pairs, human_ann)
    real_f1    = _mean_f1_per_cat(real_preds, gold_sets, categories)

    # Permuted F1 distributions
    perm_f1: dict[str, list[float]] = {cat: [] for cat in categories}

    for _ in range(n_perm):
        shuffled = pairs.copy()
        shuffled["mouse_accession"] = rng.permutation(m_accs)
        perm_preds = _build_pred_sets_from_method(method, shuffled, human_ann)
        pf1        = _mean_f1_per_cat(perm_preds, gold_sets, categories)
        for cat in categories:
            perm_f1[cat].append(pf1.get(cat, 0.0))

    results: dict[str, dict] = {}
    for cat in categories:
        obs    = real_f1.get(cat, 0.0)
        null   = np.array(perm_f1[cat])
        p_val  = float(np.mean(null >= obs))
        results[cat] = {
            "observed_f1":     round(obs, 4),
            "null_mean":       round(float(np.mean(null)), 4),
            "null_std":        round(float(np.std(null)), 4),
            "p_value":         round(p_val, 4),
            "null_distribution": null.tolist(),
        }
    return results


def _build_pred_sets_from_method(method: str,
                                   pairs:  pd.DataFrame,
                                   human_ann: dict[str, pd.DataFrame]
                                   ) -> dict[str, dict[str, set]]:
    """Quick direct-transfer predictions for permutation test."""
    preds: dict[str, dict[str, set]] = defaultdict(lambda: defaultdict(set))
    for _, pair in pairs.iterrows():
        h_acc  = str(pair.get("human_accession",""))
        m_acc  = str(pair.get("mouse_accession",""))
        anns   = human_ann.get(h_acc)
        if anns is None or anns.empty:
            continue
        for _, ann in anns.iterrows():
            cat = str(ann.get("annotation_category",""))
            val = str(ann.get("term_value",""))
            if cat and val:
                preds[m_acc][cat].add(val)
    return preds


def _mean_f1_per_cat(preds:      dict[str, dict[str, set]],
                      gold_sets:  dict,
                      categories: list[str]) -> dict[str, float]:
    f1_per_cat: dict[str, list[float]] = {cat: [] for cat in categories}
    for m_acc, cat_gold in gold_sets.items():
        for cat in categories:
            gold = cat_gold.get(cat, set())
            pred = preds.get(m_acc, {}).get(cat, set())
            f1_per_cat[cat].append(_prf(pred, gold)["f1"])
    return {cat: float(np.mean(vals)) if vals else 0.0
            for cat, vals in f1_per_cat.items()}


# ══════════════════════════════════════════════════════════════════════════════
# Difficulty-stratified evaluation
# ══════════════════════════════════════════════════════════════════════════════

def assign_difficulty(pairs: pd.DataFrame) -> pd.DataFrame:
    """Assign difficulty tier based on gene symbol similarity."""
    def _tier(row):
        gs = float(row.get("gene_similarity", 0.0))
        if str(row.get("gene_symbol_exact","")).lower() == "true" or gs >= 0.999:
            return "easy"
        elif gs >= 0.7:
            return "medium"
        else:
            return "hard"
    pairs = pairs.copy()
    pairs["difficulty"] = pairs.apply(_tier, axis=1)
    return pairs


def evaluate_by_difficulty(method:     str,
                             pred_df:    pd.DataFrame,
                             gold_sets:  dict,
                             pairs:      pd.DataFrame,
                             categories: list[str]) -> list[dict]:
    rows = []
    for tier in DIFFICULTY_ORDER:
        tier_pairs = pairs[pairs["difficulty"] == tier]
        tier_accs  = set(tier_pairs["mouse_accession"].astype(str))
        tier_gold  = {acc: cats for acc, cats in gold_sets.items()
                      if acc in tier_accs}
        if not tier_gold:
            continue
        for row in evaluate_method(method, pred_df, tier_gold, categories):
            row["difficulty"] = tier
            rows.append(row)
    return rows


# ══════════════════════════════════════════════════════════════════════════════
# Calibration: confidence tier → actual precision
# ══════════════════════════════════════════════════════════════════════════════

def compute_calibration(pred_df:   pd.DataFrame,
                          gold_sets: dict,
                          methods:   list[str]) -> pd.DataFrame:
    conf_order = ["high","medium","medium_proxy","sequence_based","low"]
    rows = []
    for method in methods:
        sub = pred_df[pred_df["method"] == method].copy()
        if sub.empty or "transfer_confidence" not in sub.columns:
            continue
        for conf in conf_order:
            csub = sub[sub["transfer_confidence"] == conf]
            if csub.empty:
                continue
            tp = fp = 0
            for _, row in csub.iterrows():
                acc  = str(row.get("mouse_accession",""))
                cat  = str(row.get("annotation_category",""))
                val  = str(row.get("term_value",""))
                gold = gold_sets.get(acc, {}).get(cat, set())
                if val in gold:
                    tp += 1
                else:
                    fp += 1
            precision = tp / (tp + fp) if (tp + fp) else 0.0
            rows.append({
                "method":    method,
                "confidence":conf,
                "n_preds":   tp + fp,
                "n_tp":      tp,
                "precision": round(precision, 4),
            })
    return pd.DataFrame(rows)


# ══════════════════════════════════════════════════════════════════════════════
# Figures (300 DPI, Nature style)
# ══════════════════════════════════════════════════════════════════════════════

def _method_color(method: str) -> str:
    for k, v in PAL.items():
        if k in method:
            return v
    return PAL["random"]


def plot_fig1_main_results(summary_df: pd.DataFrame,
                            ci_df:      pd.DataFrame,
                            out_path:   Path) -> None:
    """
    Fig 1 — Mean F1 per method × category with 95% bootstrap CI error bars.
    Methods are grouped within each category for direct comparison.
    """
    metrics   = ["recall","precision","f1"]
    methods   = [m for m in summary_df["method"].unique()
                 if m != "control_random"]
    if not methods:
        return
    categories = sorted(summary_df["category"].unique())
    n_cat      = len(categories)
    n_met      = len(methods)

    fig, axes = plt.subplots(1, 3, figsize=(4.5*n_cat, 5.5),
                              constrained_layout=True)
    if n_cat == 1:
        axes = [axes]

    x = np.arange(n_cat)
    w = min(0.75 / n_met, 0.22)

    for ax, metric in zip(axes, metrics):
        for i, method in enumerate(methods):
            vals, lo_err, hi_err = [], [], []
            for cat in categories:
                r = summary_df[(summary_df["method"]==method) &
                               (summary_df["category"]==cat)]
                v = float(r[metric].iloc[0]) if not r.empty else 0.0
                vals.append(v)

                ci_r = ci_df[(ci_df["method"]==method) &
                              (ci_df["category"]==cat) &
                              (ci_df["metric"]==metric)]
                if not ci_r.empty:
                    lo_err.append(v - float(ci_r["ci_lower"].iloc[0]))
                    hi_err.append(float(ci_r["ci_upper"].iloc[0]) - v)
                else:
                    lo_err.append(0); hi_err.append(0)

            offset = (i - n_met/2 + 0.5) * w
            bars   = ax.bar(x + offset, vals, w,
                            label=METHOD_LABELS.get(method, method),
                            color=_method_color(method),
                            edgecolor="white", linewidth=0.4, alpha=0.88)
            ax.errorbar(x + offset, vals,
                        yerr=[lo_err, hi_err],
                        fmt="none", color="#333333",
                        capsize=2.5, capthick=0.8, elinewidth=0.8, zorder=5)
            for bar, v in zip(bars, vals):
                if v > 0.04:
                    ax.text(bar.get_x() + bar.get_width()/2,
                            bar.get_height() + max(hi_err) * 0.05,
                            f"{v:.2f}", ha="center", va="bottom",
                            fontsize=6, rotation=0)

        ax.set_xlim(-0.6, n_cat - 0.4)
        ax.set_ylim(0, 1.12)
        ax.set_xticks(x)
        ax.set_xticklabels([c.replace("_"," ") for c in categories],
                            rotation=30, ha="right")
        ax.set_ylabel(metric.upper())
        ax.set_title(metric.upper(), fontweight="bold", pad=6)
        ax.axhline(0.5, color="#cccccc", lw=0.7, ls="--", alpha=0.7)
        if ax == axes[0]:
            ax.legend(fontsize=6.5, handlelength=1.2, loc="upper right")

    fig.suptitle("Transfer Evaluation: Method Comparison per Annotation Category\n"
                 "Error bars = 95% bootstrap confidence interval (n=1000)",
                 fontsize=9, fontweight="bold")
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  ✓ Fig 1 → {out_path}")


def plot_fig2_difficulty_stratified(diff_df: pd.DataFrame,
                                     out_path: Path) -> None:
    """
    Fig 2 — F1 per difficulty tier (easy / medium / hard) per method.
    Shows whether methods work beyond trivially easy orthologs.
    """
    if diff_df.empty:
        print("  ⚠  No difficulty-stratified data — skipping Fig 2")
        return

    methods    = [m for m in diff_df["method"].unique()
                  if m != "control_random"]
    categories = sorted(diff_df["category"].unique())
    n_cat      = len(categories)
    n_met      = len(methods)

    fig, axes = plt.subplots(1, len(DIFFICULTY_ORDER),
                              figsize=(5.5*len(DIFFICULTY_ORDER), 4.5),
                              constrained_layout=True)

    x = np.arange(n_cat); w = min(0.7/n_met, 0.22)

    for ax, tier in zip(axes, DIFFICULTY_ORDER):
        sub = diff_df[diff_df["difficulty"] == tier]
        if sub.empty:
            ax.text(0.5, 0.5, f"No {tier} pairs",
                    ha="center", va="center", transform=ax.transAxes, color="#999")
            ax.set_title(tier.capitalize(), fontweight="bold")
            continue
        for i, method in enumerate(methods):
            vals = []
            for cat in categories:
                r = sub[(sub["method"]==method)&(sub["category"]==cat)]
                vals.append(float(r["f1"].mean()) if not r.empty else 0.0)
            offset = (i - n_met/2 + 0.5) * w
            bars = ax.bar(x + offset, vals, w,
                          label=METHOD_LABELS.get(method, method),
                          color=_method_color(method),
                          edgecolor="white", lw=0.4, alpha=0.88)
            for bar, v in zip(bars, vals):
                if v > 0.04:
                    ax.text(bar.get_x()+bar.get_width()/2, v+0.02,
                            f"{v:.2f}", ha="center", fontsize=6)
        ax.set_xlim(-0.6, n_cat-0.4); ax.set_ylim(0, 1.1)
        ax.set_xticks(x)
        ax.set_xticklabels([c.replace("_"," ") for c in categories],
                            rotation=30, ha="right")
        ax.set_ylabel("F1")
        ax.set_title(f"{tier.capitalize()} pairs\n"
                     f"(gene sim. {'=1.0' if tier=='easy' else '≥0.7' if tier=='medium' else '<0.7'})",
                     fontweight="bold", pad=6)
        ax.axhline(0.5, color="#cccccc", lw=0.7, ls="--", alpha=0.7)
        if ax == axes[0]:
            ax.legend(fontsize=6.5, handlelength=1.2)

    fig.suptitle("Difficulty-Stratified Evaluation  (F1 score)\n"
                 "Easy = exact gene match · Medium = related name · Hard = dissimilar name",
                 fontsize=9, fontweight="bold")
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  ✓ Fig 2 → {out_path}")


def plot_fig3_calibration(calib_df: pd.DataFrame, out_path: Path) -> None:
    """
    Fig 3 — Calibration curve: transfer confidence tier vs actual precision.
    A well-calibrated method has monotone increasing precision with confidence.
    """
    if calib_df.empty:
        print("  ⚠  No calibration data — skipping Fig 3")
        return

    conf_order = ["high","medium","medium_proxy","sequence_based"]
    methods    = calib_df["method"].unique()
    fig, ax    = plt.subplots(figsize=(6, 4.5), constrained_layout=True)

    x = np.arange(len(conf_order))
    for method in methods:
        sub = calib_df[calib_df["method"] == method]
        ys  = []
        ns  = []
        for conf in conf_order:
            r = sub[sub["confidence"] == conf]
            ys.append(float(r["precision"].iloc[0]) if not r.empty else np.nan)
            ns.append(int(r["n_preds"].iloc[0]) if not r.empty else 0)
        col  = _method_color(method)
        mask = ~np.isnan(ys)
        ax.plot(x[mask], np.array(ys)[mask], "o-",
                label=METHOD_LABELS.get(method, method),
                color=col, lw=1.5, ms=6,
                markeredgecolor="white", markeredgewidth=0.6, zorder=3)
        for xi, (y, n) in zip(x, zip(ys, ns)):
            if not np.isnan(y) and n > 0:
                ax.text(xi, y + 0.025, f"n={n}", ha="center",
                        fontsize=5.5, color=col)

    ax.axhline(1.0, color="#cccccc", lw=0.7, ls=":", alpha=0.6)
    ax.axhline(0.5, color="#cccccc", lw=0.7, ls="--", alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(conf_order, rotation=15, ha="right")
    ax.set_ylim(0, 1.15)
    ax.set_xlabel("Transfer confidence tier\n(left = highest confidence)")
    ax.set_ylabel("Actual precision\n(fraction of predictions that are correct)")
    ax.set_title("Calibration: Confidence Tier vs Actual Precision\n"
                 "Well-calibrated = monotone decrease left → right",
                 fontweight="bold", pad=6)
    ax.legend(fontsize=7, handlelength=1.2)
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  ✓ Fig 3 → {out_path}")


def plot_fig4_permutation(perm_results: dict[str, dict],
                           categories:   list[str],
                           out_path:     Path) -> None:
    """
    Fig 4 — Permutation null distribution vs observed F1 for Method A.
    Each panel = one annotation category.
    Null = 1000 random human–mouse pair assignments.
    """
    n_cat  = len(categories)
    if n_cat == 0 or not perm_results:
        print("  ⚠  No permutation data — skipping Fig 4")
        return

    n_cols = min(4, n_cat)
    n_rows = int(np.ceil(n_cat / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(4.5*n_cols, 3.8*n_rows),
                              constrained_layout=True)
    axes = np.array(axes).flatten() if n_cat > 1 else [axes]

    for ax, cat in zip(axes, categories):
        res  = perm_results.get(cat, {})
        null = np.array(res.get("null_distribution", []))
        obs  = res.get("observed_f1", 0.0)
        pval = res.get("p_value", 1.0)
        if len(null) == 0:
            ax.text(0.5,0.5,"No data",ha="center",va="center",
                    transform=ax.transAxes,color="#999")
            ax.set_title(cat.replace("_"," "), fontweight="bold")
            continue
        ax.hist(null, bins=30, color=PAL["random"], edgecolor="white",
                linewidth=0.4, alpha=0.8, label="Null (permuted)")
        ax.axvline(obs, color=PAL["A"], lw=2.0, ls="--",
                   label=f"Observed F1={obs:.3f}")
        # shade p < 0.05 region
        p05_thresh = float(np.percentile(null, 95))
        ax.axvline(p05_thresh, color="#aaaaaa", lw=0.8, ls=":",
                   alpha=0.7, label=f"p=0.05 threshold")
        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "n.s."
        col = PAL["easy"] if pval < 0.05 else PAL["hard"]
        ax.text(0.97, 0.95, f"p={pval:.3f} {sig}",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=7.5, color=col, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", fc="white",
                          ec=col, lw=0.7, alpha=0.9))
        ax.set_xlabel("Permuted F1")
        ax.set_ylabel("Count")
        ax.set_title(cat.replace("_"," "), fontweight="bold", pad=4)
        if cat == categories[0]:
            ax.legend(fontsize=6, handlelength=1.2)

    # Hide unused axes
    for ax in axes[len(categories):]:
        ax.set_visible(False)

    fig.suptitle("Permutation Test: Method A vs Random Pair Assignment (n=1000)\n"
                 "* p<0.05   ** p<0.01   *** p<0.001   n.s. = not significant",
                 fontsize=9, fontweight="bold")
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  ✓ Fig 4 → {out_path}")


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main(args: argparse.Namespace) -> None:
    data_dir = Path(args.data)
    eval_dir = data_dir / "evaluation"
    eval_dir.mkdir(parents=True, exist_ok=True)

    N_BOOT = int(args.n_bootstrap)
    N_PERM = int(args.n_permutations)
    SEED   = int(args.seed)
    rng    = np.random.default_rng(SEED)

    print("=" * 60)
    print("  08 — Scientific Transfer Evaluation  (blind test set)")
    print("=" * 60)

    # ── Load files ────────────────────────────────────────────────────────────
    preds_path = data_dir / "transfer_predictions.tsv"
    gold_path  = data_dir / "mouse_annotations_gold.tsv"
    human_path = data_dir / "human_annotations_long.tsv"

    # Use blind test pairs
    test_pairs_path = data_dir / "ortholog_pairs_test.tsv"
    if test_pairs_path.exists():
        test_pairs = pd.read_csv(test_pairs_path, sep="\t", dtype=str).fillna("")
        print(f"\n  Blind test pairs: {len(test_pairs)} (from {test_pairs_path.name})")
    else:
        test_pairs_path = data_dir / "ortholog_pairs.tsv"
        test_pairs = pd.read_csv(test_pairs_path, sep="\t", dtype=str).fillna("") \
                     if test_pairs_path.exists() else pd.DataFrame()
        print(f"\n  NOTE: no test split found — evaluating on all pairs.")

    for f, label in [(preds_path,"transfer_predictions.tsv"),
                     (gold_path, "mouse_annotations_gold.tsv")]:
        if not f.exists():
            print(f"ERROR: {label} not found. Run earlier steps first.")
            return

    preds    = pd.read_csv(preds_path, sep="\t", dtype=str)
    gold_df  = pd.read_csv(gold_path,  sep="\t", dtype=str)
    human_df = pd.read_csv(human_path, sep="\t", dtype=str) \
               if human_path.exists() else pd.DataFrame()

    # Filter gold to test-set proteins only
    if not test_pairs.empty and "mouse_accession" in test_pairs.columns:
        test_accs = set(test_pairs["mouse_accession"].astype(str))
        gold_df   = gold_df[gold_df["protein_accession"].astype(str).isin(test_accs)]
        preds     = preds[preds["mouse_accession"].astype(str).isin(test_accs)]
        print(f"  Test-set mouse proteins: {gold_df['protein_accession'].nunique()}")

    # Assign difficulty to test pairs
    if not test_pairs.empty:
        test_pairs = assign_difficulty(test_pairs)

    gold_sets  = _build_gold_sets(gold_df)
    methods    = [m for m in preds["method"].unique().tolist()
                  if m != "control_random"]
    categories = sorted(gold_df["annotation_category"].unique().tolist()
                        if "annotation_category" in gold_df.columns else [])

    print(f"\n  Mouse proteins (gold):  {len(gold_sets)}")
    print(f"  Methods:                {methods}")
    print(f"  Categories:             {categories}")
    print(f"  Bootstrap resamples:    {N_BOOT}")
    print(f"  Permutation trials:     {N_PERM}")

    # ── Per-protein evaluation ────────────────────────────────────────────────
    all_rows: list[dict] = []
    for method in methods:
        print(f"\nEvaluating {method} …")
        rows = evaluate_method(method, preds, gold_sets, categories)
        all_rows.extend(rows)
        df_m = pd.DataFrame(rows)
        for cat in categories:
            sub = df_m[df_m["category"] == cat]
            if sub.empty:
                continue
            print(f"  [{cat:<20}]  R={sub['recall'].mean():.3f}  "
                  f"P={sub['precision'].mean():.3f}  "
                  f"F1={sub['f1'].mean():.3f}  "
                  f"J={sub['jaccard'].mean():.3f}")

    per_protein = pd.DataFrame(all_rows)
    per_protein.to_csv(eval_dir / "evaluation_per_protein.tsv", sep="\t", index=False)

    # ── Summary table ─────────────────────────────────────────────────────────
    if not per_protein.empty:
        summary = (per_protein
                   .groupby(["method","category"])
                   [["precision","recall","f1","jaccard"]]
                   .mean().round(4).reset_index())
        summary.to_csv(eval_dir / "evaluation_summary.tsv", sep="\t", index=False)
        print("\n=== SUMMARY ===")
        print(summary.to_string(index=False))
    else:
        summary = pd.DataFrame()

    # ── Bootstrap CIs ─────────────────────────────────────────────────────────
    print(f"\nComputing bootstrap CIs (n={N_BOOT}) …")
    ci_rows: list[dict] = []
    for method in methods:
        df_m = per_protein[per_protein["method"] == method]
        for cat in categories:
            df_c = df_m[df_m["category"] == cat]
            if df_c.empty:
                continue
            for metric in ["precision","recall","f1","jaccard"]:
                vals          = df_c[metric].astype(float).tolist()
                mean_, lo, hi = bootstrap_ci(vals, n=N_BOOT, rng=rng)
                ci_rows.append({
                    "method": method, "category": cat, "metric": metric,
                    "mean": round(mean_,4), "ci_lower": round(lo,4),
                    "ci_upper": round(hi,4),
                })
    ci_df = pd.DataFrame(ci_rows)
    ci_df.to_csv(eval_dir / "evaluation_bootstrap_ci.tsv", sep="\t", index=False)

    # ── Random control ────────────────────────────────────────────────────────
    print(f"\nRandom-pair control ({args.n_random} draws per protein) …")
    ctrl_rows = random_control(
        gold_sets, human_df, categories,
        n_random=int(args.n_random), seed=SEED)
    ctrl_df = pd.DataFrame(ctrl_rows)
    ctrl_df.to_csv(eval_dir / "evaluation_controls.tsv", sep="\t", index=False)
    if ctrl_rows:
        ctrl_summary = (ctrl_df
                        .groupby("category")[["jaccard"]]
                        .mean().round(4))
        print("  Random control Jaccard per category:")
        print(ctrl_summary.to_string())

    # ── Difficulty-stratified evaluation ──────────────────────────────────────
    print("\nDifficulty-stratified evaluation …")
    diff_rows: list[dict] = []
    if not test_pairs.empty and "difficulty" in test_pairs.columns:
        for method in methods:
            rows = evaluate_by_difficulty(method, preds, gold_sets,
                                           test_pairs, categories)
            diff_rows.extend(rows)
        diff_df = pd.DataFrame(diff_rows)
        diff_df.to_csv(eval_dir / "evaluation_difficulty.tsv", sep="\t", index=False)
        print("  Rows written:", len(diff_rows))
        for tier in DIFFICULTY_ORDER:
            sub = diff_df[diff_df["difficulty"] == tier]
            if sub.empty:
                continue
            n_pairs = len(test_pairs[test_pairs["difficulty"]==tier])
            print(f"\n  [{tier.upper()} — {n_pairs} pairs]")
            for method in methods:
                sm = sub[sub["method"]==method]
                if sm.empty:
                    continue
                print(f"    {method}: F1={sm['f1'].mean():.3f}")
    else:
        diff_df = pd.DataFrame()
        print("  Skipped (no difficulty column in test pairs — run step 05 for split)")

    # ── Permutation test ──────────────────────────────────────────────────────
    print(f"\nPermutation test, Method A (n={N_PERM}) …")
    perm_results: dict[str, dict] = {}
    human_ann_dict = {}
    if not human_df.empty:
        if "is_transferable" in human_df.columns:
            h_filt = human_df[human_df["is_transferable"].astype(str)
                              .str.lower().isin(["true","yes","1"])]
        else:
            h_filt = human_df
        for acc, grp in h_filt.groupby("protein_accession"):
            human_ann_dict[str(acc)] = grp.reset_index(drop=True)

    if not test_pairs.empty and human_ann_dict:
        perm_m_accs = test_pairs["mouse_accession"].values.copy()
        real_preds  = _build_pred_sets_from_method("A_direct_transfer",
                                                    test_pairs, human_ann_dict)
        real_f1     = _mean_f1_per_cat(real_preds, gold_sets, categories)

        perm_f1: dict[str, list[float]] = {cat: [] for cat in categories}
        perm_rng = np.random.default_rng(SEED + 1)
        for trial in range(N_PERM):
            if (trial + 1) % 200 == 0:
                print(f"  Permutation trial {trial+1}/{N_PERM} …")
            shuffled = test_pairs.copy()
            shuffled["mouse_accession"] = perm_rng.permutation(perm_m_accs)
            p_preds = _build_pred_sets_from_method("A_direct_transfer",
                                                    shuffled, human_ann_dict)
            pf1     = _mean_f1_per_cat(p_preds, gold_sets, categories)
            for cat in categories:
                perm_f1[cat].append(pf1.get(cat, 0.0))

        perm_rows = []
        for cat in categories:
            obs  = real_f1.get(cat, 0.0)
            null = np.array(perm_f1[cat])
            pval = float(np.mean(null >= obs))
            perm_results[cat] = {
                "observed_f1":       round(obs, 4),
                "null_mean":         round(float(np.mean(null)), 4),
                "null_std":          round(float(np.std(null)), 4),
                "p_value":           round(pval, 4),
                "null_distribution": null.tolist(),
            }
            sig = "***" if pval<0.001 else "**" if pval<0.01 else "*" if pval<0.05 else "n.s."
            print(f"  [{cat:<20}] obs={obs:.3f}  null_mean={np.mean(null):.3f}  "
                  f"p={pval:.4f} {sig}")

        perm_summary = pd.DataFrame([
            {"method":"A_direct_transfer","category":cat,
             "observed_f1":v["observed_f1"],"null_mean":v["null_mean"],
             "null_std":v["null_std"],"p_value":v["p_value"]}
            for cat, v in perm_results.items()
        ])
        perm_summary.to_csv(eval_dir / "evaluation_permutation.tsv",
                            sep="\t", index=False)

    # ── Calibration ───────────────────────────────────────────────────────────
    print("\nComputing calibration …")
    calib_df = compute_calibration(preds, gold_sets, methods)
    calib_df.to_csv(eval_dir / "evaluation_calibration.tsv", sep="\t", index=False)

    # ── Figures ───────────────────────────────────────────────────────────────
    print("\nGenerating publication-quality figures …")
    if not summary.empty:
        # Add control to summary for Fig 1
        full_summary = summary.copy()
        if ctrl_rows:
            ctrl_s = ctrl_df.groupby("category")[["f1","jaccard"]].mean().round(4).reset_index()
            ctrl_s["method"]    = "control_random"
            ctrl_s["precision"] = 0.0
            ctrl_s["recall"]    = 0.0
            full_summary = pd.concat([full_summary, ctrl_s], ignore_index=True)

        plot_fig1_main_results(
            full_summary, ci_df,
            eval_dir / "fig1_main_results.png")

    if not diff_df.empty:
        plot_fig2_difficulty_stratified(
            diff_df, eval_dir / "fig2_difficulty_stratified.png")

    if not calib_df.empty:
        plot_fig3_calibration(calib_df, eval_dir / "fig3_calibration.png")

    if perm_results:
        plot_fig4_permutation(perm_results, categories,
                              eval_dir / "fig4_permutation.png")

    # ── Final interpretation ───────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("  FINAL EVALUATION SUMMARY")
    print("=" * 60)
    if not summary.empty:
        for method in methods:
            sub = summary[summary["method"] == method]
            if sub.empty:
                continue
            mr = sub["recall"].mean()
            mp = sub["precision"].mean()
            mf = sub["f1"].mean()
            print(f"\n  {method}:")
            print(f"    Mean F1={mf:.3f}  Recall={mr:.3f}  Precision={mp:.3f}")
            if mf >= 0.7:
                print("    ► High performance — ontology transfers well to mouse")
            elif mf >= 0.4:
                print("    ► Moderate performance — partial coverage")
            else:
                print("    ► Low performance — ontology may need mouse expansion")

    print(f"\n  All outputs → {eval_dir}/")
    print("  Review manually: failure cases in evaluation_per_protein.tsv")
    print("  Method comparison: fig1_main_results.png")


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Scientific evaluation of transfer predictions")
    p.add_argument("--data",           default="transfer_pipeline/data")
    p.add_argument("--n_bootstrap",    default=1000, type=int,
                   help="Bootstrap resamples for 95%% CI (default 1000)")
    p.add_argument("--n_permutations", default=1000, type=int,
                   help="Permutation test trials (default 1000)")
    p.add_argument("--n_random",       default=100,  type=int,
                   help="Random-pair control draws per protein")
    p.add_argument("--seed",           default=42,   type=int)
    main(p.parse_args())
