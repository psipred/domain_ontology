"""
scripts/stage10_robustness.py
────────────────────────────────────────────────────────────────────────────
STAGE 10 — ROBUSTNESS TESTING
Bootstrap resampling, leave-one-family-out, and label-shuffling controls.
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys as _sys
import os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from scripts.utils import (benjamini_hochberg, get_input_path, load_config,
                            load_table, make_output_dir, rename_to_canonical,
                            save_tsv, setup_logging, split_cell, stage_banner)
from scripts.stage7_enrichment import run_enrichment, _build_domain_protein_map


def run(cfg: dict) -> None:
    logger = setup_logging(cfg.get("log_dir","logs"), "stage10_robustness")
    stage_banner(logger, "STAGE 10", "Robustness Testing")
    out_dir = make_output_dir(cfg, "stage10_robustness")

    rob_cfg = cfg.get("robustness", {})
    n_boot  = rob_cfg.get("bootstrap_iterations", 100)
    frac    = rob_cfg.get("bootstrap_fraction", 0.8)
    n_shuf  = rob_cfg.get("shuffles", 100)
    seed    = cfg.get("random_seed", 42)
    rng     = np.random.default_rng(seed)
    plot_cfg = cfg.get("plotting", {})
    dpi = plot_cfg.get("figure_dpi", 150)

    fg_path  = get_input_path(cfg, "foreground_protein_domains")
    bg_path  = get_input_path(cfg, "background_protein_domains")
    grp_path = get_input_path(cfg, "protein_group_labels")

    # Fallback: use occurrence table when no dedicated foreground file exists
    if not fg_path:
        occ_path = get_input_path(cfg, "domain_occurrence_table")
        if not occ_path:
            logger.warning("No foreground_protein_domains or domain_occurrence_table — skipping Stage 10")
            return
        logger.info("foreground_protein_domains not found — using domain_occurrence_table as foreground")
        fg_path = occ_path

    fg_df = load_table(fg_path)
    bg_df = load_table(bg_path) if bg_path else None
    fg_df = rename_to_canonical(fg_df, cfg.get("column_aliases", {}))

    prot_col = next((c for c in ["protein_accession","accession"] if c in fg_df.columns), None)
    dom_col  = next((c for c in ["domain_type_id","domain_id"] if c in fg_df.columns), None)
    if not prot_col or not dom_col:
        logger.error("Cannot find protein/domain columns")
        return

    fg_proteins = list(fg_df[prot_col].dropna().unique())
    if bg_df is not None:
        bg_df = rename_to_canonical(bg_df, cfg.get("column_aliases", {}))
        bg_proteins = list(set(bg_df[prot_col].dropna().unique()) | set(fg_proteins))
        d2p_bg = _build_domain_protein_map(
            pd.concat([fg_df, bg_df], ignore_index=True), dom_col, prot_col)
    else:
        bg_proteins = fg_proteins
        d2p_bg = _build_domain_protein_map(fg_df, dom_col, prot_col)

    # ── Bootstrap resampling ─────────────────────────────────────────────────
    logger.info(f"Bootstrap: {n_boot} iterations, fraction={frac} …")
    boot_rows = []
    for i in range(n_boot):
        sample = rng.choice(fg_proteins, size=int(len(fg_proteins)*frac), replace=False)
        d2p_sample = _build_domain_protein_map(
            fg_df[fg_df[prot_col].isin(sample)], dom_col, prot_col)
        res = run_enrichment(set(sample), set(bg_proteins), d2p_sample, d2p_bg, min_fg=1)
        if not res.empty:
            res["iteration"] = i
            boot_rows.append(res[["domain_type","odds_ratio","pvalue","fdr","iteration"]])
    if boot_rows:
        boot_df = pd.concat(boot_rows, ignore_index=True)
        stability = boot_df.groupby("domain_type").agg(
            mean_or=("odds_ratio","mean"),
            std_or=("odds_ratio","std"),
            n_detected=("odds_ratio","count"),
            detection_rate=("odds_ratio", lambda x: len(x)/n_boot),
        ).reset_index()
        save_tsv(stability, out_dir / "bootstrap_stability.tsv", logger)
    else:
        save_tsv(pd.DataFrame(), out_dir / "bootstrap_stability.tsv", logger)

    # ── Leave-one-family-out ─────────────────────────────────────────────────
    lofo_rows = []
    if grp_path:
        grp_df = load_table(grp_path)
        grp_col = next((c for c in ["protein_group","group","protein_family"]
                        if c in grp_df.columns), None)
        p_col2  = next((c for c in ["protein_accession","accession"]
                        if c in grp_df.columns), None)
        if grp_col and p_col2:
            for grp in grp_df[grp_col].dropna().unique():
                left_out = set(grp_df[grp_df[grp_col]==grp][p_col2].dropna())
                remaining = [p for p in fg_proteins if p not in left_out]
                d2p_r = _build_domain_protein_map(
                    fg_df[fg_df[prot_col].isin(remaining)], dom_col, prot_col)
                res = run_enrichment(set(remaining), set(bg_proteins), d2p_r, d2p_bg, min_fg=1)
                if not res.empty:
                    res["left_out_group"] = grp
                    lofo_rows.append(res.head(10)[["domain_type","odds_ratio","fdr","left_out_group"]])
    save_tsv(pd.concat(lofo_rows, ignore_index=True) if lofo_rows else pd.DataFrame(),
             out_dir / "leave_one_family_out.tsv", logger)

    # ── Label shuffling ──────────────────────────────────────────────────────
    logger.info(f"Label shuffling: {n_shuf} iterations …")
    shuf_rows = []
    all_domains = list(d2p_bg.keys())
    for i in range(n_shuf):
        shuffled_fg = fg_df.copy()
        shuffled_fg[dom_col] = rng.permutation(shuffled_fg[dom_col].values)
        d2p_shuf = _build_domain_protein_map(shuffled_fg, dom_col, prot_col)
        res = run_enrichment(set(fg_proteins), set(bg_proteins), d2p_shuf, d2p_bg, min_fg=1)
        if not res.empty:
            shuf_rows.append({"iteration": i,
                              "n_sig_fdr05": (res["fdr"]<0.05).sum(),
                              "max_or": res["odds_ratio"].max()})
    shuf_df = pd.DataFrame(shuf_rows)
    save_tsv(shuf_df, out_dir / "shuffled_controls.tsv", logger)

    # ── Robustness summary ───────────────────────────────────────────────────
    summary_rows = []
    if not shuf_df.empty:
        summary_rows.append({"metric":"shuffled_mean_sig_fdr05",
                              "value": shuf_df["n_sig_fdr05"].mean()})
        summary_rows.append({"metric":"shuffled_max_or_mean",
                              "value": shuf_df["max_or"].mean()})
    save_tsv(pd.DataFrame(summary_rows), out_dir/"robustness_summary.tsv", logger)

    # ── Plots ────────────────────────────────────────────────────────────────
    if boot_rows and "detection_rate" in stability.columns:
        top_stable = stability.nlargest(20, "detection_rate")
        fig, ax = plt.subplots(figsize=(8, max(4, len(top_stable)*0.3)), dpi=dpi)
        ax.barh(top_stable["domain_type"], top_stable["detection_rate"],
                color="#4C72B0", edgecolor="white")
        ax.axvline(0.8, color="red", linestyle="--", linewidth=0.8,
                   label="80% threshold")
        ax.set_xlabel("Bootstrap Detection Rate")
        ax.set_title("Bootstrap Stability (top 20 domains)",
                     fontweight="bold")
        ax.legend()
        ax.spines[["top","right"]].set_visible(False)
        plt.tight_layout()
        fig.savefig(out_dir/"robustness_plots.png", dpi=dpi, bbox_inches="tight")
        plt.close(fig)

    logger.info("Stage 10 complete.")


if __name__ == "__main__":
    import sys, os
    # Always run from the ontology_pipeline project root so that
    # relative paths (input/, output/, config/) resolve correctly
    # regardless of where PyCharm / the terminal launched the script from.
    _project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(_project_root)
    _cfg_path = sys.argv[1] if len(sys.argv) > 1 else "config/pipeline_config.yaml"
    run(load_config(_cfg_path))
