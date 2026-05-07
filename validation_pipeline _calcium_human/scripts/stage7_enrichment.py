"""
scripts/stage7_enrichment.py
────────────────────────────────────────────────────────────────────────────
STAGE 7 — DOMAIN ENRICHMENT ANALYSIS
Fisher's exact test + BH-FDR for cAMP foreground vs background proteins.
Produces a publication-style horizontal bar plot (matplotlib only).
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys as _sys
import os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from scripts.utils import (benjamini_hochberg, get_input_path, load_config,
                            load_table, make_output_dir, rename_to_canonical,
                            save_tsv, setup_logging, split_cell, stage_banner)


def _build_domain_protein_map(df: pd.DataFrame,
                               domain_col: str = "domain_type_id",
                               protein_col: str = "protein_accession",
                               delimiter: str = ";") -> dict[str, set]:
    """Map domain_type_id → set of protein accessions."""
    d2p: dict[str, set] = {}
    for _, row in df.iterrows():
        prot = str(row.get(protein_col, "")).strip()
        for dt in split_cell(row.get(domain_col, ""), delimiter):
            d2p.setdefault(dt, set()).add(prot)
    return d2p


def run_enrichment(fg_proteins: set, bg_proteins: set,
                   d2p_fg: dict, d2p_bg: dict,
                   min_fg: int = 2) -> pd.DataFrame:
    """
    Fisher's exact test for each domain type.
    fg_proteins / bg_proteins: sets of accessions.
    d2p_fg / d2p_bg: {domain → set of proteins} for each group.
    """
    all_domains = set(d2p_fg.keys()) | set(d2p_bg.keys())
    n_fg = len(fg_proteins)
    n_bg = len(bg_proteins)

    rows = []
    for domain in all_domains:
        fg_hit = len(d2p_fg.get(domain, set()) & fg_proteins)
        bg_hit = len(d2p_bg.get(domain, set()) & bg_proteins)
        fg_miss = n_fg - fg_hit
        bg_miss = n_bg - bg_hit

        if fg_hit < min_fg:
            continue

        table = [[fg_hit, fg_miss], [bg_hit, bg_miss]]
        try:
            odds_ratio, pvalue = fisher_exact(table, alternative="greater")
        except Exception:
            odds_ratio, pvalue = float("nan"), 1.0

        rows.append({
            "domain_type": domain,
            "fg_count": fg_hit, "fg_total": n_fg,
            "bg_count": bg_hit, "bg_total": n_bg,
            "fg_freq": round(fg_hit / n_fg, 4) if n_fg else 0,
            "bg_freq": round(bg_hit / n_bg, 4) if n_bg else 0,
            "odds_ratio": round(odds_ratio, 4),
            "pvalue": pvalue,
        })

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df["fdr"] = benjamini_hochberg(df["pvalue"].tolist())
    df = df.sort_values("fdr").reset_index(drop=True)
    return df


def _enrichment_barplot(df: pd.DataFrame, out_path: Path,
                         top_n: int = 20, fdr_thr: float = 0.05,
                         dpi: int = 150, font_size: int = 9) -> None:
    sig = df[df["fdr"] < fdr_thr].head(top_n).copy()
    if sig.empty:
        return
    sig = sig.sort_values("odds_ratio")
    labels = sig["domain_type"].tolist()
    ors    = sig["odds_ratio"].tolist()
    fdrs   = sig["fdr"].tolist()

    fig, ax = plt.subplots(figsize=(9, max(4, len(labels) * 0.4)), dpi=dpi)
    colors = ["#E05C5C" if f < 0.01 else "#F0A050" for f in fdrs]
    y_pos = range(len(labels))
    bars = ax.barh(list(y_pos), ors, color=colors, edgecolor="white", height=0.65)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(labels, fontsize=font_size)
    ax.axvline(x=1, color="#666666", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Odds Ratio", fontsize=font_size)
    ax.set_title("Top Enriched Domain Types (FDR < 0.05)",
                 fontsize=font_size + 2, fontweight="bold")
    from matplotlib.patches import Patch
    legend = [Patch(facecolor="#E05C5C", label="FDR < 0.01"),
              Patch(facecolor="#F0A050", label="FDR < 0.05")]
    ax.legend(handles=legend, fontsize=font_size - 1, loc="lower right")
    for bar, fdr in zip(bars, fdrs):
        ax.text(bar.get_width() * 1.01,
                bar.get_y() + bar.get_height() / 2,
                f"FDR={fdr:.2e}", va="center", fontsize=font_size - 2)
    ax.spines[["top","right"]].set_visible(False)
    plt.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def run(cfg: dict) -> None:
    logger = setup_logging(cfg.get("log_dir","logs"), "stage7_enrichment")
    stage_banner(logger, "STAGE 7", "Domain Enrichment Analysis")
    out_dir = make_output_dir(cfg, "stage7_enrichment")

    enr_cfg = cfg.get("enrichment", {})
    fdr_thr  = enr_cfg.get("fdr_threshold", 0.05)
    min_fg   = enr_cfg.get("min_foreground_count", 2)
    top_n    = enr_cfg.get("top_n_plot", 20)
    plot_cfg = cfg.get("plotting", {})

    fg_path = get_input_path(cfg, "foreground_protein_domains")
    bg_path = get_input_path(cfg, "background_protein_domains")
    occ_path = get_input_path(cfg, "domain_occurrence_table")

    # ── Resolve input tables ───────────────────────────────────────────────
    if fg_path and bg_path:
        fg_df = load_table(fg_path)
        bg_df = load_table(bg_path)
        prot_col = next((c for c in ["protein_accession","accession","uniprot_id"]
                         if c in fg_df.columns), None)
        dom_col  = next((c for c in ["domain_type_id","domain_id"]
                         if c in fg_df.columns), None)
        if not prot_col or not dom_col:
            logger.error(f"Cannot find protein/domain columns in {fg_path}")
            return
        fg_proteins = set(fg_df[prot_col].dropna().str.strip())
        bg_proteins = (set(bg_df.get(prot_col, pd.Series()).dropna().str.strip())
                       | fg_proteins)
        d2p_fg = _build_domain_protein_map(fg_df, dom_col, prot_col)
        d2p_bg = _build_domain_protein_map(
            pd.concat([fg_df, bg_df], ignore_index=True), dom_col, prot_col)
        logger.info(f"Foreground: {len(fg_proteins)} proteins  "
                    f"Background: {len(bg_proteins)} proteins")
    elif occ_path:
        logger.warning(
            "foreground_protein_domains / background_protein_domains not found. "
            "Falling back to occurrence table — all proteins treated as foreground. "
            "Provide foreground_protein_domains.tsv + background_protein_domains.tsv "
            "for biologically meaningful enrichment."
        )
        occ_full = load_table(occ_path)
        occ_full = rename_to_canonical(occ_full, cfg.get("column_aliases", {}))
        prot_col = "protein_accession"
        dom_col  = "domain_type_id"
        fg_proteins = set(occ_full[prot_col].dropna().str.strip())
        bg_proteins = fg_proteins
        d2p_fg = _build_domain_protein_map(occ_full, dom_col, prot_col)
        d2p_bg = d2p_fg
        logger.info(f"Occurrence table fallback: {len(fg_proteins)} proteins  "
                    f"{len(d2p_fg)} domain types")
    else:
        logger.warning("No input data found — skipping Stage 7")
        return

    results = run_enrichment(fg_proteins, bg_proteins, d2p_fg, d2p_bg, min_fg)
    if results.empty:
        logger.warning("No enrichment results (insufficient data)")
        save_tsv(results, out_dir / "enrichment_results.tsv", logger)
        return

    save_tsv(results, out_dir / "enrichment_results.tsv", logger)
    n_sig = (results["fdr"] < fdr_thr).sum()
    logger.info(f"Significant at FDR<{fdr_thr}: {n_sig}/{len(results)} domains")

    _enrichment_barplot(results, out_dir / "top_enriched_domains.png",
                        top_n, fdr_thr,
                        plot_cfg.get("figure_dpi", 150),
                        plot_cfg.get("font_size", 9))
    logger.info("Stage 7 complete.")


if __name__ == "__main__":
    import sys, os
    # Always run from the ontology_pipeline project root so that
    # relative paths (input/, output/, config/) resolve correctly
    # regardless of where PyCharm / the terminal launched the script from.
    _project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(_project_root)
    _cfg_path = sys.argv[1] if len(sys.argv) > 1 else "config/pipeline_config.yaml"
    run(load_config(_cfg_path))
