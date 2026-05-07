"""
scripts/stage12_benchmark.py
────────────────────────────────────────────────────────────────────────────
STAGE 12 — KNOWN-BIOLOGY BENCHMARK  (complete fixed version)

Fixes applied
─────────────
  FIX 1  _REL_ALIAS  — maps OWL property names to pipeline column names
  FIX 2  _norm_val   — strips leading DOT:/DOP: ID tokens from gold values
  FIX 3  All 6 relationship types benchmarked (was only 2)
  FIX 4  Clear warning when auto-generating gold standard (trivial scores)
  FIX 5  Improved figure with value annotations and interpretation note
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations
import re
from collections import defaultdict
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from scripts.utils import (get_input_path, load_config, load_table, make_output_dir,
                            rename_to_canonical, save_tsv, setup_logging,
                            split_cell, stage_banner)


# FIX 1 — maps OWL property names → pipeline relationship_type names
_REL_ALIAS: dict[str, str] = {
    "hasFunctionTag": "function_tag", "function_tags": "function_tag",
    "hasFunctionalTag": "function_tag", "functionalTag": "function_tag",
    "hasBindingTargetCategory": "binding_target", "hasBindingTarget": "binding_target",
    "binding_parnter": "binding_target", "binding_partner": "binding_target",
    "bindingTarget": "binding_target",
    "hasPositionalCategory": "positional_category", "positionalCategory": "positional_category",
    "hasCopyNumberProperty": "copy_number_property", "copyNumberProperty": "copy_number_property",
    "copy_number": "copy_number_property",
    "hasTopologyContext": "topology_context", "topologyContext": "topology_context",
    "topology": "topology_context",
    "hasProximityCategory": "proximity_categories", "proximityCategory": "proximity_categories",
    "proximity_category": "proximity_categories",
}

# FIX 2 — strip leading "DOT:000079 " or "DOP:000003 " from gold values
_ID_PREFIX_RE = re.compile(r"^(?:DOT|DOP|DOF|GO|CHEBI|R-[A-Z]{3}):[0-9A-Za-z\-]+\s+")

def _norm_val(v: str) -> str:
    return _ID_PREFIX_RE.sub("", str(v).strip()).strip()


# FIX 3 — all 6 relationship types
def _extract_predictions(occ: pd.DataFrame, cfg: dict) -> dict[str, set[tuple]]:
    occ = rename_to_canonical(occ, cfg.get("column_aliases", {}))
    bind_col = next((c for c in ["binding_partner","binding_parnter"] if c in occ.columns), "binding_partner")
    column_map = {
        "function_tag":         "function_tags",
        "binding_target":       bind_col,
        "positional_category":  "positional_category",
        "copy_number_property": "copy_number_property",
        "topology_context":     "topology_context",
        "proximity_categories": "proximity_categories",
    }
    preds: dict[str, set] = defaultdict(set)
    dom_col = "domain_type_id"
    if dom_col not in occ.columns:
        return preds
    for rel_type, col in column_map.items():
        if col not in occ.columns:
            continue
        for _, row in occ.iterrows():
            for dt in split_cell(row.get(dom_col, "")):
                for val in split_cell(row.get(col, "")):
                    # Normalise to bare ID so pred matches gold format
                    preds[rel_type].add((dt.strip(), _norm_val(val)))
    return preds


def _auto_gold_standard(occ: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Fallback: derive gold standard from occurrence table itself (gives trivial P=R=F1=1.0)."""
    occ = rename_to_canonical(occ, cfg.get("column_aliases", {}))
    dom_col  = "domain_type_id"
    bind_col = next((c for c in ["binding_partner","binding_parnter"] if c in occ.columns), None)
    column_map = {
        "function_tag": "function_tags", "binding_target": bind_col,
        "positional_category": "positional_category",
        "copy_number_property": "copy_number_property", "topology_context": "topology_context",
    }
    rows: list[dict] = []
    if dom_col not in occ.columns:
        return pd.DataFrame(columns=["domain_type","relationship_type","expected_value"])
    for rel_type, col in column_map.items():
        if not col or col not in occ.columns:
            continue
        seen: set[tuple] = set()
        for _, row in occ.iterrows():
            for dt in split_cell(row.get(dom_col, "")):
                for val in split_cell(row.get(col, "")):
                    key = (dt.strip(), rel_type, val.strip())
                    if key not in seen:
                        seen.add(key)
                        rows.append({"domain_type": dt.strip(),
                                     "relationship_type": rel_type,
                                     "expected_value": val.strip()})
    return pd.DataFrame(rows)


def _prf(n_tp: int, n_pred: int, n_gold: int) -> tuple[float, float, float]:
    p = n_tp / n_pred if n_pred else 0.0
    r = n_tp / n_gold if n_gold else 0.0
    f = 2*p*r/(p+r) if (p+r) else 0.0
    return round(p, 4), round(r, 4), round(f, 4)


def run_benchmark(gold_df: pd.DataFrame, preds: dict[str, set[tuple]]) -> dict[str, pd.DataFrame]:
    required = {"domain_type", "relationship_type", "expected_value"}
    if not required.issubset(gold_df.columns):
        raise ValueError(f"gold_standard must have columns: {required}")
    matched_rows, missed_rows, unexpected_rows, summary_rows = [], [], [], []
    for rel_type, grp in gold_df.groupby("relationship_type"):
        gold_set = {(str(r["domain_type"]).strip(), str(r["expected_value"]).strip())
                    for _, r in grp.iterrows()}
        pred_set = preds.get(rel_type, set())
        tp, fn, fp = gold_set & pred_set, gold_set - pred_set, pred_set - gold_set
        prec, rec, f1 = _prf(len(tp), len(pred_set), len(gold_set))
        summary_rows.append({"relationship_type": rel_type, "n_gold": len(gold_set),
                              "n_predicted": len(pred_set), "n_true_pos": len(tp),
                              "n_false_neg": len(fn), "n_false_pos": len(fp),
                              "precision": prec, "recall": rec, "f1": f1})
        for dt, val in sorted(tp):
            matched_rows.append({"relationship_type": rel_type, "domain_type": dt, "value": val, "status": "MATCHED"})
        for dt, val in sorted(fn):
            missed_rows.append({"relationship_type": rel_type, "domain_type": dt, "value": val, "status": "MISSED"})
        for dt, val in sorted(fp):
            unexpected_rows.append({"relationship_type": rel_type, "domain_type": dt, "value": val, "status": "UNEXPECTED"})
    return {"matched": pd.DataFrame(matched_rows), "missed": pd.DataFrame(missed_rows),
            "unexpected": pd.DataFrame(unexpected_rows), "summary": pd.DataFrame(summary_rows)}


def _benchmark_plot(summary: pd.DataFrame, out_path: Path, dpi: int = 150, font_size: int = 9,
                    matched_df=None, missed_df=None, unexpected_df=None) -> None:
    if summary.empty:
        return
    rel_types = summary["relationship_type"].tolist()
    n = len(rel_types)
    x = list(range(n))
    fig, axes = plt.subplots(1, 2, figsize=(max(14, n*3.0), 6), dpi=dpi)

    # Left: P/R/F1
    ax = axes[0]
    w  = 0.25
    for offset, (metric, col, colour) in enumerate([
        ("Precision","precision","#4C72B0"),
        ("Recall","recall","#55A868"),
        ("F1","f1","#C44E52")
    ]):
        vals = summary[col].tolist()
        bars = ax.bar([i+(offset-1)*w for i in x], vals, w,
                      label=metric, color=colour, alpha=0.88)
        for bar, v in zip(bars, vals):
            if v > 0.01:
                ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.02,
                        f"{v:.2f}", ha="center", va="bottom",
                        fontsize=max(6, font_size-2))
    ax.set_xticks(x)
    ax.set_xticklabels(rel_types, fontsize=font_size, rotation=20, ha="right")
    ax.set_ylim(0, 1.30)
    ax.set_ylabel("Score", fontsize=font_size)
    ax.set_title("Precision / Recall / F1\nvs gold-standard relationships",
                 fontsize=font_size+1, fontweight="bold")
    ax.legend(fontsize=font_size)
    ax.spines[["top","right"]].set_visible(False)
    ax.text(0.02, 0.97,
            "High Recall = ontology covers curated biology\n"
            "Low Precision = ontology makes additional predictions\n"
            "(may be correct but not yet in gold standard)",
            transform=ax.transAxes, fontsize=max(6, font_size-2), va="top",
            color="#555", bbox=dict(boxstyle="round,pad=0.4", fc="#f9f9f9", ec="#ccc", alpha=0.9))

    # Right: stacked coverage
    ax2 = axes[1]
    mc = summary["n_true_pos"].tolist()
    sc = summary["n_false_neg"].tolist()
    uc = summary["n_false_pos"].tolist()
    b1 = ax2.bar(x, mc, label="Matched (True Positive)", color="#55A868", alpha=0.88)
    b2 = ax2.bar(x, sc, bottom=mc, label="Missed (in gold, not predicted)", color="#C44E52", alpha=0.88)
    bot2 = [m+s for m,s in zip(mc, sc)]
    b3 = ax2.bar(x, uc, bottom=bot2, label="Unexpected (predicted, not in gold)", color="#DD8452", alpha=0.88)
    for bars, bots in [(b1,[0]*n),(b2,mc),(b3,bot2)]:
        for bar, bot in zip(bars, bots):
            h = bar.get_height()
            if h >= 5:
                ax2.text(bar.get_x()+bar.get_width()/2, bot+h/2, str(int(h)),
                         ha="center", va="center", fontsize=max(6,font_size-1),
                         color="white", fontweight="bold")
    ax2.set_xticks(x)
    ax2.set_xticklabels(rel_types, fontsize=font_size, rotation=20, ha="right")
    ax2.set_ylabel("Relationship count", fontsize=font_size)
    ax2.set_title("Coverage Breakdown\n(matched · missed · unexpected)",
                  fontsize=font_size+1, fontweight="bold")
    ax2.legend(fontsize=max(6,font_size-1), loc="upper right")
    ax2.spines[["top","right"]].set_visible(False)

    plt.suptitle("Known-Biology Benchmark: cAMP Pathway Domain Ontology",
                 fontsize=font_size+3, fontweight="bold", y=1.01)
    plt.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    if unexpected_df is not None and not unexpected_df.empty and "value" in unexpected_df.columns:
        unex_path = out_path.parent / "benchmark_unexpected_plot.png"
        ux = (unexpected_df.groupby(["relationship_type","domain_type"])
              .size().reset_index(name="count")
              .sort_values("count", ascending=False).head(30))
        if not ux.empty:
            labels = [f"{r['domain_type'].split(':')[-1][:18]} [{r['relationship_type']}]"
                      for _, r in ux.iterrows()]
            fig3, ax3 = plt.subplots(figsize=(12, max(5, len(labels)*0.38)), dpi=dpi)
            ax3.barh(labels[::-1], ux["count"].tolist()[::-1], color="#DD8452", edgecolor="white", alpha=0.88)
            ax3.set_xlabel("Number of unexpected relationship instances", fontsize=font_size)
            ax3.set_title("Unexpected Relationships\n(predicted by ontology, absent from gold — may be correct biology)",
                          fontsize=font_size+1, fontweight="bold")
            ax3.spines[["top","right"]].set_visible(False)
            plt.tight_layout()
            fig3.savefig(unex_path, dpi=dpi, bbox_inches="tight")
            plt.close(fig3)


def run(cfg: dict) -> None:
    logger   = setup_logging(cfg.get("log_dir","logs"), "stage12_benchmark")
    stage_banner(logger, "STAGE 12", "Known-Biology Benchmark")
    out_dir  = make_output_dir(cfg, "stage12_benchmark")
    plot_cfg = cfg.get("plotting", {})

    occ_path = get_input_path(cfg, "domain_occurrence_table")
    if not occ_path:
        logger.warning("domain_occurrence_table not found — skipping Stage 12")
        return
    occ = load_table(occ_path)

    gold_path  = get_input_path(cfg, "gold_standard_relationships")
    using_auto = False

    if not gold_path or not Path(gold_path).exists():
        logger.warning(
            "gold_standard_relationships.tsv not found.\n"
            "  AUTO-GENERATING from occurrence table — P=R=F1=1.0 trivially.\n"
            "  This is NOT a real benchmark. Place a curated gold standard\n"
            "  in input/ and set gold_standard_relationships in config."
        )
        gold_df    = _auto_gold_standard(occ, cfg)
        using_auto = True
        logger.info(f"Auto-generated {len(gold_df)} gold-standard relationships")
    else:
        gold_df = load_table(gold_path)
        logger.info(f"Gold standard: {len(gold_df)} rows from {gold_path}")
        # FIX 1 — normalise relationship_type
        if "relationship_type" in gold_df.columns:
            gold_df["relationship_type"] = (gold_df["relationship_type"]
                                             .str.strip()
                                             .map(lambda x: _REL_ALIAS.get(x, x)))
        # FIX 2 — normalise expected_value
        if "expected_value" in gold_df.columns:
            gold_df["expected_value"] = gold_df["expected_value"].apply(_norm_val)
        logger.info(f"  Types: {gold_df['relationship_type'].value_counts().to_dict()}")

    preds = _extract_predictions(occ, cfg)
    logger.info("Predicted: " + str({k: len(v) for k, v in preds.items()}))

    results = run_benchmark(gold_df, preds)
    save_tsv(results["matched"],    out_dir/"benchmark_matched.tsv",    logger)
    save_tsv(results["missed"],     out_dir/"benchmark_missed.tsv",     logger)
    save_tsv(results["unexpected"], out_dir/"benchmark_unexpected.tsv", logger)
    save_tsv(results["summary"],    out_dir/"benchmark_summary.tsv",    logger)
    _benchmark_plot(results["summary"], out_dir/"benchmark_summary_plot.png",
                    plot_cfg.get("figure_dpi",150), plot_cfg.get("font_size",9),
                    matched_df=results["matched"], missed_df=results["missed"],
                    unexpected_df=results["unexpected"])

    if not results["summary"].empty:
        logger.info("")
        for _, row in results["summary"].iterrows():
            logger.info(f"  [{row['relationship_type']:<25}]  "
                        f"P={row['precision']:.3f}  R={row['recall']:.3f}  "
                        f"F1={row['f1']:.3f}  "
                        f"(gold={row['n_gold']}, pred={row['n_predicted']}, TP={row['n_true_pos']})")
        if using_auto:
            logger.warning("  Auto-generated gold standard used. Scores are trivially 1.0.")
        else:
            avg_r = results["summary"]["recall"].mean()
            avg_p = results["summary"]["precision"].mean()
            logger.info(f"\n  Average Recall={avg_r:.3f}  Precision={avg_p:.3f}")
            if avg_r >= 0.8:
                logger.info("  High recall — ontology covers curated biology")
            logger.info("  Low precision is expected (ontology predicts beyond gold standard)")
    logger.info("Stage 12 complete.")


if __name__ == "__main__":
    import sys, os
    _r = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(_r)
    _c = sys.argv[1] if len(sys.argv) > 1 else "config/pipeline_config.yaml"
    run(load_config(_c))
