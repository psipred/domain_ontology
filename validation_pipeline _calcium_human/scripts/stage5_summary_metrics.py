"""
scripts/stage5_summary_metrics.py
────────────────────────────────────────────────────────────────────────────
STAGE 5 — ONTOLOGY SUMMARY METRICS
Counts classes, properties, individuals, mapped external IDs, and rule stats.
Produces a publication-style bar plot using matplotlib only.
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations
from pathlib import Path
from collections import defaultdict
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys as _sys
import os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from scripts.utils import (get_input_path, load_config, load_table,
                            make_output_dir, save_tsv, setup_logging,
                            split_cell, stage_banner)


def _count_rdf_entities(path: Path) -> dict:
    """Count OWL classes, object properties, data properties via rdflib."""
    try:
        from rdflib import Graph, RDF, OWL
        g = Graph().parse(str(path))
        counts = {
            "owl_classes": sum(1 for _ in g.subjects(RDF.type, OWL.Class)),
            "object_properties": sum(1 for _ in g.subjects(RDF.type, OWL.ObjectProperty)),
            "data_properties": sum(1 for _ in g.subjects(RDF.type, OWL.DatatypeProperty)),
        }
        return counts
    except Exception:
        return {"owl_classes": "N/A", "object_properties": "N/A",
                "data_properties": "N/A"}


def _count_external_ids(df: pd.DataFrame, col_patterns: dict) -> dict:
    import re
    counts = {}
    for label, (col, pattern) in col_patterns.items():
        if col not in df.columns:
            counts[label] = 0
            continue
        rx = re.compile(pattern)
        seen: set = set()
        for val in df[col].dropna():
            for tok in split_cell(val):
                if rx.fullmatch(tok):
                    seen.add(tok)
        counts[label] = len(seen)
    return counts


def run(cfg: dict) -> None:
    logger = setup_logging(cfg.get("log_dir","logs"), "stage5_summary_metrics")
    stage_banner(logger, "STAGE 5", "Ontology Summary Metrics")
    out_dir = make_output_dir(cfg, "stage5_metrics")
    plot_cfg = cfg.get("plotting", {})
    dpi = plot_cfg.get("figure_dpi", 150)

    metrics: dict[str, object] = {}

    # OWL file metrics
    ont_path = get_input_path(cfg, "core_ontology")
    if ont_path:
        metrics.update(_count_rdf_entities(ont_path))
    else:
        metrics.update({"owl_classes": "N/A", "object_properties": "N/A",
                        "data_properties": "N/A"})

    # RDF instance metrics
    inst_path = get_input_path(cfg, "instance_data")
    if inst_path:
        try:
            from rdflib import Graph, RDF
            g = Graph().parse(str(inst_path))
            types: dict = defaultdict(int)
            for _, _, obj in g.triples((None, RDF.type, None)):
                local = str(obj).split("#")[-1] if "#" in str(obj) else str(obj)
                types[local] += 1
            for cls in ["Protein","DomainOccurrence","DomainType",
                        "GOTerm","FunctionalCategoryContext",
                        "BindingTarget","PositionalCategory"]:
                metrics[f"rdf_{cls}"] = types.get(cls, 0)
        except Exception as e:
            logger.warning(f"RDF instance counting failed: {e}")

    # Table-based metrics
    dt_path = get_input_path(cfg, "domain_type_summary")
    occ_path = get_input_path(cfg, "domain_occurrence_table")
    dt  = load_table(dt_path) if dt_path else pd.DataFrame()
    occ = load_table(occ_path) if occ_path else pd.DataFrame()

    metrics["domain_types"] = len(dt)
    metrics["domain_occurrences"] = len(occ)
    if "protein_accession" in occ.columns:
        metrics["proteins"] = occ["protein_accession"].nunique()

    ext_patterns = {
        "go_terms":        ("go_terms_all",    r"GO:\d{7}"),
        "binding_targets": ("binding_partner", r"DOT:\d+(?:\s+.+)?"),
    }
    metrics.update(_count_external_ids(occ, ext_patterns))

    # Validation rule stats from stage 2 output (if exists)
    rule_summary = out_dir.parent / "stage2_conformance" / "rule_violation_summary.tsv"
    if rule_summary.exists():
        rs = pd.read_csv(rule_summary, sep="\t")
        metrics["n_rules"] = len(rs)
        metrics["rules_failing"] = int((rs["status"] == "FAIL").sum())
        metrics["rules_passing"] = int((rs["status"] == "PASS").sum())
    else:
        metrics["n_rules"] = "N/A"

    # Save TSV
    summary_df = pd.DataFrame([
        {"metric": k, "value": v} for k, v in metrics.items()
    ])
    save_tsv(summary_df, out_dir / "ontology_summary_metrics.tsv", logger)

    # Bar plot
    numeric = {k: v for k, v in metrics.items() if isinstance(v, (int, float))}
    if numeric:
        fig, ax = plt.subplots(figsize=(10, max(4, len(numeric) * 0.45)), dpi=dpi)
        labels = list(numeric.keys())
        values = [numeric[k] for k in labels]
        y_pos = range(len(labels))
        bars = ax.barh(list(y_pos), values, color="#4C72B0", edgecolor="white",
                       height=0.65)
        ax.set_yticks(list(y_pos))
        ax.set_yticklabels(labels, fontsize=plot_cfg.get("font_size", 9))
        ax.set_xlabel("Count", fontsize=plot_cfg.get("font_size", 9))
        ax.set_title("Ontology Summary Metrics",
                     fontsize=plot_cfg.get("title_font_size", 11), fontweight="bold")
        for bar, val in zip(bars, values):
            ax.text(bar.get_width() + max(values) * 0.01,
                    bar.get_y() + bar.get_height() / 2,
                    str(val), va="center",
                    fontsize=plot_cfg.get("font_size", 9) - 1)
        ax.spines[["top", "right"]].set_visible(False)
        plt.tight_layout()
        fig_path = out_dir / "ontology_summary_barplot.png"
        fig.savefig(fig_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        logger.info(f"Bar plot saved → {fig_path}")

    logger.info("Stage 5 complete.")


if __name__ == "__main__":
    import sys, os
    # Always run from the ontology_pipeline project root so that
    # relative paths (input/, output/, config/) resolve correctly
    # regardless of where PyCharm / the terminal launched the script from.
    _project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(_project_root)
    _cfg_path = sys.argv[1] if len(sys.argv) > 1 else "config/pipeline_config.yaml"
    run(load_config(_cfg_path))
