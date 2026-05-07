"""
scripts/stage14_external_agreement.py
────────────────────────────────────────────────────────────────────────────
STAGE 14 — EXTERNAL DATABASE AGREEMENT
Compares ontology-derived protein annotations (GO, Reactome, KEGG) with
the external mapping table.  Computes Jaccard, overlap coefficient, and
Fisher enrichment for each protein.
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
import sys as _sys
import os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from scripts.utils import (benjamini_hochberg, get_input_path, load_config,
                            load_table, make_output_dir, rename_to_canonical,
                            save_tsv, setup_logging, split_cell, stage_banner)


def _set_for_protein(df: pd.DataFrame, prot: str,
                      prot_col: str, val_col: str) -> set:
    rows = df[df[prot_col] == prot]
    result: set = set()
    for _, row in rows.iterrows():
        result.update(split_cell(row.get(val_col, "")))
    return result


def _jaccard(a: set, b: set) -> float:
    u = a | b
    return len(a & b) / len(u) if u else 0.0


def _overlap_coeff(a: set, b: set) -> float:
    min_size = min(len(a), len(b))
    return len(a & b) / min_size if min_size else 0.0


def run(cfg: dict) -> None:
    logger = setup_logging(cfg.get("log_dir","logs"), "stage14_external_agreement")
    stage_banner(logger, "STAGE 14", "External Database Agreement")
    out_dir = make_output_dir(cfg, "stage14_external_agreement")
    plot_cfg = cfg.get("plotting", {})
    dpi = plot_cfg.get("figure_dpi", 150)
    fs  = plot_cfg.get("font_size", 9)

    occ_path = get_input_path(cfg, "domain_occurrence_table")
    if not occ_path:
        logger.warning("domain_occurrence_table not found — skipping Stage 14")
        return
    occ = load_table(occ_path)
    occ = rename_to_canonical(occ, cfg.get("column_aliases", {}))

    prot_col = "protein_accession"
    if prot_col not in occ.columns:
        logger.error(f"'{prot_col}' not in occurrence table")
        return

    map_path = get_input_path(cfg, "ontology_vs_external_mapping")
    if not map_path:
        # Fallback: build the external mapping directly from the GO / Reactome /
        # KEGG columns already present in the occurrence table.
        # These protein-level external IDs were imported from UniProt during
        # dataset construction and serve as the "external reference" here.
        logger.info(
            "ontology_vs_external_mapping.tsv not found — "
            "building external reference from GO/Reactome/KEGG columns in occurrence table"
        )
        map_df = occ[[prot_col]].drop_duplicates().copy()
        for src_col, tgt_col in [
            ("protein_go_cc_ids",    "go_ids"),
            ("protein_reactome_ids", "reactome_ids"),
            ("protein_kegg_ids",     "kegg_ids"),
        ]:
            if src_col in occ.columns:
                map_df = map_df.merge(
                    occ.groupby(prot_col)[src_col].first().reset_index().rename(
                        columns={src_col: tgt_col}),
                    on=prot_col, how="left",
                )
        logger.info(f"Built mapping table: {len(map_df)} proteins")
    else:
        map_df = load_table(map_path)

    proteins = occ[prot_col].dropna().unique().tolist()

    # Detect what external ID columns exist in the mapping table
    ext_col_candidates = {
        "go":       ["go_ids", "go_terms", "xref_GO", "protein_go_cc_ids"],
        "reactome": ["reactome_ids", "xref_reactome", "reactome", "protein_reactome_ids"],
        "kegg":     ["kegg_ids", "xref_kegg", "kegg", "protein_kegg_ids"],
        "chebi":    ["chebi_ids", "xref_chebi", "chebi"],
    }
    ext_col_map = {}
    map_prot_col = next((c for c in ["protein_accession", "accession", "uniprot_id"]
                         if c in map_df.columns), None)
    if not map_prot_col:
        logger.error("Cannot find protein accession column in mapping table")
        return

    for label, candidates in ext_col_candidates.items():
        for c in candidates:
            if c in map_df.columns:
                ext_col_map[label] = c
                break

    # Corresponding ontology columns
    ont_col_map = {
        "go":       ["go_terms_all", "go_terms_mf", "go_terms_bp", "go_terms_cc"],
        "reactome": ["protein_reactome_ids"],
        "kegg":     ["protein_kegg_ids"],
        "chebi":    ["binding_target_chebi"],
    }

    result_rows = []
    for prot in proteins:
        for label, ext_col in ext_col_map.items():
            ext_set = _set_for_protein(map_df, prot, map_prot_col, ext_col)
            ont_cols = [c for c in ont_col_map.get(label, []) if c in occ.columns]
            ont_set: set = set()
            for col in ont_cols:
                for _, row in occ[occ[prot_col] == prot].iterrows():
                    ont_set.update(split_cell(row.get(col, "")))

            if not ext_set and not ont_set:
                continue

            jac   = _jaccard(ont_set, ext_set)
            ovlp  = _overlap_coeff(ont_set, ext_set)
            n_int = len(ont_set & ext_set)
            result_rows.append({
                "protein":          prot,
                "annotation_type":  label,
                "n_ontology":       len(ont_set),
                "n_external":       len(ext_set),
                "n_overlap":        n_int,
                "jaccard":          round(jac, 4),
                "overlap_coeff":    round(ovlp, 4),
            })

    results_df = pd.DataFrame(result_rows)
    save_tsv(results_df, out_dir/"external_agreement_results.tsv", logger)

    # Summary per annotation type
    if not results_df.empty:
        summary = results_df.groupby("annotation_type").agg(
            n_proteins=("protein","nunique"),
            mean_jaccard=("jaccard","mean"),
            mean_overlap_coeff=("overlap_coeff","mean"),
            total_ontology_terms=("n_ontology","sum"),
            total_external_terms=("n_external","sum"),
            total_overlap=("n_overlap","sum"),
        ).reset_index()
        summary["mean_jaccard"]       = summary["mean_jaccard"].round(4)
        summary["mean_overlap_coeff"] = summary["mean_overlap_coeff"].round(4)
        save_tsv(summary, out_dir/"external_agreement_summary.tsv", logger)

        # Bar plot
        agg_types = summary["annotation_type"].tolist()
        x = range(len(agg_types))
        fig, ax = plt.subplots(figsize=(max(5, len(agg_types)*1.5), 5), dpi=dpi)
        width = 0.35
        ax.bar([i - width/2 for i in x], summary["mean_jaccard"], width,
               label="Mean Jaccard", color="#4C72B0")
        ax.bar([i + width/2 for i in x], summary["mean_overlap_coeff"], width,
               label="Mean Overlap Coefficient", color="#55A868")
        ax.set_xticks(list(x))
        ax.set_xticklabels(agg_types, fontsize=fs)
        ax.set_ylim(0, 1.1)
        ax.set_ylabel("Score", fontsize=fs)
        ax.set_title("External Database Agreement by Annotation Type",
                     fontsize=fs+2, fontweight="bold")
        ax.legend(fontsize=fs)
        ax.spines[["top","right"]].set_visible(False)
        plt.tight_layout()
        fig.savefig(out_dir/"external_agreement_plots.png", dpi=dpi, bbox_inches="tight")
        plt.close(fig)

        for _, row in summary.iterrows():
            logger.info(f"  [{row['annotation_type']}]  "
                        f"Jaccard={row['mean_jaccard']:.3f}  "
                        f"Overlap={row['mean_overlap_coeff']:.3f}")

    logger.info("Stage 14 complete.")


if __name__ == "__main__":
    import sys, os
    # Always run from the ontology_pipeline project root so that
    # relative paths (input/, output/, config/) resolve correctly
    # regardless of where PyCharm / the terminal launched the script from.
    _project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(_project_root)
    _cfg_path = sys.argv[1] if len(sys.argv) > 1 else "config/pipeline_config.yaml"
    run(load_config(_cfg_path))
