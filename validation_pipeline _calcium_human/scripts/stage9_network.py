"""
scripts/stage9_network.py
────────────────────────────────────────────────────────────────────────────
STAGE 9 — ONTOLOGY-DERIVED NETWORK
Builds a heterogeneous graph of proteins, domains, functions, binding targets,
and processes using networkx.  Exports node/edge tables and a figure.
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

try:
    import networkx as nx
    _HAS_NX = True
except ImportError:
    _HAS_NX = False

import sys as _sys
import os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from scripts.utils import (get_input_path, load_config, load_table,
                            make_output_dir, rename_to_canonical,
                            save_tsv, setup_logging, split_cell, stage_banner)

# Node type colour palette (matplotlib-friendly, no seaborn)
NODE_COLORS = {
    "protein":       "#4C72B0",
    "domain_type":   "#DD8452",
    "function_tag":  "#55A868",
    "binding_target":"#C44E52",
    "process":       "#8172B2",
}


def _add_nodes(G, node_id: str, node_type: str,
               label: str = "") -> None:
    if not G.has_node(node_id):
        G.add_node(node_id, node_type=node_type,
                   label=label or node_id)


def build_network(occ: pd.DataFrame, cfg: dict) -> "nx.Graph":
    G = nx.Graph()
    aliases = cfg.get("column_aliases", {})
    occ = rename_to_canonical(occ, aliases)

    prot_col   = "protein_accession"
    dom_col    = "domain_type_id"
    func_col   = "function_tags"
    # After rename_to_canonical the typo is corrected — check both forms.
    bind_col   = next(
        (c for c in ["binding_partner", "binding_parnter"] if c in occ.columns),
        "binding_partner",
    )
    proc_col   = "process"

    for _, row in occ.iterrows():
        # Protein node
        prot = str(row.get(prot_col, "")).strip()
        if not prot or prot == "nan":
            continue
        _add_nodes(G, prot, "protein")

        # Domain type nodes and protein–domain edges
        for dt in split_cell(row.get(dom_col, "")):
            _add_nodes(G, dt, "domain_type")
            G.add_edge(prot, dt, edge_type="has_domain")

            # Domain–function edges
            for ft in split_cell(row.get(func_col, "")):
                _add_nodes(G, ft, "function_tag")
                G.add_edge(dt, ft, edge_type="has_function")

            # Domain–binding edges
            for bt in split_cell(row.get(bind_col, "")):
                _add_nodes(G, bt, "binding_target")
                G.add_edge(dt, bt, edge_type="has_binding_target")

            # Protein–process edges
            if proc_col in occ.columns:
                for proc in split_cell(row.get(proc_col, "")):
                    _add_nodes(G, proc, "process")
                    G.add_edge(prot, proc, edge_type="involved_in_process")

    return G


def _draw_network(G: "nx.Graph", out_path: Path,
                  min_degree: int = 2, dpi: int = 150,
                  font_size: int = 7, cmap_name: str = "tab10",
                  logger=None) -> None:
    """Produce two PNG files: full network + top-hub subgraph."""
    import logging
    _log = logger or logging.getLogger("stage9_network")

    edge_style = {"has_domain":           ("#AAAAAA", 0.20, 0.5),
                  "has_function":          ("#55A868", 0.55, 1.0),
                  "has_binding_target":    ("#C44E52", 0.55, 1.0),
                  "involved_in_process":   ("#8172B2", 0.45, 0.8)}

    def _render(H, ax, pos, node_type_list):
        nc = [NODE_COLORS.get(t, "#AAAAAA") for t in node_type_list]
        ns = [max(40, H.degree(n) * 18) for n in H.nodes]
        for etype, (ec, alpha, width) in edge_style.items():
            elist = [(u, v) for u, v, d in H.edges(data=True)
                     if d.get("edge_type") == etype]
            if elist:
                nx.draw_networkx_edges(H, pos, edgelist=elist,
                                       edge_color=ec, alpha=alpha,
                                       width=width, ax=ax)
        nx.draw_networkx_nodes(H, pos, node_color=nc, node_size=ns,
                               alpha=0.88, linewidths=0.4,
                               edgecolors="white", ax=ax)
        top_n  = max(20, len(H.nodes) // 6)
        top_nd = sorted(H.degree, key=lambda x: x[1], reverse=True)[:top_n]
        labels = {n: H.nodes[n].get("label", n)
                     .replace("INTERPRO:", "").replace("PFAM:", "")
                     .replace("CDD:", "").replace("SMART:", "")
                  for n, _ in top_nd}
        nx.draw_networkx_labels(H, pos, labels=labels,
                                font_size=max(5, font_size - 1), ax=ax)
        handles = [mpatches.Patch(facecolor=c,
                                  label=t.replace("_", " ").title())
                   for t, c in NODE_COLORS.items()
                   if t in set(node_type_list)]
        ax.legend(handles=handles, loc="upper left", fontsize=font_size,
                  framealpha=0.85, title="Node type")
        ax.axis("off")

    # ── Plot 1: Full network (spring layout — fast, always works) ────────────
    nodes_keep = [n for n, d in G.degree() if d >= min_degree]
    H = G.subgraph(nodes_keep).copy()
    if len(H.nodes) == 0:
        _log.warning("No nodes with degree >= %d — network plot skipped", min_degree)
        return

    node_type_list = [H.nodes[n].get("node_type", "domain_type") for n in H.nodes]
    _log.info("Generating full network layout (%d nodes) …", len(H.nodes))
    k_val = 1.8 / (len(H.nodes) ** 0.5)
    pos_full = nx.spring_layout(H, seed=42, k=k_val, iterations=60)

    fig, ax = plt.subplots(figsize=(18, 14), dpi=dpi)
    _render(H, ax, pos_full, node_type_list)
    ax.set_title(f"cAMP Pathway Ontology Network — Full Graph\n"
                 f"({len(H.nodes)} nodes · {len(H.edges)} edges · "
                 f"min_degree≥{min_degree})",
                 fontsize=font_size + 3, fontweight="bold", pad=12)
    plt.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    _log.info("Saved full network → %s", out_path.name)

    # ── Plot 2: Top-hub subgraph (top 60 nodes by degree, Kamada-Kawai) ─────
    hub_path = out_path.parent / "ontology_network_hubs.png"
    top60 = [n for n, _ in sorted(H.degree, key=lambda x: -x[1])[:60]]
    HH = H.subgraph(top60).copy()
    nt_hubs = [HH.nodes[n].get("node_type", "domain_type") for n in HH.nodes]
    _log.info("Generating hub subgraph layout (%d nodes) …", len(HH.nodes))
    try:
        pos_hub = nx.kamada_kawai_layout(HH)
    except Exception:
        pos_hub = nx.spring_layout(HH, seed=42, k=0.8)

    fig2, ax2 = plt.subplots(figsize=(16, 13), dpi=dpi)
    _render(HH, ax2, pos_hub, nt_hubs)
    ax2.set_title(f"cAMP Pathway Ontology Network — Top-60 Hub Nodes\n"
                  f"(Kamada-Kawai layout · {len(HH.nodes)} nodes · "
                  f"{len(HH.edges)} edges)",
                  fontsize=font_size + 3, fontweight="bold", pad=12)
    plt.tight_layout()
    fig2.savefig(hub_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig2)
    _log.info("Saved hub network → %s", hub_path.name)


def run(cfg: dict) -> None:
    logger = setup_logging(cfg.get("log_dir","logs"), "stage9_network")
    stage_banner(logger, "STAGE 9", "Ontology-Derived Network")
    out_dir = make_output_dir(cfg, "stage9_network")

    if not _HAS_NX:
        logger.error("networkx not installed — pip install networkx")
        return

    occ_path = get_input_path(cfg, "domain_occurrence_table")
    if not occ_path:
        logger.warning("domain_occurrence_table not found — skipping Stage 9")
        return

    occ = load_table(occ_path)
    plot_cfg = cfg.get("plotting", {})
    min_deg  = plot_cfg.get("network_min_degree", 2)
    dpi      = plot_cfg.get("figure_dpi", 150)

    G = build_network(occ, cfg)
    logger.info(f"Network: {G.number_of_nodes()} nodes, "
                f"{G.number_of_edges()} edges")

    # Export node and edge tables
    node_rows = [{"node_id": n, **G.nodes[n]} for n in G.nodes]
    edge_rows = [{"source": u, "target": v, **d}
                 for u, v, d in G.edges(data=True)]
    save_tsv(pd.DataFrame(node_rows), out_dir / "network_nodes.tsv", logger)
    save_tsv(pd.DataFrame(edge_rows), out_dir / "network_edges.tsv", logger)

    logger.info("Generating network plots …")
    _draw_network(G, out_dir / "ontology_network.png", min_deg, dpi,
                  plot_cfg.get("font_size", 7), logger=logger)
    logger.info("Stage 9 complete.")


if __name__ == "__main__":
    import sys, os
    # Always run from the ontology_pipeline project root so that
    # relative paths (input/, output/, config/) resolve correctly
    # regardless of where PyCharm / the terminal launched the script from.
    _project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(_project_root)
    _cfg_path = sys.argv[1] if len(sys.argv) > 1 else "config/pipeline_config.yaml"
    run(load_config(_cfg_path))
