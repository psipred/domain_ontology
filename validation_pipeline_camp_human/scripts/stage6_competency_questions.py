"""
scripts/stage6_competency_questions.py
────────────────────────────────────────────────────────────────────────────
STAGE 6 — COMPETENCY QUESTION ENGINE
Executes 5 built-in competency questions against the tabular occurrence data.
If rdflib is available, also runs SPARQL queries against instance_data.ttl.
Works fully offline without rdflib using pandas over the occurrence TSV.
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations
import sys as _sys
import os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from scripts.utils import (get_input_path, load_config, load_table,
                            make_output_dir, rename_to_canonical, save_tsv,
                            setup_logging, split_cell, stage_banner)
import pandas as pd


def _load_occ(cfg: dict) -> pd.DataFrame | None:
    path = get_input_path(cfg, "domain_occurrence_table")
    if not path:
        return None
    df = load_table(path)
    return rename_to_canonical(df, cfg.get("column_aliases", {}))


# ── Built-in competency questions (pandas) ────────────────────────────────────

def _cq01(occ: pd.DataFrame) -> pd.DataFrame:
    """CQ01: For each protein, which domain types are present and what are their function tags?"""
    rows = []
    for _, r in occ.iterrows():
        prot = str(r.get("protein_accession", "")).strip()
        dom  = str(r.get("domain_type_id", "")).strip()
        for ft in split_cell(r.get("function_tags", "")):
            rows.append({"protein": prot, "domain_type": dom, "function_tag": ft})
    df = pd.DataFrame(rows).drop_duplicates().sort_values(["protein", "domain_type"])
    return df


def _cq02(occ: pd.DataFrame) -> pd.DataFrame:
    """CQ02: Which pairs of domain types co-occur in ≥2 proteins?"""
    dom_col, prot_col = "domain_type_id", "protein_accession"
    # protein → set of domains
    p2d: dict[str, set] = {}
    for _, r in occ.iterrows():
        p = str(r.get(prot_col, "")).strip()
        for d in split_cell(r.get(dom_col, "")):
            p2d.setdefault(p, set()).add(d)
    # count co-occurrence pairs
    from collections import Counter
    pair_count: Counter = Counter()
    for domains in p2d.values():
        dl = sorted(domains)
        for i in range(len(dl)):
            for j in range(i + 1, len(dl)):
                pair_count[(dl[i], dl[j])] += 1
    rows = [{"domain_a": a, "domain_b": b, "protein_count": n}
            for (a, b), n in pair_count.items() if n >= 2]
    df = pd.DataFrame(rows).sort_values("protein_count", ascending=False).reset_index(drop=True)
    return df


def _cq03(occ: pd.DataFrame) -> pd.DataFrame:
    """CQ03: Which Reactome pathways are associated with each domain type, and how many proteins link them?"""
    rows = []
    for _, r in occ.iterrows():
        dom  = str(r.get("domain_type_id", "")).strip()
        prot = str(r.get("protein_accession", "")).strip()
        for rp in split_cell(r.get("protein_reactome_ids", "")):
            rows.append({"reactome_pathway": rp, "domain_type": dom, "protein": prot})
    if not rows:
        return pd.DataFrame(columns=["reactome_pathway", "domain_type", "protein_count"])
    df = (pd.DataFrame(rows).drop_duplicates()
          .groupby(["reactome_pathway", "domain_type"])
          .agg(protein_count=("protein", "nunique"))
          .reset_index()
          .sort_values("protein_count", ascending=False))
    return df


def _cq04(occ: pd.DataFrame) -> pd.DataFrame:
    """CQ04: Which ChEBI binding targets are linked to which proteins and domain types?"""
    chebi_col = next((c for c in ["binding_target_chebi", "binding_target_chebi_ids",
                                   "chebi_ids", "chebi_binding_ids"]
                      if c in occ.columns), None)
    rows = []
    if chebi_col:
        for _, r in occ.iterrows():
            prot = str(r.get("protein_accession", "")).strip()
            dom  = str(r.get("domain_type_id", "")).strip()
            for chebi in split_cell(r.get(chebi_col, "")):
                if chebi.strip():
                    rows.append({"binding_target_chebi": chebi.strip(),
                                 "protein": prot, "domain_type": dom})
    df = pd.DataFrame(rows).drop_duplicates().sort_values("binding_target_chebi") \
        if rows else pd.DataFrame(columns=["binding_target_chebi", "protein", "domain_type"])
    return df


def _cq05(occ: pd.DataFrame) -> pd.DataFrame:
    """CQ05: Distribution of domain occurrences across positional categories."""
    rows = []
    for v in occ.get("positional_category", pd.Series(dtype=str)).dropna():
        for cat in split_cell(v):
            if cat.strip():
                rows.append(cat.strip())
    if not rows:
        return pd.DataFrame(columns=["positional_category", "count"])
    df = (pd.Series(rows).value_counts()
          .reset_index()
          .rename(columns={"index": "positional_category", 0: "count",
                            "count": "positional_category", "proportion": "count"}))
    # pandas 2.x compat
    df.columns = ["positional_category", "count"]
    return df.sort_values("count", ascending=False)


BUILTIN_QUERIES = {
    "CQ01_protein_domain_function":       _cq01,
    "CQ02_domain_cooccurrence":           _cq02,
    "CQ03_process_linked_domains":        _cq03,
    "CQ04_binding_target_proteins":       _cq04,
    "CQ05_positional_category_distribution": _cq05,
}


def run(cfg: dict) -> None:
    logger = setup_logging(cfg.get("log_dir", "logs"), "stage6_competency")
    stage_banner(logger, "STAGE 6", "Competency Question Engine")
    out_dir = make_output_dir(cfg, "stage6_competency")

    # Try rdflib/SPARQL first (bonus path)
    rdflib_ok = False
    try:
        from rdflib import Graph
        inst_path = get_input_path(cfg, "instance_data")
        if inst_path:
            g = Graph()
            g.parse(str(inst_path))
            logger.info(f"rdflib available — RDF graph loaded: {len(g):,} triples")
            rdflib_ok = True
    except ImportError:
        logger.info("rdflib not installed — running pandas-based competency queries")
    except Exception as exc:
        logger.warning(f"rdflib failed ({exc}) — running pandas-based queries")

    # Always run pandas queries (these are the primary deliverable)
    occ = _load_occ(cfg)
    if occ is None:
        logger.error("domain_occurrence_table not found — cannot run competency questions")
        return

    logger.info(f"Occurrence table loaded: {len(occ)} rows")
    summary_rows = []

    for qname, fn in BUILTIN_QUERIES.items():
        logger.info(f"Running {qname} …")
        try:
            df = fn(occ)
            n = len(df)
            save_tsv(df, out_dir / f"{qname}.tsv", logger)
            status = "OK" if n > 0 else "EMPTY"
            summary_rows.append({"query_id": qname, "result_rows": n,
                                  "status": status, "engine": "pandas"})
            if n == 0:
                logger.warning(f"  {qname} → 0 rows — check data coverage")
            else:
                logger.info(f"  {qname} → {n} rows")
        except Exception as exc:
            logger.error(f"  {qname} failed: {exc}")
            summary_rows.append({"query_id": qname, "result_rows": 0,
                                  "status": f"ERROR: {exc}", "engine": "pandas"})

    # Also load and run any external .rq files if rdflib available
    if rdflib_ok:
        from pathlib import Path
        query_dir = Path(cfg.get("query_dir", "queries"))
        if query_dir.exists():
            for rq_file in sorted(query_dir.glob("*.rq")):
                qname = rq_file.stem
                logger.info(f"Running external SPARQL query: {qname}")
                try:
                    results = g.query(rq_file.read_text(encoding="utf-8"))
                    rows = [{str(var): str(val)
                             for var, val in zip(results.vars, row)}
                            for row in results]
                    df = pd.DataFrame(rows)
                    n = len(df)
                    save_tsv(df, out_dir / f"{qname}.tsv", logger)
                    summary_rows.append({"query_id": qname, "result_rows": n,
                                          "status": "OK" if n > 0 else "EMPTY",
                                          "engine": "sparql"})
                    logger.info(f"  {qname} → {n} rows")
                except Exception as exc:
                    logger.error(f"  {qname} failed: {exc}")
                    summary_rows.append({"query_id": qname, "result_rows": 0,
                                          "status": f"ERROR: {exc}",
                                          "engine": "sparql"})

    summary_df = pd.DataFrame(summary_rows)
    save_tsv(summary_df, out_dir / "competency_question_summary.tsv", logger)

    n_ok    = sum(1 for r in summary_rows if r["status"] == "OK")
    n_empty = sum(1 for r in summary_rows if r["status"] == "EMPTY")
    n_err   = sum(1 for r in summary_rows if r["status"] not in ("OK", "EMPTY"))
    logger.info(f"Stage 6 complete — {len(summary_rows)} queries: "
                f"{n_ok} OK, {n_empty} empty, {n_err} errors")


if __name__ == "__main__":
    import sys, os
    _project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(_project_root)
    _cfg_path = sys.argv[1] if len(sys.argv) > 1 else "config/pipeline_config.yaml"
    run(load_config(_cfg_path))
