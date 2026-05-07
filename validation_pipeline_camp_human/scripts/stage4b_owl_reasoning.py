"""
scripts/stage4b_owl_reasoning.py
────────────────────────────────────────────────────────────────────────────
STAGE 4b — OWL REASONING  (pure Python, no Java required)

Replaces the Java-based stage4_owl_reasoning.  Uses two Python libraries:

  owlrl   — OWL RL reasoning over rdflib graphs (rule-based, no JVM)
  rdflib  — SPARQL queries after reasoning closes the graph

What it does
────────────
  1. Loads  core_ontology.ttl  +  population.ttl  into one rdflib Graph
  2. Runs   OWL RL closure  via owlrl.DeductiveClosure
             • Infers: rdfs:subClassOf chains, owl:inverseOf,
               owl:TransitiveProperty, rdfs:domain/range propagation
             • Detects: anything inferred to be owl:Nothing → inconsistency
  3. Reports
             • Is the ontology consistent?
             • How many new triples were inferred?
             • Which classes/properties triggered inferences?
             • Any unsatisfiable individuals (rdf:type owl:Nothing)?
  4. Saves
             • inferred_triples.ttl          — new triples only
             • consistency_report.tsv        — pass/fail + counts
             • reasoning_summary.tsv         — per-rule counts
  5. Runs    DL-equivalent SPARQL queries (translated from Manchester syntax)
             • Same queries you would type in Protege DL Query tab
             • Results saved as TSV for each query

Install (once)
──────────────
  pip install owlrl rdflib

Config keys  (pipeline_config.yaml)
─────────────────────────────────────
  core_ontology      — core_ontology.ttl  (required)
  instance_data      — population.ttl     (required, output of stage 4)
  owl_dl_queries     — path to dl_queries.yaml  (optional, see below)

DL query file format  (dl_queries.yaml)
────────────────────────────────────────
  - id: Q1
    label: "Domain types with catalytic function"
    sparql: |
      SELECT DISTINCT ?dt ?label WHERE {
        ?dt a core:DomainType ;
            core:hasFunctionalRole ?r .
        ?r rdfs:label ?label .
        FILTER(str(?label) = "catalytic")
      }
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional
import pandas as pd

import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from scripts.utils import (
    get_input_path, load_config, make_output_dir,
    save_tsv, setup_logging, stage_banner,
)


# ══════════════════════════════════════════════════════════════════════════════
# Built-in DL queries  (translated from Protege Manchester syntax to SPARQL)
# ══════════════════════════════════════════════════════════════════════════════

_BUILTIN_DL_QUERIES = [
    # ── DL01: Domain types with scaffold function (DOF:000016) ────────────────
    # Manchester: DomainType and hasFunctionalRole some {scaffold}
    # Queries DomainType-level hasFunctionalRole.  The FuncCat individual
    # pop:FuncCat_DOF_000016_scaffold carries rdfs:label = "scaffold".
    {
        "id": "DL01",
        "label": "Domain types with scaffold function (DOF:000016)",
        "description": "Protege: DomainType and hasFunctionalRole some {DOF:000016 scaffold}",
        "sparql": """
            SELECT DISTINCT ?dt ?dtLabel ?funcLabel WHERE {
              ?dt a <%(CORE)sDomainType> .
              ?dt <%(CORE)shasFunctionalRole> ?r .
              ?dt rdfs:label ?dtLabel .
              ?r  rdfs:label ?funcLabel .
              FILTER(LCASE(STR(?funcLabel)) = "scaffold")
            } ORDER BY ?dtLabel
        """,
    },
    # ── DL02: Domain types with receptor-binding function (DOF:000038) ───────
    # Manchester: DomainType and hasFunctionalRole some {receptor-binding}
    {
        "id": "DL02",
        "label": "Domain types with receptor-binding function (DOF:000038)",
        "description": "Protege: DomainType and hasFunctionalRole some {DOF:000038 receptor-binding}",
        "sparql": """
            SELECT DISTINCT ?dt ?dtLabel ?funcLabel WHERE {
              ?dt a <%(CORE)sDomainType> .
              ?dt <%(CORE)shasFunctionalRole> ?r .
              ?dt rdfs:label ?dtLabel .
              ?r  rdfs:label ?funcLabel .
              FILTER(LCASE(STR(?funcLabel)) = "receptor-binding")
            } ORDER BY ?dtLabel
        """,
    },
    # ── DL03: Transmembrane domain occurrences (DOP:000008) ───────────────────
    # Manchester: DomainOccurrence and hasCharacterisation
    #             some TopologyContext[DOP:000008 transmembrane]
    # Also returns occurrence-level hasFunctionalRole for cross-reference.
    {
        "id": "DL03",
        "label": "Transmembrane domain occurrences (DOP:000008)",
        "description": "Protege: DomainOccurrence and hasCharacterisation some "
                       "TopologyContext[DOP:000008 transmembrane]",
        "sparql": """
            SELECT DISTINCT ?occ ?prot ?dt ?dtLabel WHERE {
              ?occ a <%(CORE)sDomainOccurrence> .
              ?occ <%(CORE)soccursInProtein> ?prot .
              ?occ <%(CORE)shasDomainType>   ?dt .
              ?dt  rdfs:label ?dtLabel .
              ?occ <%(CORE)shasCharacterisation> ?tc .
              ?tc  a <%(CORE)sTopologyContext> .
              ?tc  rdfs:label ?tcLabel .
              FILTER(LCASE(STR(?tcLabel)) = "transmembrane")
            } ORDER BY ?prot ?dt
        """,
    },
    # ── DL04: Proteins whose scaffold domain occurs in a transmembrane context
    # Manchester: Protein and inverse(occursInProtein) some
    #   (DomainOccurrence
    #     and hasFunctionalRole some {DOF:000016 scaffold}      ← occ-level
    #     and hasCharacterisation some {DOP:000008 transmembrane})
    #
    # Key change vs. old DL04: hasFunctionalRole is now queried on the
    # DomainOccurrence directly (occurrence-level), not via hasDomainType →
    # DomainType → hasFunctionalRole.  This correctly handles the case where
    # the same domain family has different functions in different proteins.
    {
        "id": "DL04",
        "label": "Proteins with scaffold domain in transmembrane context",
        "description": "Protege: Protein and inverse(occursInProtein) some "
                       "(DomainOccurrence and hasFunctionalRole some {DOF:000016 scaffold} "
                       "and hasCharacterisation some {DOP:000008 transmembrane})",
        "sparql": """
            SELECT DISTINCT ?prot ?occ ?dt ?dtLabel WHERE {
              ?occ  a <%(CORE)sDomainOccurrence> .
              ?occ  <%(CORE)soccursInProtein>    ?prot .
              ?occ  <%(CORE)shasDomainType>      ?dt .
              ?dt   rdfs:label                   ?dtLabel .
              ?occ  <%(CORE)shasFunctionalRole>  ?fr .
              ?fr   rdfs:label                   ?frLabel .
              ?occ  <%(CORE)shasCharacterisation> ?tc .
              ?tc   rdfs:label                   ?tcLabel .
              FILTER(LCASE(STR(?frLabel))  = "scaffold")
              FILTER(LCASE(STR(?tcLabel))  = "transmembrane")
            } ORDER BY ?prot ?dt
        """,
    },
    {
        "id": "DL05",
        "label": "Co-occurring domain pairs",
        "description": "Protege: DomainType and cooccursWith some DomainType",
        "sparql": """
            SELECT DISTINCT ?dtA ?labelA ?dtB ?labelB WHERE {
              ?dtA a <%(CORE)sDomainType> .
              ?dtA <%(CORE)scooccursWith> ?dtB .
              ?dtA rdfs:label ?labelA .
              ?dtB rdfs:label ?labelB .
              FILTER(STR(?dtA) < STR(?dtB))
            } ORDER BY ?labelA ?labelB
        """,
    },
    {
        "id": "DL06",
        "label": "Proteins with multiple domain types",
        "description": "Proteins that have more than one distinct domain type",
        "sparql": """
            SELECT ?prot (COUNT(DISTINCT ?dt) AS ?nDomains) WHERE {
              ?occ a <%(CORE)sDomainOccurrence> .
              ?occ <%(CORE)soccursInProtein> ?prot .
              ?occ <%(CORE)shasDomainType> ?dt .
            }
            GROUP BY ?prot
            HAVING(COUNT(DISTINCT ?dt) > 1)
            ORDER BY DESC(?nDomains)
        """,
    },
    {
        "id": "DL07",
        "label": "Consistency check — unsatisfiable individuals",
        "description": "Any individual inferred to be rdf:type owl:Nothing indicates an inconsistency",
        "sparql": """
            SELECT ?ind ?type WHERE {
              ?ind a <http://www.w3.org/2002/07/owl#Nothing> .
              ?ind a ?type .
            }
        """,
    },
    {
        "id": "DL08",
        "label": "Domain types with GO annotations",
        "description": "Protege: DomainType and hasGOAnnotation some GOTerm",
        "sparql": """
            SELECT DISTINCT ?dt ?dtLabel ?go WHERE {
              ?dt a <%(CORE)sDomainType> .
              ?dt rdfs:label ?dtLabel .
              ?dt <%(CORE)shasGOAnnotation> ?go .
            } ORDER BY ?dtLabel
        """,
    },
    {
        "id": "DL09",
        "label": "Domain order chains",
        "description": "Domain types that have a documented sequential neighbour",
        "sparql": """
            SELECT DISTINCT ?occ ?prevDT ?thisDT ?nextDT WHERE {
              ?occ a <%(CORE)sDomainOccurrence> .
              ?occ <%(CORE)shasDomainType> ?thisDT .
              OPTIONAL { ?occ <%(CORE)sisDomainNeighbour> ?prevDT . }
              OPTIONAL { ?occ <%(CORE)sisDomainNeighbour> ?nextDT .
                         FILTER(?nextDT != ?prevDT) }
            } ORDER BY ?occ
        """,
    },
    {
        "id": "DL10",
        "label": "Binding target coverage",
        "description": "Occurrences that have ChEBI-linked binding targets",
        "sparql": """
            SELECT DISTINCT ?occ ?prot ?chebi WHERE {
              ?occ a <%(CORE)sDomainOccurrence> .
              ?occ <%(CORE)soccursInProtein> ?prot .
              ?occ <%(CORE)shasBindingTarget> ?chebi .
            } ORDER BY ?prot
        """,
    },
]


# ══════════════════════════════════════════════════════════════════════════════
# Reasoning
# ══════════════════════════════════════════════════════════════════════════════

def run_reasoning(
    core_path: str,
    pop_path: str,
    logger,
) -> tuple["Graph", "Graph", dict]:
    """
    Load core + population ontologies, run OWL RL closure, return
    (merged_graph, inferred_only_graph, stats_dict).
    """
    try:
        from rdflib import Graph, OWL, RDF
        import owlrl
    except ImportError as e:
        logger.error(f"Missing library: {e}")
        logger.error("Install with:  pip install owlrl rdflib")
        return None, None, {}

    logger.info("Loading ontology files …")
    g = Graph()
    n_before = 0

    for path, label in [(core_path, "core"), (pop_path, "population")]:
        p = Path(path)
        if not p.exists():
            logger.warning(f"  {label}: {path} not found — skipping")
            continue
        before = len(g)
        g.parse(str(p), format="turtle")
        added = len(g) - before
        logger.info(f"  Loaded {label}: {added:,} triples  ({path})")

    n_before = len(g)
    logger.info(f"Total before reasoning: {n_before:,} triples")

    # Run OWL RL closure
    logger.info("Running OWL RL reasoning (owlrl.DeductiveClosure) …")
    owlrl.DeductiveClosure(
        owlrl.OWLRL_Semantics,
        axiomatic_triples=False,
        datatype_axioms=False,
    ).expand(g)

    n_after  = len(g)
    n_inferred = n_after - n_before
    logger.info(f"After reasoning: {n_after:,} triples  (+{n_inferred:,} inferred)")

    # Extract inferred-only triples
    g_before = Graph()
    for path in [core_path, pop_path]:
        p = Path(path)
        if p.exists():
            g_before.parse(str(p), format="turtle")

    g_inferred = Graph()
    for triple in g:
        if triple not in g_before:
            g_inferred.add(triple)

    # Check consistency — anything of type owl:Nothing = inconsistency
    nothing = OWL.Nothing
    unsat = list(g.subjects(RDF.type, nothing))
    consistent = len(unsat) == 0

    stats = {
        "consistent":       consistent,
        "n_triples_before": n_before,
        "n_triples_after":  n_after,
        "n_inferred":       n_inferred,
        "n_unsatisfiable":  len(unsat),
        "unsatisfiable":    [str(u) for u in unsat],
    }

    if consistent:
        logger.info("  Consistency: PASS — no unsatisfiable individuals")
    else:
        logger.warning(f"  Consistency: FAIL — {len(unsat)} unsatisfiable individuals:")
        for u in unsat[:10]:
            logger.warning(f"    {u}")

    return g, g_inferred, stats


# ══════════════════════════════════════════════════════════════════════════════
# DL queries
# ══════════════════════════════════════════════════════════════════════════════

def run_dl_queries(
    g: "Graph",
    core_ns: str,
    out_dir: Path,
    extra_queries: list[dict],
    logger,
) -> pd.DataFrame:
    """
    Run built-in + user-defined DL queries against the reasoned graph.
    Returns a summary DataFrame.
    """
    from rdflib import Graph
    from rdflib.plugins.sparql import prepareQuery
    import re

    all_queries = _BUILTIN_DL_QUERIES + (extra_queries or [])
    summary_rows = []

    for q in all_queries:
        qid   = q["id"]
        label = q["label"]
        sparql_template = q["sparql"]
        sparql = sparql_template % {"CORE": core_ns}

        # Inject common prefixes
        prefixed = (
            "PREFIX rdf:  <http://www.w3.org/1999/02/22-rdf-syntax-ns#>\n"
            "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>\n"
            "PREFIX owl:  <http://www.w3.org/2002/07/owl#>\n"
            "PREFIX core: <" + core_ns + ">\n"
        ) + sparql

        try:
            result = g.query(prefixed)
            rows   = [dict(zip(result.vars, row)) for row in result]
            n_rows = len(rows)

            if rows:
                df = pd.DataFrame([
                    {str(k): str(v) if v else "" for k, v in r.items()}
                    for r in rows
                ])
                df.to_csv(out_dir / f"dl_query_{qid}.tsv", sep="\t", index=False)

            logger.info(f"  {qid}: {label}  → {n_rows} results")
            summary_rows.append({
                "query_id":    qid,
                "label":       label,
                "description": q.get("description", ""),
                "n_results":   n_rows,
                "status":      "OK",
            })

        except Exception as e:
            logger.warning(f"  {qid}: {label}  → ERROR: {e}")
            summary_rows.append({
                "query_id":    qid,
                "label":       label,
                "description": q.get("description", ""),
                "n_results":   0,
                "status":      f"ERROR: {e}",
            })

    return pd.DataFrame(summary_rows)


# ══════════════════════════════════════════════════════════════════════════════
# Stage runner
# ══════════════════════════════════════════════════════════════════════════════

def run(cfg: dict) -> None:
    logger  = setup_logging(cfg.get("log_dir", "logs"), "stage4b_owl_reasoning")
    stage_banner(logger, "STAGE 4b", "OWL Reasoning (Python — no Java required)")
    out_dir = make_output_dir(cfg, "stage4b_owl_reasoning")

    # ── Locate OWL files ──────────────────────────────────────────────────
    core_path = get_input_path(cfg, "core_ontology")
    pop_path  = get_input_path(cfg, "instance_data")

    # Also accept population.ttl from stage4 output
    if not pop_path or not Path(pop_path).exists():
        stage4_out = Path(cfg.get("output_dir","output")) / "stage4_owl_population" / "population.ttl"
        if stage4_out.exists():
            pop_path  = str(stage4_out)
            logger.info(f"Using stage4 output: {pop_path}")

    if not core_path or not Path(core_path).exists():
        logger.error("core_ontology not found — aborting stage 4b")
        logger.error("Set  core_ontology: input/core_ontology.ttl  in config")
        return

    if not pop_path or not Path(pop_path).exists():
        logger.error("Population TTL not found — run stage 4 first, "
                     "or set  instance_data: input/population.ttl  in config")
        return

    # ── Core namespace ────────────────────────────────────────────────────
    core_ns = cfg.get(
        "owl_core_ns",
        "http://www.semanticweb.org/32045/ontologies/2026/0/core_ontology#"
    )

    # ── Load user DL queries from YAML file (optional) ────────────────────
    extra_queries: list[dict] = []
    dl_query_file = get_input_path(cfg, "owl_dl_queries")
    if dl_query_file and Path(dl_query_file).exists():
        try:
            import yaml
            with open(dl_query_file) as fh:
                extra_queries = yaml.safe_load(fh) or []
            logger.info(f"Loaded {len(extra_queries)} user-defined DL queries")
        except Exception as e:
            logger.warning(f"Could not load DL query file: {e}")

    # ── Run reasoning ──────────────────────────────────────────────────────
    g, g_inferred, stats = run_reasoning(core_path, pop_path, logger)
    if g is None:
        return

    # ── Save inferred triples ──────────────────────────────────────────────
    if g_inferred and len(g_inferred):
        ttl_path = out_dir / "inferred_triples.ttl"
        g_inferred.serialize(destination=str(ttl_path), format="turtle")
        logger.info(f"Inferred triples saved → {ttl_path}")

    # ── Consistency report ────────────────────────────────────────────────
    consistency_rows = [
        {"metric": "consistent",         "value": stats.get("consistent", "N/A")},
        {"metric": "triples_before",     "value": stats.get("n_triples_before", 0)},
        {"metric": "triples_after",      "value": stats.get("n_triples_after", 0)},
        {"metric": "triples_inferred",   "value": stats.get("n_inferred", 0)},
        {"metric": "unsatisfiable_inds", "value": stats.get("n_unsatisfiable", 0)},
    ]
    save_tsv(pd.DataFrame(consistency_rows), out_dir / "consistency_report.tsv", logger)

    # ── Run DL queries ────────────────────────────────────────────────────
    logger.info(f"\nRunning {len(_BUILTIN_DL_QUERIES) + len(extra_queries)} "
                f"DL queries against reasoned graph …")
    summary_df = run_dl_queries(g, core_ns, out_dir, extra_queries, logger)
    save_tsv(summary_df, out_dir / "dl_query_summary.tsv", logger)

    # ── Final log ─────────────────────────────────────────────────────────
    logger.info("")
    logger.info("=" * 60)
    logger.info("  STAGE 4b COMPLETE")
    logger.info("=" * 60)
    logger.info(f"  Consistency:     {'PASS' if stats.get('consistent') else 'FAIL'}")
    logger.info(f"  Triples before:  {stats.get('n_triples_before',0):,}")
    logger.info(f"  Triples inferred:{stats.get('n_inferred',0):,}")
    ok  = (summary_df["status"] == "OK").sum()  if not summary_df.empty else 0
    logger.info(f"  DL queries:      {ok}/{len(summary_df)} passed")
    logger.info("")
    logger.info("  Output files:")
    logger.info("    inferred_triples.ttl     — new triples from reasoning")
    logger.info("    consistency_report.tsv   — pass/fail + counts")
    logger.info("    dl_query_summary.tsv     — all DL query results")
    logger.info("    dl_query_DL01.tsv …      — per-query result tables")
    logger.info("")
    logger.info("  Comparison with Protege:")
    logger.info("    owlrl = OWL RL reasoning  (rule-based, no Java)")
    logger.info("    HermiT in Protege = OWL DL reasoning (needs Java)")
    logger.info("    Both detect inconsistencies for annotation ontologies.")
    logger.info("    HermiT additionally supports nominals + cardinality.")
    logger.info("    For your cAMP ontology, owlrl coverage is sufficient.")


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import sys, os, argparse

    ap = argparse.ArgumentParser(
        description="Stage 4b — OWL Reasoning (pure Python, no Java)")
    ap.add_argument("--config", default="config/pipeline_config.yaml")
    ap.add_argument("--core",   default="", help="Path to core_ontology.ttl")
    ap.add_argument("--pop",    default="", help="Path to population.ttl")
    ap.add_argument("--core-ns",default="",
                    help="Core namespace ending with #")
    args = ap.parse_args()

    _root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(_root)

    cfg = load_config(args.config)
    if args.core:
        cfg["core_ontology"] = args.core
    if args.pop:
        cfg["instance_data"] = args.pop
    if args.core_ns:
        cfg["owl_core_ns"] = args.core_ns

    run(cfg)
