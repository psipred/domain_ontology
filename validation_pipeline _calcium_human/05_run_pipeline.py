#!/usr/bin/env python3
"""
run_pipeline.py
════════════════════════════════════════════════════════════════════════════════
cAMP Pathway Domain Ontology Testing Pipeline — Master Runner

Usage
─────
  python run_pipeline.py                                    # all enabled stages
  python run_pipeline.py --stage 4                         # single stage
  python run_pipeline.py --stages 1,2,3,4,4b,5,6          # selected stages
  python run_pipeline.py --skip 8,9                        # skip descriptive stages
  python run_pipeline.py --list                            # show all stages

Stage overview
──────────────
  1   Table Validation             — required, data quality
  2   Rule Conformance             — required, cross-column rules
  3   SHACL Validation             — required, OWL shape constraints
  4   OWL Population Builder       — required, builds population.ttl (no Java)
  4b  OWL Reasoning                — required, consistency + DL queries (no Java)
  5   Summary Metrics              — fast, entity counts for Table 1
  6   Competency Questions         — required, SPARQL + pandas queries
  7   Domain Enrichment            — required, Fisher exact + BH-FDR
  8   Heatmaps                     — optional (supplementary figures)
  9   Network Analysis             — optional (supplementary figures)
  10  Robustness Testing           — required, bootstrap + label shuffling
  11  Jena SPARQL Runner           — DISABLED (requires Java)
  12  Known-Biology Benchmark      — required, Precision/Recall/F1
  13  Protein Module Separation    — required, clustering + PCoA
  14  External Database Agreement  — required, GO/Reactome/KEGG Jaccard
  15  Negative Controls            — required, real vs random comparison
  16  Sequence-Ontology Congruence — required (needs MMseqs2 clusters)
════════════════════════════════════════════════════════════════════════════════
"""

from __future__ import annotations

import argparse
import importlib
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from scripts.utils import load_config, make_output_dir, setup_logging

# ══════════════════════════════════════════════════════════════════════════════
# Stage registry
# (module_name, display_name, is_java)
# Stage "4b" is stored as key 40 internally to allow integer sorting
# ══════════════════════════════════════════════════════════════════════════════

STAGES: dict[int, tuple[str, str, bool]] = {
    1:  ("scripts.stage1_table_validation",             "Table Validation",                  False),
    2:  ("scripts.stage2_rule_conformance",             "Rule-Based Conformance",             False),
    3:  ("scripts.stage3_shacl",                        "SHACL Validation",                   False),
    4:  ("scripts.stage4_owl_population",               "OWL Population Builder (Python)",    False),
    40: ("scripts.stage4b_owl_reasoning",               "OWL Reasoning — Python/owlrl",       False),
    5:  ("scripts.stage5_summary_metrics",              "Ontology Summary Metrics",           False),
    6:  ("scripts.stage6_competency_questions",         "Competency Questions (SPARQL)",      False),
    7:  ("scripts.stage7_enrichment",                   "Domain Enrichment Analysis",         False),
    8:  ("scripts.stage8_heatmaps",                     "Domain-Function Heatmaps",           False),
    9:  ("scripts.stage9_network",                      "Ontology-Derived Network",           False),
    10: ("scripts.stage10_robustness",                  "Robustness Testing",                 False),
    11: ("",                                            "Jena SPARQL Runner (Java — disabled)",True),
    12: ("scripts.stage12_benchmark",                   "Known-Biology Benchmark",            False),
    13: ("scripts.stage13_protein_modules",             "Protein Module Separation",          False),
    14: ("scripts.stage14_external_agreement",          "External Database Agreement",        False),
    15: ("scripts.stage15_negative_controls",           "Negative-Control Validation",        False),
    16: ("scripts.stage16_sequence_ontology_validation","Sequence-Ontology Congruence",       False),
}

# Human-readable stage names including "4b"
STAGE_LABELS: dict[int, str] = {
    **{k: str(k) for k in STAGES if k != 40},
    40: "4b",
}

JAVA_STAGES: dict[int, dict] = {
    11: {
        "jar": "java_jena/target/ontology-jena.jar",
        "main": "org.ontology.JenaSparqlRunner",
        "args_fn": lambda cfg: [
            str(Path(cfg.get("input_dir","input")) / cfg.get("instance_data","instance_data.ttl")),
            cfg.get("query_dir","queries"),
            str(Path(cfg.get("output_dir","output")) / "stage11_jena"),
        ],
    },
}

CONFIG_KEY_MAP: dict[int, str] = {
    1:  "stage1_table_validation",
    2:  "stage2_rule_conformance",
    3:  "stage3_shacl",
    4:  "stage4_owl_population",
    40: "stage4b_owl_reasoning",
    5:  "stage5_summary_metrics",
    6:  "stage6_competency_questions",
    7:  "stage7_enrichment",
    8:  "stage8_heatmaps",
    9:  "stage9_network",
    10: "stage10_robustness",
    11: "stage11_jena_queries",
    12: "stage12_benchmark",
    13: "stage13_protein_modules",
    14: "stage14_external_agreement",
    15: "stage15_negative_controls",
    16: "stage16_sequence_ontology_validation",
}


def _parse_stage_arg(s: str) -> int:
    """Convert stage string to internal int key (handles '4b' → 40)."""
    s = s.strip().lower()
    if s == "4b":
        return 40
    return int(s)


def _stage_label(stage_id: int) -> str:
    return STAGE_LABELS.get(stage_id, str(stage_id))


def _stage_enabled(stage_id: int, cfg: dict) -> bool:
    stages_cfg = cfg.get("stages", {})
    key = CONFIG_KEY_MAP.get(stage_id)
    return stages_cfg.get(key, True) if key else True


def _run_python_stage(stage_id: int, module_name: str,
                       cfg: dict, logger) -> tuple[str, float]:
    start = time.time()
    try:
        mod = importlib.import_module(module_name)
        mod.run(cfg)
        elapsed = time.time() - start
        return "OK", elapsed
    except Exception as exc:
        elapsed = time.time() - start
        logger.error(f"Stage {_stage_label(stage_id)} raised exception: {exc}",
                     exc_info=True)
        return f"ERROR: {exc}", elapsed


def _run_java_stage(stage_id: int, cfg: dict, logger) -> tuple[str, float]:
    jcfg = JAVA_STAGES[stage_id]
    jar_path = Path(jcfg["jar"])
    if not jar_path.exists():
        msg = f"JAR not found: {jar_path}"
        logger.warning(f"Stage {_stage_label(stage_id)} skipped — {msg}")
        return f"SKIPPED ({msg})", 0.0
    java_args = jcfg["args_fn"](cfg)
    if len(java_args) >= 2:
        Path(java_args[-1]).mkdir(parents=True, exist_ok=True)
    cmd = ["java", "-cp", str(jar_path), jcfg["main"]] + java_args
    logger.info(f"Running Java: {' '.join(cmd)}")
    start = time.time()
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        elapsed = time.time() - start
        if result.returncode != 0:
            logger.error(f"Java stage {_stage_label(stage_id)} stderr:\n{result.stderr[:2000]}")
            return f"ERROR (exit {result.returncode})", elapsed
        logger.info(result.stdout[:2000])
        return "OK", elapsed
    except subprocess.TimeoutExpired:
        return "ERROR (timeout)", time.time() - start
    except FileNotFoundError:
        return "ERROR (java not found — ensure JDK 11+ is on PATH)", 0.0


def run_pipeline(cfg: dict, selected_stages: list[int],
                 skip_stages: list[int]) -> None:
    log_dir = cfg.get("log_dir", "logs")
    logger  = setup_logging(log_dir, "run_pipeline")

    stage_labels = [_stage_label(s) for s in selected_stages]
    logger.info("=" * 66)
    logger.info("  cAMP Pathway Domain Ontology Testing Pipeline")
    logger.info(f"  Project: {cfg.get('project_name','')}")
    logger.info(f"  Species: {cfg.get('species','')}")
    logger.info(f"  Stages : {stage_labels}")
    logger.info("=" * 66)

    make_output_dir(cfg)
    summary_rows: list[dict] = []

    for stage_id in sorted(selected_stages):
        label = _stage_label(stage_id)

        if stage_id in skip_stages:
            logger.info(f"[Stage {label}] SKIPPED (--skip)")
            summary_rows.append({"stage": label, "name": STAGES[stage_id][1],
                                  "status": "SKIPPED", "elapsed_s": 0})
            continue

        if not _stage_enabled(stage_id, cfg):
            logger.info(f"[Stage {label}] DISABLED in config")
            summary_rows.append({"stage": label, "name": STAGES[stage_id][1],
                                  "status": "DISABLED", "elapsed_s": 0})
            continue

        module_name, display_name, is_java = STAGES[stage_id]
        logger.info("")
        logger.info(f"┌─ STAGE {label}: {display_name} " + "─" * max(0, 40-len(display_name)))

        if is_java:
            status, elapsed = _run_java_stage(stage_id, cfg, logger)
        else:
            status, elapsed = _run_python_stage(stage_id, module_name, cfg, logger)

        logger.info(f"└─ Stage {label} {status}  ({elapsed:.1f}s)")
        summary_rows.append({"stage": label, "name": display_name,
                              "status": status, "elapsed_s": round(elapsed, 1)})

    import pandas as pd
    summary_df   = pd.DataFrame(summary_rows)
    summary_path = Path(cfg.get("output_dir", "output")) / "pipeline_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)

    logger.info("")
    logger.info("=" * 66)
    logger.info("  PIPELINE COMPLETE")
    logger.info("=" * 66)
    ok   = sum(1 for r in summary_rows if r["status"] == "OK")
    errs = sum(1 for r in summary_rows if "ERROR" in str(r["status"]))
    skip = len(summary_rows) - ok - errs
    logger.info(f"  {len(summary_rows)} stages  |  {ok} OK  |  {errs} errors  |  {skip} skipped/disabled")
    logger.info(f"  Summary: {summary_path}")
    if errs:
        logger.warning("  Some stages failed — check logs/ for details")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="cAMP Pathway Ontology Testing Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--config", "-c", default="config/pipeline_config.yaml")
    parser.add_argument("--stage",  "-s", type=str,
                        help="Run a single stage (e.g. 1, 4, 4b, 13)")
    parser.add_argument("--stages",       type=str,
                        help="Comma-separated stages (e.g. 1,2,3,4,4b)")
    parser.add_argument("--skip",         type=str, default="",
                        help="Comma-separated stages to skip (e.g. 8,9,11)")
    parser.add_argument("--list", "-l",   action="store_true",
                        help="List all available stages and exit")
    args = parser.parse_args()

    if args.list:
        print("\nAvailable stages:")
        for sid in sorted(STAGES):
            label = _stage_label(sid)
            _, name, is_java = STAGES[sid]
            tag = "[Java — disabled]" if is_java else "[Python]"
            print(f"  {label:>3}  {tag:<20}  {name}")
        print()
        return

    cfg_path = Path(args.config)
    if not cfg_path.exists():
        print(f"[ERROR] Config not found: {cfg_path}")
        sys.exit(1)
    cfg = load_config(cfg_path)

    if args.stage:
        selected = [_parse_stage_arg(args.stage)]
    elif args.stages:
        selected = [_parse_stage_arg(s) for s in args.stages.split(",")]
    else:
        selected = sorted(STAGES.keys())

    skip = [_parse_stage_arg(s) for s in args.skip.split(",") if s.strip()]

    invalid = [s for s in selected if s not in STAGES]
    if invalid:
        print(f"[ERROR] Unknown stages: {[_stage_label(s) for s in invalid]}")
        print(f"Valid: {[_stage_label(s) for s in sorted(STAGES)]}")
        sys.exit(1)

    run_pipeline(cfg, selected, skip)


if __name__ == "__main__":
    main()
