"""
scripts/stage2_rule_conformance.py
────────────────────────────────────────────────────────────────────────────
STAGE 2 — RULE-BASED CONFORMANCE CHECKING

Implements an extensible rule registry.  Each rule is a Python function
with signature:  rule_fn(occ, dt, ann, cfg) -> list[dict violation]

Rules are registered via the @register_rule decorator and can be added
without touching the main runner.

Built-in rules
──────────────
  R001  hasPositionalCategory points to DOP terms
  R002  hasFunctionTag points to DOF terms
  R003  hasBindingTargetCategory points to DOT terms
  R004  topology terms not in function_tags
  R005  function tags not in topology_context
  R006  binding-target terms not in function_tags
  R007  each occurrence maps to exactly one protein and one domain type
  R008  start < end
  R009  start and end are positive integers
  R010  blank primary keys
  R011  duplicate IDs
  R012  invalid namespace prefixes
  R013  copy_number_property in controlled vocabulary
  R014  topology_context in controlled vocabulary
  R015  proximity_categories in controlled vocabulary

Outputs
───────
  rule_violation_summary.tsv   — one row per rule, PASS/FAIL + count
  rule_violation_rows.tsv      — every violating row with rule + reason
────────────────────────────────────────────────────────────────────────────
"""

from __future__ import annotations

import re
from collections import defaultdict
from typing import Callable

import pandas as pd

import sys as _sys
import os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from scripts.utils import (
    build_cv_lookup,
    get_delimiter,
    get_input_path,
    is_cv_backed,
    load_config,
    load_table,
    make_output_dir,
    normalise_id_label,
    rename_to_canonical,
    safe_int,
    save_tsv,
    setup_logging,
    split_cell,
    stage_banner,
)

# ══════════════════════════════════════════════════════════════════════════
# Rule registry
# ══════════════════════════════════════════════════════════════════════════

_RULES: dict[str, dict] = {}


def register_rule(rule_id: str, name: str, description: str):
    """Decorator that registers a rule function."""
    def decorator(fn: Callable):
        _RULES[rule_id] = {
            "rule_id": rule_id,
            "name": name,
            "description": description,
            "fn": fn,
        }
        return fn
    return decorator


def get_all_rules() -> dict:
    return _RULES


# ══════════════════════════════════════════════════════════════════════════
# Helpers
# ══════════════════════════════════════════════════════════════════════════

def _viol(rule_id: str, row_idx, reason: str,
          source: str, row_data: dict) -> dict:
    """Create a violation record."""
    return {
        "_rule_id": rule_id,
        "_source": source,
        "_reason": reason,
        "_row_index": row_idx,
        **{k: str(v)[:200] for k, v in row_data.items()},
    }


def _row(df: pd.DataFrame, idx) -> dict:
    return df.loc[idx].to_dict()


# ── CV token helper ────────────────────────────────────────────────────────

# Matches any DOx namespace prefix pattern
_DOX_RE = re.compile(r"^(DO[A-Z]):\d+")

def _token_ok(tok: str, namespace: str,
               valid_set: frozenset[str]) -> bool:
    """
    Return True if a token is valid for a given namespace.

    Accepts ALL formats:
      "NearTransmembraneHelix"            plain name in CV
      "DOP:000001"                        bare ID with correct namespace
      "DOP:000001 NearTransmembraneHelix" combined — checks namespace prefix
                                          AND name independently
      Any of the above is sufficient.

    Args:
        tok:       the raw token from the data cell
        namespace: expected namespace prefix, e.g. "DOP", "DOF", "DOT"
        valid_set: frozenset of valid names/IDs from build_cv_lookup()
    """
    tok = tok.strip()
    if not tok:
        return True   # skip empty tokens
    id_part, lbl_part = normalise_id_label(tok)

    # 1. Check if the name part (or full token) is in the CV
    if is_cv_backed(tok, valid_set):
        return True

    # 2. Accept bare correct-namespace ID even if name not in CV yet
    #    e.g. "DOP:000099" — ID exists but name not registered yet
    ns_re = re.compile(rf"^{re.escape(namespace)}:\d+$")
    if ns_re.match(id_part):
        return True

    # 3. Accept combined "NSP:NNNNNN TermName" if namespace matches
    #    even if TermName not yet in CV — new terms added to ontology
    if ns_re.match(id_part) and lbl_part:
        return True

    return False


def _wrong_namespace(tok: str, expected_ns: str) -> bool:
    """
    Return True if the token has an explicit DOx: prefix that is NOT
    the expected namespace.
    e.g. _wrong_namespace("DOF:000009 Binding", "DOP") → True
         _wrong_namespace("NearTransmembraneHelix", "DOP") → False (no prefix)
         _wrong_namespace("DOP:000001", "DOP") → False
    """
    id_part, _ = normalise_id_label(tok.strip())
    m = _DOX_RE.match(id_part)
    if m:
        return m.group(1) != expected_ns
    return False


# ══════════════════════════════════════════════════════════════════════════
# Built-in rules
# ══════════════════════════════════════════════════════════════════════════

@register_rule("R001", "positional_category_DOP",
               "positional_category values must use DOP namespace or be a valid positional_or_topology_term")
def rule_r001(occ: pd.DataFrame, dt: pd.DataFrame,
              ann: pd.DataFrame, cfg: dict) -> list:
    """
    Accepts all of:
      "Cterminal"                  plain CV name
      "DOP:000003"                 bare DOP ID
      "DOP:000003 Cterminal"       combined ID + name
    Rejects:
      "DOF:000009 Binding"         wrong namespace (DOF instead of DOP)
      "randomtext"                 not in CV and no valid ID prefix
    """
    viols = []
    col = "positional_category"
    if col not in occ.columns:
        return viols
    cv = build_cv_lookup(ann)
    valid = cv.get("positional_or_topology_term", frozenset())
    for idx, val in occ[col].items():
        for tok in split_cell(val, get_delimiter(cfg, col)):
            if _wrong_namespace(tok, "DOP"):
                viols.append(_viol("R001", idx,
                                   f"wrong namespace in positional_category "
                                   f"(expected DOP:): '{tok}'",
                                   "occurrence", _row(occ, idx)))
            elif not _token_ok(tok, "DOP", valid):
                viols.append(_viol("R001", idx,
                                   f"unrecognised positional_category: '{tok}' "
                                   f"(expected DOP:NNNNNN, DOP:NNNNNN Name, or plain CV name)",
                                   "occurrence", _row(occ, idx)))
    return viols


@register_rule("R002", "function_tag_DOF",
               "function_tags values must use DOF namespace or be a valid functional_term")
def rule_r002(occ: pd.DataFrame, dt: pd.DataFrame,
              ann: pd.DataFrame, cfg: dict) -> list:
    """
    Accepts all of:
      "kinase"                     plain CV name
      "DOF:000009"                 bare DOF ID
      "DOF:000009 Binding"         combined ID + name
    Rejects:
      "DOP:000003 Cterminal"       wrong namespace (DOP instead of DOF)
      "randomtext"                 not in CV and no valid ID prefix
    """
    viols = []
    cv = build_cv_lookup(ann)
    valid = cv.get("functional_term", frozenset())
    for tbl, src in [(occ, "occurrence"), (dt, "domain_type")]:
        col = "function_tags"
        if col not in tbl.columns:
            continue
        for idx, val in tbl[col].items():
            for tok in split_cell(val, get_delimiter(cfg, col)):
                if _wrong_namespace(tok, "DOF"):
                    viols.append(_viol("R002", idx,
                                       f"wrong namespace in function_tags "
                                       f"(expected DOF:): '{tok}'",
                                       src, _row(tbl, idx)))
                elif not _token_ok(tok, "DOF", valid):
                    viols.append(_viol("R002", idx,
                                       f"unrecognised function_tag: '{tok}' "
                                       f"(expected DOF:NNNNNN, DOF:NNNNNN Name, or plain CV name)",
                                       src, _row(tbl, idx)))
    return viols


@register_rule("R003", "binding_target_DOT",
               "binding_parnter values must use DOT namespace or be a valid binding_target_term")
def rule_r003(occ: pd.DataFrame, dt: pd.DataFrame,
              ann: pd.DataFrame, cfg: dict) -> list:
    """
    Accepts all of:
      "Gs alpha subunit"           plain CV name
      "DOT:000079"                 bare DOT ID
      "DOT:000079 Mg2+"            combined ID + name
    Rejects:
      "DOP:000003 Cterminal"       wrong namespace (DOP instead of DOT)
      "randomtext"                 not in CV and no valid ID prefix
    """
    viols = []
    cv = build_cv_lookup(ann)
    valid = cv.get("binding_target_term", frozenset())
    # Handle the original typo 'binding_parnter' and any renamed form
    col = next(
        (c for c in ["binding_partner", "binding_parnter"] if c in occ.columns),
        None,
    )
    if not col:
        return viols
    for idx, val in occ[col].items():
        for tok in split_cell(val, get_delimiter(cfg, col)):
            if _wrong_namespace(tok, "DOT"):
                viols.append(_viol("R003", idx,
                                   f"wrong namespace in binding_partner "
                                   f"(expected DOT:): '{tok}'",
                                   "occurrence", _row(occ, idx)))
            elif not _token_ok(tok, "DOT", valid):
                viols.append(_viol("R003", idx,
                                   f"unrecognised binding_partner: '{tok}' "
                                   f"(expected DOT:NNNNNN, DOT:NNNNNN Name, or plain CV name)",
                                   "occurrence", _row(occ, idx)))
    return viols


@register_rule("R004", "topology_not_in_function",
               "topology terms must not appear in function_tags")
def rule_r004(occ: pd.DataFrame, dt: pd.DataFrame,
              ann: pd.DataFrame, cfg: dict) -> list:
    viols = []
    if "topology_context" not in occ.columns or "function_tags" not in occ.columns:
        return viols
    cv = build_cv_lookup(ann)
    topo_terms = cv.get("positional_or_topology_term", frozenset())
    func_delim = get_delimiter(cfg, "function_tags")
    for idx, row in occ[["topology_context", "function_tags"]].iterrows():
        func_toks = set(split_cell(row["function_tags"], func_delim))
        for tok in func_toks:
            if is_cv_backed(tok, topo_terms):
                viols.append(_viol("R004", idx,
                                   f"topology term '{tok}' in function_tags",
                                   "occurrence", _row(occ, idx)))
    return viols


@register_rule("R005", "function_not_in_topology",
               "function tags must not appear in topology_context")
def rule_r005(occ: pd.DataFrame, dt: pd.DataFrame,
              ann: pd.DataFrame, cfg: dict) -> list:
    viols = []
    if "function_tags" not in occ.columns or "topology_context" not in occ.columns:
        return viols
    cv = build_cv_lookup(ann)
    func_terms = cv.get("functional_term", frozenset())
    topo_delim = get_delimiter(cfg, "topology_context")
    for idx, row in occ[["function_tags", "topology_context"]].iterrows():
        topo_toks = set(split_cell(row["topology_context"], topo_delim))
        for tok in topo_toks:
            if is_cv_backed(tok, func_terms):
                viols.append(_viol("R005", idx,
                                   f"function term '{tok}' in topology_context",
                                   "occurrence", _row(occ, idx)))
    return viols


@register_rule("R006", "binding_not_in_function",
               "binding-target terms must not appear in function_tags")
def rule_r006(occ: pd.DataFrame, dt: pd.DataFrame,
              ann: pd.DataFrame, cfg: dict) -> list:
    viols = []
    # After rename_to_canonical the typo 'binding_parnter' → 'binding_partner'
    bind_col = next(
        (c for c in ["binding_partner", "binding_parnter"] if c in occ.columns),
        None,
    )
    if not bind_col or "function_tags" not in occ.columns:
        return viols
    cv = build_cv_lookup(ann)
    bind_terms = cv.get("binding_target_term", frozenset())
    func_delim = get_delimiter(cfg, "function_tags")
    for idx, row in occ[[bind_col, "function_tags"]].iterrows():
        func_toks = set(split_cell(row["function_tags"], func_delim))
        for tok in func_toks:
            # Skip tokens that explicitly carry a non-binding namespace prefix.
            # DOF: = function ontology  →  never a binding target
            # DOP: = positional/topology →  never a binding target
            # Only DOT: prefixed tokens (or plain names) can be binding targets.
            id_part, _ = normalise_id_label(tok)
            ns = id_part.split(":")[0].upper() if ":" in id_part else ""
            if ns in ("DOF", "DOP"):
                continue
            if is_cv_backed(tok, bind_terms):
                viols.append(_viol("R006", idx,
                                   f"binding term '{tok}' in function_tags",
                                   "occurrence", _row(occ, idx)))
    return viols


@register_rule("R007", "one_protein_one_domain",
               "each occurrence must have exactly one protein_accession and one domain_type_id")
def rule_r007(occ: pd.DataFrame, dt: pd.DataFrame,
              ann: pd.DataFrame, cfg: dict) -> list:
    viols = []
    for col in ["protein_accession", "domain_type_id"]:
        if col not in occ.columns:
            continue
        for idx, val in occ[col].items():
            toks = split_cell(val, ";")
            if len(toks) != 1:
                viols.append(_viol("R007", idx,
                                   f"{col} has {len(toks)} values (expected 1)",
                                   "occurrence", _row(occ, idx)))
    return viols


@register_rule("R008", "start_less_than_end",
               "start position must be strictly less than end position")
def rule_r008(occ: pd.DataFrame, dt: pd.DataFrame,
              ann: pd.DataFrame, cfg: dict) -> list:
    viols = []
    if "start" not in occ.columns or "end" not in occ.columns:
        return viols
    for idx, row in occ[["start", "end"]].iterrows():
        si = safe_int(row["start"])
        ei = safe_int(row["end"])
        if si is not None and ei is not None and si >= ei:
            viols.append(_viol("R008", idx,
                               f"start({si}) >= end({ei})",
                               "occurrence", _row(occ, idx)))
    return viols


@register_rule("R009", "start_end_positive_integers",
               "start and end must be positive integers (> 0)")
def rule_r009(occ: pd.DataFrame, dt: pd.DataFrame,
              ann: pd.DataFrame, cfg: dict) -> list:
    viols = []
    for col in ["start", "end"]:
        if col not in occ.columns:
            continue
        for idx, val in occ[col].items():
            if pd.isna(val) or str(val).strip() == "":
                continue
            n = safe_int(val)
            if n is None:
                viols.append(_viol("R009", idx,
                                   f"{col} not int: {val}",
                                   "occurrence", _row(occ, idx)))
            elif n <= 0:
                viols.append(_viol("R009", idx,
                                   f"{col} not positive: {n}",
                                   "occurrence", _row(occ, idx)))
    return viols


@register_rule("R010", "blank_primary_keys",
               "primary key columns must not be blank")
def rule_r010(occ: pd.DataFrame, dt: pd.DataFrame,
              ann: pd.DataFrame, cfg: dict) -> list:
    viols = []
    checks = [
        (occ, "occurrence_id", "occurrence"),
        (dt,  "domain_type_id", "domain_type"),
        (ann, "ID", "annotation"),
    ]
    for tbl, col, src in checks:
        if col not in tbl.columns:
            continue
        blank = tbl[col].isna() | (tbl[col].astype(str).str.strip() == "")
        for idx in tbl[blank].index:
            viols.append(_viol("R010", idx,
                               f"blank {col}", src, _row(tbl, idx)))
    return viols


@register_rule("R011", "duplicate_ids",
               "primary key IDs must be unique within each table")
def rule_r011(occ: pd.DataFrame, dt: pd.DataFrame,
              ann: pd.DataFrame, cfg: dict) -> list:
    viols = []
    checks = [
        (occ, "occurrence_id", "occurrence"),
        (dt,  "domain_type_id", "domain_type"),
        (ann, "ID", "annotation"),
    ]
    for tbl, col, src in checks:
        if col not in tbl.columns:
            continue
        dup = tbl.duplicated(subset=col, keep=False)
        for idx in tbl[dup].index:
            viols.append(_viol("R011", idx,
                               f"dup {col}: {tbl.loc[idx, col]}",
                               src, _row(tbl, idx)))
    return viols


# Built-in domain database namespace prefixes that are ALWAYS valid,
# regardless of what namespace_patterns is set to in pipeline_config.yaml.
# All domain_type_id values use the full prefix form:
#   INTERPRO:IPR003113, PFAM:PF00001, CDD:cd00001,
#   SMART:SM00003, NCBIFAM:TIGR00231, PRINTS:PR00008, PROFILE:PS50001
_BUILTIN_VALID_PREFIXES: frozenset = frozenset({
    # Domain database prefixes (full form as stored in domain_type_id)
    "INTERPRO", "PFAM", "CDD", "SMART", "NCBIFAM", "PRINTS", "PROFILE",
    # Domain ontology namespaces
    "DOT", "DOF", "DOP",
    # External cross-reference namespaces
    "GO", "CHEBI", "IPR",
})


@register_rule("R012", "invalid_namespace_prefixes",
               "all IDs must use recognised namespace prefixes")
def rule_r012(occ: pd.DataFrame, dt: pd.DataFrame,
              ann: pd.DataFrame, cfg: dict) -> list:
    viols = []

    # Merge built-in prefixes with any additional ones from config
    valid_prefixes = (
        _BUILTIN_VALID_PREFIXES
        | set(cfg.get("namespace_patterns", {}).keys())
    )

    # Columns to skip entirely (e.g. pipe-delimited occurrence_id)
    exclude_cols = set(
        cfg.get("validation", {}).get("namespace_check_exclude_columns", [])
    )

    id_cols = [
        (occ, "domain_type_id", "occurrence"),
        (dt,  "domain_type_id", "domain_type"),
        (ann, "ID", "annotation"),
    ]
    prefix_re = re.compile(r"^([A-Z][A-Z0-9_]+):")
    for tbl, col, table_src in id_cols:
        if col not in tbl.columns or col in exclude_cols:
            continue
        for idx, val in tbl[col].items():
            if pd.isna(val) or str(val).strip() == "":
                continue
            m = prefix_re.match(str(val).strip())
            if m and m.group(1) not in valid_prefixes:
                viols.append(_viol("R012", idx,
                                   f"unknown prefix in {col}: {val}",
                                   table_src, _row(tbl, idx)))
    return viols


@register_rule("R013", "copy_number_cv",
               "copy_number_property must be in CopyNumberCategory (plain name, DOP:NNNNNN, or DOP:NNNNNN Name)")
def rule_r013(occ: pd.DataFrame, dt: pd.DataFrame,
              ann: pd.DataFrame, cfg: dict) -> list:
    viols = []
    col = "copy_number_property"
    if col not in occ.columns:
        return viols
    cv = build_cv_lookup(ann)
    valid = cv.get("CopyNumberCategory", frozenset())
    for idx, val in occ[col].items():
        for tok in split_cell(val, get_delimiter(cfg, col)):
            if _wrong_namespace(tok, "DOP"):
                viols.append(_viol("R013", idx,
                                   f"wrong namespace in copy_number_property "
                                   f"(expected DOP:): '{tok}'",
                                   "occurrence", _row(occ, idx)))
            elif not _token_ok(tok, "DOP", valid):
                viols.append(_viol("R013", idx,
                                   f"invalid copy_number_property: '{tok}'",
                                   "occurrence", _row(occ, idx)))
    return viols


@register_rule("R014", "topology_context_cv",
               "topology_context must be in positional_or_topology_term (plain name, DOP:NNNNNN, or DOP:NNNNNN Name)")
def rule_r014(occ: pd.DataFrame, dt: pd.DataFrame,
              ann: pd.DataFrame, cfg: dict) -> list:
    viols = []
    col = "topology_context"
    if col not in occ.columns:
        return viols
    cv = build_cv_lookup(ann)
    valid = cv.get("positional_or_topology_term", frozenset())
    for idx, val in occ[col].items():
        for tok in split_cell(val, get_delimiter(cfg, col)):
            if _wrong_namespace(tok, "DOP"):
                viols.append(_viol("R014", idx,
                                   f"wrong namespace in topology_context "
                                   f"(expected DOP:): '{tok}'",
                                   "occurrence", _row(occ, idx)))
            elif not _token_ok(tok, "DOP", valid):
                viols.append(_viol("R014", idx,
                                   f"invalid topology_context: '{tok}'",
                                   "occurrence", _row(occ, idx)))
    return viols


@register_rule("R015", "proximity_categories_cv",
               "proximity_categories must be in positional_or_topology_term (plain name, DOP:NNNNNN, or DOP:NNNNNN Name)")
def rule_r015(occ: pd.DataFrame, dt: pd.DataFrame,
              ann: pd.DataFrame, cfg: dict) -> list:
    viols = []
    col = "proximity_categories"
    if col not in occ.columns:
        return viols
    cv = build_cv_lookup(ann)
    valid = cv.get("positional_or_topology_term", frozenset())
    for idx, val in occ[col].items():
        for tok in split_cell(val, get_delimiter(cfg, col)):
            if _wrong_namespace(tok, "DOP"):
                viols.append(_viol("R015", idx,
                                   f"wrong namespace in proximity_categories "
                                   f"(expected DOP:): '{tok}'",
                                   "occurrence", _row(occ, idx)))
            elif not _token_ok(tok, "DOP", valid):
                viols.append(_viol("R015", idx,
                                   f"invalid proximity_categories: '{tok}'",
                                   "occurrence", _row(occ, idx)))
    return viols


# ══════════════════════════════════════════════════════════════════════════
# Runner
# ══════════════════════════════════════════════════════════════════════════

def run(cfg: dict) -> None:
    log_dir = cfg.get("log_dir", "logs")
    logger = setup_logging(log_dir, "stage2_rule_conformance")
    stage_banner(logger, "STAGE 2", "Rule-Based Conformance Checking")
    out_dir = make_output_dir(cfg, "stage2_conformance")

    # Load tables
    def _load(key):
        p = get_input_path(cfg, key)
        return load_table(p) if p else pd.DataFrame()

    ann = _load("annotation_table")
    dt  = _load("domain_type_summary")
    occ = _load("domain_occurrence_table")

    if ann.empty and dt.empty and occ.empty:
        logger.error("No input tables found — cannot run conformance checks")
        return

    # Rename to canonical
    aliases = cfg.get("column_aliases", {})
    for tbl in [dt, occ]:
        rename_to_canonical(tbl, aliases)

    # Run all registered rules
    all_violations: list[dict] = []
    summary_rows: list[dict] = []

    for rule_id, rule in sorted(_RULES.items()):
        logger.info(f"Running {rule_id}: {rule['name']} …")
        try:
            viols = rule["fn"](occ, dt, ann, cfg)
        except Exception as exc:
            logger.error(f"  Rule {rule_id} raised: {exc}")
            viols = []
        n = len(viols)
        status = "FAIL" if n > 0 else "PASS"
        summary_rows.append({
            "rule_id": rule_id,
            "name": rule["name"],
            "description": rule["description"],
            "status": status,
            "n_violations": n,
        })
        all_violations.extend(viols)
        if n:
            logger.warning(f"  {rule_id} FAIL — {n} violations")
        else:
            logger.info(f"  {rule_id} PASS")

    summary_df = pd.DataFrame(summary_rows)
    save_tsv(summary_df, out_dir / "rule_violation_summary.tsv", logger)

    if all_violations:
        viol_df = pd.DataFrame(all_violations)
        meta = [c for c in ["_rule_id", "_source", "_reason", "_row_index"]
                if c in viol_df.columns]
        other = [c for c in viol_df.columns if c not in meta]
        save_tsv(viol_df[meta + other], out_dir / "rule_violation_rows.tsv", logger)
    else:
        save_tsv(pd.DataFrame(), out_dir / "rule_violation_rows.tsv", logger)

    fails = sum(1 for r in summary_rows if r["status"] == "FAIL")
    logger.info(f"Stage 2 complete — {fails}/{len(summary_rows)} rules failed.")


if __name__ == "__main__":
    import sys, os
    # Always run from the ontology_pipeline project root so that
    # relative paths (input/, output/, config/) resolve correctly
    # regardless of where PyCharm / the terminal launched the script from.
    _project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(_project_root)
    _cfg_path = sys.argv[1] if len(sys.argv) > 1 else "config/pipeline_config.yaml"
    run(load_config(_cfg_path))
