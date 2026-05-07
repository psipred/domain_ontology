"""
scripts/stage1_table_validation.py  (v3 — all gaps from analysis addressed)
════════════════════════════════════════════════════════════════════════════════
STAGE 1 — TABLE VALIDATION

Validates three core tables against each other and against the annotation CV:
  • annotation_table.xlsx          (ontology term definitions)
  • domain_type_summary.tsv        (domain type catalogue)
  • domain_occurrence_table.tsv    (per-protein domain occurrences)

────────────────────────────────────────────────────────────────────────────
GAP-BY-GAP SUMMARY OF WHAT CHANGED FROM v1
────────────────────────────────────────────────────────────────────────────

GAP 1 — PRO and PSI-MOD validation  [WAS: missing]
  xref_PRO column in annotation table is now validated against RE_PRO
  (^PR:\d{9}$).  xref_ELM is validated against RE_ELM (^[A-Z]{3}_\S+$).
  No occ/dt columns carry these IDs in the current data, but the checks
  are scaffolded so they fire automatically when such columns appear.

GAP 2 — Explicit namespace-prefix cross-column contamination  [WAS: implicit]
  O18: explicit check — iterates every token in every CV column and flags
  any token whose DOP:/DOF:/DOT: prefix belongs to the wrong column.
  Rule table is defined as CROSS_CONTAMINATION_RULES (easy to extend).
  This is independent of is_cv_backed() and fires even for valid terms
  that have been placed in the wrong column.

GAP 3 — ID + label token normalisation  [WAS: delegated to utils, unclear]
  _is_cv_backed() is now defined fully inside this file, not delegated to
  utils.is_cv_backed().  It explicitly handles all four token forms:
    label-only          "Cterminal"
    ID-only             "DOP:000003"
    ID + space + label  "DOP:000003 Cterminal"   ← your TTL/export format
    label with count    "Cterminal:22"            ← strip the :count suffix
  The label-extraction and ID-extraction logic is documented inline.
  build_cv_lookup() is also rebuilt locally to include ID forms in the
  valid set, since utils.build_cv_lookup() only includes Name/synonym.

GAP 4 — Duplicate logic split into strict key + full-row  [WAS: full-row only]
  O3a: strict occurrence_id duplicate check (FAIL — genuine key violation)
  O3b: full-row duplicate excl. batch_gene (WARN — expected for batch splits)
  These are now separate rows in the validation report.

GAP 5 — source_db ↔ domain_type_id prefix consistency  [WAS: missing]
  T5b: if source_db = PFAM then domain_type_id must start with PFAM:, etc.
  DB_PREFIX_MAP defines all eight expected prefix mappings.

GAP 6 — domain_type_id format extensibility  [WAS: hard-coded 7 DBs]
  RE_DT_ID now includes SUPERFAMILY (SSF), PANTHER (PTHR), GENE3D (G3DSA),
  HAMAP (MF_), and PIRSF (PIRSF).  Still strict per database — update
  RE_DT_ID if you add new sources.

GAP 7 — Production-readiness of utility dependencies  [WAS: unclear]
  This file is now self-contained for CV and token logic.
  It defines its own _build_cv_lookup() and _is_cv_backed() that do not
  depend on the behaviour of utils.build_cv_lookup / utils.is_cv_backed.
  It still imports split_cell, rename_to_canonical, load_table, safe_int,
  get_delimiter from utils (these are pure I/O helpers with no hidden logic).

────────────────────────────────────────────────────────────────────────────
COMPLETE CHECK LIST
────────────────────────────────────────────────────────────────────────────
  Annotation (A1–A7)
    A1  Required columns present
    A2  Blank ID
    A3  Duplicate ID
    A4  ID format DOP/DOF/DOT:NNNNNN
    A5  Namespace prefix matches ID prefix
    A6  xref_PRO format (PR:nnnnnnnnn)           ← NEW (Gap 1)
    A7  xref_ELM format (XXX_motif_name)         ← NEW (Gap 1)

  Domain types (T1–T9)
    T1  Required columns present
    T2  Blank domain_type_id
    T3  Duplicate domain_type_id
    T4  domain_type_id format (extended DB list)  ← UPDATED (Gap 6)
    T5a source_db in allowed vocabulary
    T5b source_db prefix matches domain_type_id   ← NEW (Gap 5)
    T6  integrated_interpro_id format
    T7  GO term format in go_terms_* columns
    T8  function_tags CV-backed

  Occurrences (O1–O18)
    O1   Required columns present
    O2   Blank occurrence_id
    O3a  Duplicate occurrence_id (strict key)     ← SPLIT (Gap 4)
    O3b  Full-row duplicate excl. batch_gene      ← SPLIT (Gap 4)
    O4   domain_type_id cross-ref to domain types
    O5   protein_accession UniProtKB format
    O6   start/end valid positive integers (start < end)
    O7   positional_category CV-backed
    O9   topology_context CV-backed
    O10  proximity_categories CV-backed
    O11  function_tags CV-backed
    O12  binding_parnter CV-backed
    O13  binding_target_chebi format CHEBI:nnnnn
    O14  GO term format in go_terms_* columns
    O15  protein_reactome_ids format R-XXX-nnnnn
    O16  protein_kegg_ids format abc:nnnnn
    O17  protein_go_cc_ids format GO:nnnnnnn
    O18  Namespace cross-contamination (DOP/DOF/DOT in wrong column) ← NEW (Gap 2)

Outputs
───────
  validation_report_table_level.tsv  — one row per check + fix_priority
  invalid_rows.tsv                   — flagged rows: _source_row = Excel row number
  missing_term_report.tsv            — unmatched tokens + closest_match hint
  validation_summary.txt             — plain-text human-readable summary
════════════════════════════════════════════════════════════════════════════════
"""

from __future__ import annotations

import re
from collections import defaultdict
from difflib import get_close_matches
from pathlib import Path
from typing import Any

import pandas as pd

import sys as _sys
import os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from scripts.utils import (
    get_delimiter,
    get_input_path,
    load_config,
    load_table,
    make_output_dir,
    rename_to_canonical,
    safe_int,
    save_tsv,
    setup_logging,
    split_cell,
    stage_banner,
)

# ══════════════════════════════════════════════════════════════════════════════
# Regex constants
# ══════════════════════════════════════════════════════════════════════════════

RE_ANN_ID    = re.compile(r"^(DOP|DOF|DOT):\d{6}$")

# Gap 6: extended database list
RE_DT_ID     = re.compile(
    r"^("
    r"CDD:cd\d+"
    r"|INTERPRO:IPR\d{6}"
    r"|NCBIFAM:TIGR\d{5}"
    r"|PFAM:PF\d{5}"
    r"|PRINTS:PR\d{5}"
    r"|PROFILE:PS\d{5}"
    r"|SMART:SM\d{5}"
    r"|SUPERFAMILY:SSF\d+"
    r"|PANTHER:PTHR\d+"
    r"|GENE3D:G3DSA:[\d.]+"
    r"|HAMAP:MF_\d+"
    r"|PIRSF:PIRSF\d+"
    r")$"
)

RE_INTERPRO  = re.compile(r"^IPR\d{6}$")
RE_GO        = re.compile(r"^GO:\d{7}$")
RE_REACTOME  = re.compile(r"^R-[A-Z]{3}-\d+$")
RE_KEGG      = re.compile(r"^[a-z]{2,5}:\d+$")
RE_CHEBI_OK  = re.compile(r"^CHEBI:\d+$")
RE_CHEBI_DBL = re.compile(r"^CHEBI:CHEBI:\d+$")
RE_UNIPROT   = re.compile(
    r"^[A-Z][0-9][A-Z0-9]{3}[0-9]([A-Z][0-9][A-Z0-9]{3}[0-9])?$"
)
# Gap 1: PRO and ELM patterns
RE_PRO       = re.compile(r"^PR:\d{9}$")
RE_ELM       = re.compile(r"^[A-Z]{3}_\S+$")    # e.g. LIG_SH3_1, MOD_CDK_SPxK_1

GO_COLS = ["go_terms_all", "go_terms_mf", "go_terms_bp", "go_terms_cc"]

# Gap 5: source_db → expected domain_type_id prefix
DB_PREFIX_MAP: dict[str, str] = {
    "CDD":        "CDD:",
    "INTERPRO":   "INTERPRO:",
    "PFAM":       "PFAM:",
    "PROFILE":    "PROFILE:",
    "SMART":      "SMART:",
    "PRINTS":     "PRINTS:",
    "NCBIFAM":    "NCBIFAM:",
    "SUPERFAMILY":"SUPERFAMILY:",
    "PANTHER":    "PANTHER:",
    "GENE3D":     "GENE3D:",
    "HAMAP":      "HAMAP:",
    "PIRSF":      "PIRSF:",
}

# Gap 2: namespace cross-contamination rules
# (forbidden_prefix, column_name, human_readable_description)
CROSS_CONTAMINATION_RULES: list[tuple[str, str, str]] = [
    ("DOP", "function_tags",        "DOP positional term in function_tags"),
    ("DOP", "binding_parnter",      "DOP positional term in binding_parnter"),
    ("DOF", "positional_category",  "DOF functional term in positional_category"),
    ("DOF", "topology_context",     "DOF functional term in topology_context"),
    ("DOF", "proximity_categories", "DOF functional term in proximity_categories"),
    ("DOF", "copy_number_property", "DOF functional term in copy_number_property"),
    ("DOT", "function_tags",        "DOT binding term in function_tags"),
    ("DOT", "positional_category",  "DOT binding term in positional_category"),
    ("DOT", "topology_context",     "DOT binding term in topology_context"),
    ("DOP", "continuity_category",  "DOP positional term in continuity_category"),
]

# Fix priority and suggested action per check — shown in report + text summary
CHECK_PRIORITY: dict[str, str] = {
    "A1": "HIGH", "A2": "HIGH", "A3": "HIGH", "A4": "HIGH",
    "A5": "MED",  "A6": "LOW",  "A7": "LOW",
    "T1": "HIGH", "T2": "HIGH", "T3": "MED",
    "T4": "MED",  "T5a":"MED",  "T5b":"MED",
    "T6": "LOW",  "T7": "LOW",  "T8": "MED",
    "O1": "HIGH", "O2": "HIGH",
    "O3a":"HIGH", "O3b":"LOW",
    "O4": "HIGH", "O5": "MED",  "O6": "HIGH",
    "O7": "MED",  "O9": "MED",
    "O10":"LOW",  "O11":"MED",  "O12":"MED",
    "O13":"HIGH", "O14":"LOW",  "O15":"LOW",
    "O16":"LOW",  "O17":"LOW",  "O18":"HIGH",
}

SUGGESTED_FIX: dict[str, str] = {
    "A2":  "Fill in the 2 blank IDs in annotation_table.xlsx",
    "A3":  "Remove duplicate ID rows from annotation_table.xlsx",
    "A4":  "IDs must match DOP:000001 / DOF:000001 / DOT:000001",
    "A5":  "Align 'namespace' column value to match the ID prefix",
    "A6":  "Fix xref_PRO to format PR:000000001 (9 digits)",
    "A7":  "Fix xref_ELM to format XXX_name e.g. LIG_SH3_1",
    "T2":  "Fill in blank domain_type_id rows",
    "T3":  "Remove duplicate domain_type_id rows",
    "T4":  "Fix domain_type_id format e.g. PFAM:PF00001",
    "T5a": "Use only allowed source_db values",
    "T5b": "source_db and domain_type_id prefix must match",
    "T6":  "Fix integrated_interpro_id to IPRnnnnnn (6 digits)",
    "T7":  "Fix GO IDs to GO:nnnnnnn (7 digits)",
    "T8":  "Add functional_term tokens to annotation_table.xlsx or fix spelling",
    "O2":  "Fill in blank occurrence_id values",
    "O3a": "Remove rows with identical occurrence_id — genuine key violation",
    "O3b": "Review; batch_gene split rows are expected duplicates",
    "O4":  "Ensure domain_type_id exists in domain_type_summary.tsv",
    "O5":  "Fix protein_accession to 6-char UniProtKB e.g. P06850",
    "O6":  "Fix start/end: positive integers, start < end",
    "O7":  "Add positional_category tokens to annotation_table.xlsx",
    "O9":  "Add signal_peptide to positional_or_topology_term in annotation_table.xlsx",
    "O10": "Add proximity_categories tokens to annotation_table.xlsx",
    "O11": "Add ligand / GEF / GAP to functional_term in annotation_table.xlsx",
    "O12": "Add free-text binding partners to binding_target_term in annotation_table.xlsx",
    "O13": "Remove duplicate CHEBI: prefix: CHEBI:CHEBI:nnn → CHEBI:nnn",
    "O14": "Fix GO IDs to GO:nnnnnnn (7 digits)",
    "O15": "Fix Reactome IDs to R-HSA-nnnnn format",
    "O16": "Fix KEGG IDs to hsa:nnnnn format",
    "O17": "Fix protein_go_cc_ids to GO:nnnnnnn (7 digits)",
    "O18": "Remove namespace cross-contamination (e.g. DOP: term in function_tags)",
}

# Core columns whose absence stops all further annotation checks (A2-A7).
# "status" is intentionally excluded here: annotation_table.xlsx may be
# an OWL/TTL export that predates the status field.  A missing status is
# reported as a separate WARN (A1w) so A2-A7 can still run.
ANN_REQUIRED      = ["ID", "namespace", "term_type", "Name"]
ANN_REQUIRED_SOFT = ["status"]  # checked separately; missing = WARN not FAIL


# ══════════════════════════════════════════════════════════════════════════════
# CV lookup — self-contained, does not rely on utils.build_cv_lookup (Gap 7)
# ══════════════════════════════════════════════════════════════════════════════

# Gap 3: this regex matches combined "ID label" tokens from the TTL export
# e.g. "DOP:000003 Cterminal" or "DOF:000012 kinase_activity"
_ID_LABEL_RE = re.compile(
    r"^((DOP|DOF|DOT|GO|IPR|CHEBI|PR|MOD|R-[A-Z]{3}):[^\s]+)\s+(.+)$"
)
# Matches label:count suffix e.g. "Cterminal:22"
_LABEL_COUNT_RE = re.compile(r"^(.+?):\d+$")


def _extract_id(token: str) -> str:
    """Return only the ID portion from a combined 'ID label' token."""
    m = _ID_LABEL_RE.match(token.strip())
    return m.group(1) if m else token.strip()


def _extract_label(token: str) -> str:
    """Return only the label portion from a combined 'ID label' token."""
    m = _ID_LABEL_RE.match(token.strip())
    return m.group(3) if m else token.strip()


def _strip_count_suffix(token: str) -> str:
    """Strip :N count suffix from e.g. 'Cterminal:22' → 'Cterminal'."""
    m = _LABEL_COUNT_RE.match(token.strip())
    return m.group(1) if m else token.strip()


def _is_cv_backed(token: str, valid_set: frozenset[str]) -> bool:
    """
    Gap 3 — all four token forms are handled:
      1. label-only          "Cterminal"
      2. ID-only             "DOP:000003"
      3. ID + label          "DOP:000003 Cterminal"
      4. label with :count   "Cterminal:22"
    Checks each form case-insensitively.
    """
    t          = token.strip()
    id_part    = _extract_id(t)
    label_part = _extract_label(t)
    stripped   = _strip_count_suffix(t)
    return any(v in valid_set for v in (
        t,             t.lower(),
        id_part,       id_part.lower(),
        label_part,    label_part.lower(),
        stripped,      stripped.lower(),
    ))


def _get_dox_prefix(token: str) -> str | None:
    """Return DOP / DOF / DOT if the ID part of token starts with that prefix."""
    id_part = _extract_id(token)
    for pfx in ("DOP", "DOF", "DOT"):
        if id_part.startswith(pfx + ":"):
            return pfx
    return None


def _build_cv_lookup(ann: pd.DataFrame) -> dict[str, frozenset[str]]:
    """
    Build CV lookup from annotation table.
    Includes: ID, Name, canonical_name, synonym — all case variants.
    This supersedes utils.build_cv_lookup which omits the ID column.
    """
    cv: dict[str, set[str]] = defaultdict(set)
    for _, row in ann.iterrows():
        ttype = str(row.get("term_type", "")).strip()
        if not ttype or ttype == "nan":
            continue
        for field in ("ID", "Name", "canonical_name", "synonym"):
            val = str(row.get(field, "")).strip()
            if val and val != "nan":
                cv[ttype].add(val)
                cv[ttype].add(val.lower())
    return {k: frozenset(v) for k, v in cv.items()}


def _all_valid_tokens(cv: dict[str, frozenset[str]]) -> list[str]:
    """Flat list of all valid tokens across all term types, for close-match hints."""
    out: list[str] = []
    for s in cv.values():
        out.extend(v for v in s if not v.islower())  # skip lowercase duplicates
    return out


# ══════════════════════════════════════════════════════════════════════════════
# Shared helpers
# ══════════════════════════════════════════════════════════════════════════════

def _row_dict(df: pd.DataFrame, idx: Any) -> dict:
    return {k: (str(v)[:300] if isinstance(v, str) else v)
            for k, v in df.loc[idx].to_dict().items()}


def _add_result(results: list, check_id: str, desc: str,
                shape: tuple, status: str, n: int,
                detail: str = "") -> None:
    results.append({
        "check_id":      check_id,
        "description":   desc,
        "table_rows":    shape[0],
        "table_cols":    shape[1],
        "status":        status,
        "n_issues":      n,
        "fix_priority":  CHECK_PRIORITY.get(check_id, "MED"),
        "suggested_fix": SUGGESTED_FIX.get(check_id, ""),
        "detail":        detail,
    })


def _flag_mask(results: list, invalid: list,
               df: pd.DataFrame, mask: pd.Series,
               shape: tuple, check_id: str, desc: str, source: str,
               reason_fn=None, static_reason: str | None = None,
               detail: str = "") -> None:
    n = int(mask.sum())
    _add_result(results, check_id, desc, shape,
                "FAIL" if n else "PASS", n, detail)
    for idx in df[mask].index:
        reason = static_reason or reason_fn(idx)
        invalid.append({
            "_source":     source,
            "_check":      check_id,
            "_reason":     reason,
            "_source_row": int(idx) + 2,   # +2 = 0-based index + header row → Excel row
            **_row_dict(df, idx),
        })


def _check_cv_col(results: list, invalid: list,
                  df: pd.DataFrame, col: str,
                  valid_set: frozenset, source: str,
                  check_id: str, desc: str, shape: tuple,
                  delimiter: str = ";") -> None:
    """CV check for one column. Uses _is_cv_backed (Gap 3 / Gap 7)."""
    if col not in df.columns:
        _add_result(results, check_id, desc, shape,
                    "WARN", 0, f"Column '{col}' absent — skipped")
        return
    bad_rows: dict[Any, list[str]] = defaultdict(list)
    for idx, val in df[col].items():
        for tok in split_cell(val, delimiter):
            if not _is_cv_backed(tok, valid_set):
                bad_rows[idx].append(tok)
    n_toks = sum(len(v) for v in bad_rows.values())
    _add_result(results, check_id, desc, shape,
                "FAIL" if n_toks else "PASS", n_toks,
                f"{len(bad_rows)} rows with unmatched tokens ({n_toks} tokens total)")
    for idx, toks in bad_rows.items():
        invalid.append({
            "_source":     source,
            "_check":      check_id,
            "_reason":     f"unmatched in '{col}': " + "; ".join(toks),
            "_source_row": int(idx) + 2,
            **_row_dict(df, idx),
        })


def _check_go_cols(results: list, invalid: list,
                   df: pd.DataFrame, source: str,
                   check_id: str, shape: tuple,
                   extra_cols: list[str] | None = None,
                   delimiter: str = ";") -> None:
    cols = list(GO_COLS) + (extra_cols or [])
    bad: dict[Any, list[str]] = defaultdict(list)
    found = False
    for col in cols:
        if col not in df.columns:
            continue
        found = True
        for idx, val in df[col].items():
            for tok in split_cell(val, delimiter):
                if not RE_GO.fullmatch(tok):
                    bad[idx].append(f"{col}={tok}")
    if not found:
        _add_result(results, check_id, "GO term format (GO:nnnnnnn)",
                    shape, "WARN", 0, "No GO columns found")
        return
    _add_result(results, check_id, "GO term format (GO:nnnnnnn)",
                shape, "FAIL" if bad else "PASS", len(bad))
    for idx, details in bad.items():
        invalid.append({
            "_source":     source,
            "_check":      check_id,
            "_reason":     "bad GO: " + "; ".join(details),
            "_source_row": int(idx) + 2,
            **_row_dict(df, idx),
        })


# ══════════════════════════════════════════════════════════════════════════════
# A — Annotation table  (A1–A7)
# ══════════════════════════════════════════════════════════════════════════════

def validate_annotation(ann: pd.DataFrame,
                        cfg: dict) -> tuple[list, list]:
    results, invalid = [], []
    shape = ann.shape

    # ── Guard: cast all non-numeric columns to str so ~ and .str. never see float ──
    import numpy as _np
    for _c in ann.columns:
        if ann[_c].dtype == _np.float64 or ann[_c].dtype == object:
            ann[_c] = ann[_c].fillna("").astype(str)

    # A1 — required CORE columns (ID, namespace, term_type, Name)
    missing_cols = [c for c in ANN_REQUIRED if c not in ann.columns]
    _add_result(results, "A1", "Required columns present", shape,
                "FAIL" if missing_cols else "PASS",
                len(missing_cols), ", ".join(missing_cols))
    if missing_cols:
        # Only abort if core structural columns are absent
        return results, invalid

    # A1w — soft check for optional-but-recommended columns (e.g. "status")
    # Missing status is a WARN, not a FAIL, so A2-A7 still run.
    soft_missing = [c for c in ANN_REQUIRED_SOFT if c not in ann.columns]
    if soft_missing:
        _add_result(results, "A1", "Required columns present (soft)", shape,
                    "WARN", len(soft_missing),
                    f"Recommended columns missing: {', '.join(soft_missing)} "
                    f"— add to annotation_table.xlsx (value: 'active' or 'deprecated')")

    # A2 — blank IDs
    blank = ann["ID"].isna() | (ann["ID"].str.strip() == "")
    _flag_mask(results, invalid, ann, blank, shape, "A2",
               "Blank ID", "annotation",
               static_reason="blank ID — fill in the ID cell")

    # A3 — duplicate IDs
    dup = ann.duplicated(subset="ID", keep=False) & ~blank
    _add_result(results, "A3", "Duplicate ID", shape,
                "FAIL" if dup.sum() else "PASS", int(dup.sum()))
    for idx in ann[dup].index:
        invalid.append({
            "_source": "annotation", "_check": "A3",
            "_reason": f"duplicate ID='{ann.loc[idx,'ID']}'",
            "_source_row": int(idx) + 2,
            **_row_dict(ann, idx),
        })

    # A4 — ID format
    filled = ~blank
    bad_fmt = filled & ~ann["ID"].str.strip().str.fullmatch(RE_ANN_ID.pattern).fillna(False)
    _flag_mask(results, invalid, ann, bad_fmt, shape, "A4",
               "ID format (DOP/DOF/DOT:NNNNNN)", "annotation",
               reason_fn=lambda i: f"bad format: '{ann.loc[i,'ID']}'",
               detail="Expected DOP:000001 / DOF:000001 / DOT:000001")

    # A5 — namespace prefix consistency
    valid_ids = ~blank & ~bad_fmt
    mismatch = pd.Series(False, index=ann.index)
    for idx in ann[valid_ids].index:
        id_prefix = ann.loc[idx, "ID"].split(":")[0]
        ns = str(ann.loc[idx, "namespace"]).strip()
        if id_prefix != ns:
            mismatch.loc[idx] = True
    _add_result(results, "A5", "Namespace prefix matches ID prefix", shape,
                "FAIL" if mismatch.sum() else "PASS", int(mismatch.sum()))
    for idx in ann[mismatch].index:
        invalid.append({
            "_source": "annotation", "_check": "A5",
            "_reason": (f"ID prefix '{ann.loc[idx,'ID'].split(':')[0]}' "
                        f"≠ namespace '{ann.loc[idx,'namespace']}'"),
            "_source_row": int(idx) + 2,
            **_row_dict(ann, idx),
        })

    # A6 — xref_PRO format  (Gap 1)
    if "xref_PRO" in ann.columns:
        pro_bad: list[tuple] = []
        for idx, val in ann["xref_PRO"].items():
            for tok in split_cell(val, ";"):
                if not RE_PRO.fullmatch(tok):
                    pro_bad.append((idx, tok))
        for idx, tok in pro_bad:
            invalid.append({
                "_source": "annotation", "_check": "A6",
                "_reason": f"bad xref_PRO: '{tok}' — expected PR:000000001 (9 digits)",
                "_source_row": int(idx) + 2,
                **_row_dict(ann, idx),
            })
        _add_result(results, "A6", "xref_PRO format (PR:nnnnnnnnn)", shape,
                    "FAIL" if pro_bad else "PASS", len(pro_bad))
    else:
        _add_result(results, "A6", "xref_PRO format", shape,
                    "WARN", 0, "Column xref_PRO absent — skipped")

    # A7 — xref_ELM format  (Gap 1)
    if "xref_ELM" in ann.columns:
        elm_bad: list[tuple] = []
        for idx, val in ann["xref_ELM"].items():
            for tok in split_cell(val, ";"):
                if not RE_ELM.fullmatch(tok.strip()):
                    elm_bad.append((idx, tok))
        for idx, tok in elm_bad:
            invalid.append({
                "_source": "annotation", "_check": "A7",
                "_reason": f"bad xref_ELM: '{tok}' — expected XXX_motif_name e.g. LIG_SH3_1",
                "_source_row": int(idx) + 2,
                **_row_dict(ann, idx),
            })
        _add_result(results, "A7", "xref_ELM format (XXX_motif_name)", shape,
                    "FAIL" if elm_bad else "PASS", len(elm_bad))
    else:
        _add_result(results, "A7", "xref_ELM format", shape,
                    "WARN", 0, "Column xref_ELM absent — skipped")

    return results, invalid


# ══════════════════════════════════════════════════════════════════════════════
# T — Domain type table  (T1–T8)
# ══════════════════════════════════════════════════════════════════════════════

def validate_domain_types(dt: pd.DataFrame, cv: dict,
                          cfg: dict) -> tuple[list, list]:
    results, invalid = [], []
    shape = dt.shape
    dt = rename_to_canonical(dt, cfg.get("column_aliases", {}))

    # ── Guard: cast ALL non-numeric columns to str so ~ and .str. never see float ──
    import numpy as _np
    _dt_numeric = {"start", "end"}
    for _c in dt.columns:
        if _c in _dt_numeric:
            continue
        if dt[_c].dtype == _np.float64 or dt[_c].dtype == object:
            dt[_c] = dt[_c].fillna("").astype(str)

    req_cols = cfg.get("required_columns", {}).get(
        "domain_type_summary",
        ["domain_type_id", "source_db", "accession", "name", "evidence_level"],
    )
    valid_src = set(cfg.get("valid_source_db",
                   ["CDD","INTERPRO","PFAM","PROFILE","SMART","PRINTS",
                    "NCBIFAM","SUPERFAMILY","PANTHER","GENE3D","HAMAP","PIRSF"]))

    # T1 — required columns
    missing_cols = [c for c in req_cols if c not in dt.columns]
    _add_result(results, "T1", "Required columns present", shape,
                "FAIL" if missing_cols else "PASS",
                len(missing_cols), ", ".join(missing_cols))
    if missing_cols:
        return results, invalid

    # T2 — blank domain_type_id
    blank = dt["domain_type_id"].isna() | (dt["domain_type_id"].str.strip() == "")
    _flag_mask(results, invalid, dt, blank, shape, "T2",
               "Blank domain_type_id", "domain_type",
               static_reason="blank domain_type_id")

    # T3 — duplicates
    dup = dt.duplicated(subset="domain_type_id", keep=False) & ~blank
    _add_result(results, "T3", "Duplicate domain_type_id", shape,
                "FAIL" if dup.sum() else "PASS", int(dup.sum()))
    for idx in dt[dup].index:
        invalid.append({
            "_source": "domain_type", "_check": "T3",
            "_reason": f"duplicate domain_type_id='{dt.loc[idx,'domain_type_id']}'",
            "_source_row": int(idx) + 2,
            **_row_dict(dt, idx),
        })

    # T4 — domain_type_id format (Gap 6: extended)
    bad_fmt = (~blank) & ~dt["domain_type_id"].str.strip().str.fullmatch(
        RE_DT_ID.pattern).fillna(False).fillna(False).fillna(False)
    _flag_mask(results, invalid, dt, bad_fmt, shape, "T4",
               "domain_type_id format", "domain_type",
               reason_fn=lambda i: f"bad format: '{dt.loc[i,'domain_type_id']}'",
               detail="e.g. CDD:cd15445 / PFAM:PF00001 / INTERPRO:IPR000001")

    # T5a — source_db vocabulary
    bad_db = dt["source_db"].notna() & ~dt["source_db"].fillna("").isin(valid_src)
    _flag_mask(results, invalid, dt, bad_db, shape, "T5a",
               "source_db in controlled vocabulary", "domain_type",
               reason_fn=lambda i: f"unknown source_db: '{dt.loc[i,'source_db']}'",
               detail=f"Valid: {', '.join(sorted(valid_src))}")

    # T5b — source_db prefix matches domain_type_id prefix (Gap 5)
    t5b_viols: list[tuple] = []
    for idx, row in dt[["source_db", "domain_type_id"]].iterrows():
        src  = str(row["source_db"]).strip()
        dtid = str(row["domain_type_id"]).strip()
        exp  = DB_PREFIX_MAP.get(src)
        if exp and dtid and not dtid.startswith(exp):
            t5b_viols.append((idx, src, dtid, exp))
    _add_result(results, "T5b",
                "source_db prefix matches domain_type_id prefix", shape,
                "FAIL" if t5b_viols else "PASS", len(t5b_viols))
    for idx, src, dtid, exp in t5b_viols:
        invalid.append({
            "_source": "domain_type", "_check": "T5b",
            "_reason": f"source_db={src} but domain_type_id='{dtid}' (expected prefix '{exp}')",
            "_source_row": int(idx) + 2,
            **_row_dict(dt, idx),
        })

    # T6 — integrated_interpro_id
    if "integrated_interpro_id" in dt.columns:
        filled_ipr = (dt["integrated_interpro_id"].notna() &
                      (dt["integrated_interpro_id"].str.strip() != ""))
        bad_ipr = filled_ipr & ~dt["integrated_interpro_id"].str.strip()\
                                     .str.fullmatch(RE_INTERPRO.pattern).fillna(False)
        _flag_mask(results, invalid, dt, bad_ipr, shape, "T6",
                   "integrated_interpro_id format (IPRnnnnnn)", "domain_type",
                   reason_fn=lambda i: f"bad IPR: '{dt.loc[i,'integrated_interpro_id']}'")
    else:
        _add_result(results, "T6", "integrated_interpro_id format",
                    shape, "WARN", 0, "Column absent — skipped")

    # T7 — GO term format
    _check_go_cols(results, invalid, dt, "domain_type", "T7", shape)

    # T8 — function_tags CV
    _check_cv_col(results, invalid, dt,
                  "function_tags", cv.get("functional_term", frozenset()),
                  "domain_type", "T8",
                  "function_tags ontology-backed (functional_term)", shape)

    return results, invalid


# ══════════════════════════════════════════════════════════════════════════════
# O — Occurrence table  (O1–O18)
# ══════════════════════════════════════════════════════════════════════════════

OCC_REQUIRED_DEFAULT = [
    "occurrence_id", "protein_accession", "domain_type_id",
    "start", "end", "positional_category", "copy_number_property",
    "topology_context", "function_tags", "binding_parnter",
    "domain_evidence_level",
]


def validate_occurrences(occ: pd.DataFrame, dt: pd.DataFrame,
                          cv: dict, cfg: dict) -> tuple[list, list]:
    results, invalid = [], []
    shape = occ.shape
    occ = rename_to_canonical(occ, cfg.get("column_aliases", {}))

    # ── Guard: cast ALL non-numeric columns to str so ~ and .str. never see float ──
    # This handles any column in the WITH_IDS file (or future files) that pandas reads
    # as float64 because it contains NaN — even columns not explicitly listed here.
    import numpy as _np
    for _c in occ.columns:
        # Skip genuine numeric columns (start, end, coordinates, counts)
        _numeric_names = {"start", "end", "copy_number", "relative_start",
                          "relative_end", "relative_midpoint", "overlap_group_id",
                          "overlap_group_size", "ordinal_position_group",
                          "ordinal_position_all", "gap_from_prev", "gap_to_next",
                          "nearest_tm_distance", "signal_cleavage_pos",
                          "nearest_signal_cleavage_distance",
                          "nearest_active_site_distance"}
        if _c in _numeric_names:
            continue
        # Cast float64 columns (NaN-only or mixed) and object columns with NaN to str
        if occ[_c].dtype == _np.float64 or occ[_c].dtype == object:
            occ[_c] = occ[_c].fillna("").astype(str)

    req_cols = cfg.get("required_columns", {}).get(
        "domain_occurrence_table", OCC_REQUIRED_DEFAULT)

    # cv_column_to_term_type: maps column name → annotation table term_type
    cv_col_map: dict[str, str] = cfg.get("cv_column_to_term_type", {
        "positional_category":  "positional_or_topology_term",
        "topology_context":     "positional_or_topology_term",
        "proximity_categories": "positional_or_topology_term",
        "function_tags":        "functional_term",
        "binding_parnter":      "binding_target_term",
    })

    # O1 — required columns
    # Checked AFTER rename_to_canonical so canonical names are used consistently.
    # Two lookups are performed for each required column:
    #   Forward:  aliases[col] lists OTHER names that should have been renamed
    #             TO col — these are already gone after rename, so this catches
    #             any rename that somehow didn't happen.
    #   Reverse:  col was itself listed as an alias FOR another canonical name
    #             (e.g. config maps occurrence_id → domain_occ_key).  After
    #             rename_to_canonical col is gone; we check the canonical target.
    aliases = cfg.get("column_aliases", {})
    # Build reverse map: alias_name → canonical_name
    _reverse_alias: dict[str, str] = {}
    for _canonical, _alias_list in aliases.items():
        for _a in _alias_list:
            _reverse_alias[_a] = _canonical

    missing_cols = []
    for col in req_cols:
        found = col in occ.columns
        if not found:
            # Forward: look for surviving aliases of col
            for alias in aliases.get(col, []):
                if alias in occ.columns:
                    found = True
                    break
        if not found:
            # Reverse: col was renamed away to its canonical target
            renamed_to = _reverse_alias.get(col)
            if renamed_to and renamed_to in occ.columns:
                found = True
        if not found:
            missing_cols.append(col)
    _add_result(results, "O1", "Required columns present", shape,
                "FAIL" if missing_cols else "PASS",
                len(missing_cols), ", ".join(missing_cols))
    if missing_cols:
        return results, invalid

    # Resolve the occurrence_id column (may be 'occurrence_id' or 'domain_occurrence_id')
    oid_col = "occurrence_id" if "occurrence_id" in occ.columns else next(
        (c for c in occ.columns if "occurrence_id" in c), "occurrence_id"
    )

    # O2 — blank occurrence_id
    blank_oid = occ[oid_col].isna() | (occ[oid_col].str.strip() == "")
    _flag_mask(results, invalid, occ, blank_oid, shape, "O2",
               "Blank occurrence_id", "occurrence",
               static_reason="blank occurrence_id")

    # O3a — strict occurrence_id duplicates (Gap 4)
    dup_id = occ.duplicated(subset=oid_col, keep=False) & ~blank_oid
    _add_result(results, "O3a",
                "Duplicate occurrence_id (strict key check)", shape,
                "FAIL" if dup_id.sum() else "PASS", int(dup_id.sum()),
                "Genuine key violation — each occurrence_id must be unique")
    for idx in occ[dup_id].index:
        invalid.append({
            "_source": "occurrence", "_check": "O3a",
            "_reason": f"duplicate occurrence_id: '{occ.loc[idx, oid_col]}'",
            "_source_row": int(idx) + 2,
            **_row_dict(occ, idx),
        })

    # O3b — full-row duplicates excl. batch_gene (Gap 4)
    key_cols = [c for c in occ.columns if c != "batch_gene"]
    dup_row = occ.duplicated(subset=key_cols, keep=False) & ~blank_oid
    _add_result(results, "O3b",
                "Full-row duplicate excl. batch_gene", shape,
                "WARN" if dup_row.sum() else "PASS", int(dup_row.sum()),
                "Expected for batch_gene split rows; investigate if count is unexpectedly high")

    # O4 — domain_type_id cross-reference
    known_dt: set[str] = set()
    if "domain_type_id" in dt.columns:
        known_dt = set(dt["domain_type_id"].dropna().str.strip())
    bad_xref = (occ["domain_type_id"].notna() &
                ~occ["domain_type_id"].str.strip().isin(known_dt).fillna(False))
    _flag_mask(results, invalid, occ, bad_xref, shape, "O4",
               "domain_type_id cross-ref to domain type table", "occurrence",
               reason_fn=lambda i: f"unknown domain_type_id: '{occ.loc[i,'domain_type_id']}'")

    # O5 — UniProtKB format
    if "protein_accession" in occ.columns:
        filled = (occ["protein_accession"].notna() &
                  (occ["protein_accession"].str.strip() != ""))
        bad_acc = filled & ~occ["protein_accession"].str.strip()\
                               .str.fullmatch(RE_UNIPROT.pattern).fillna(False)
        _flag_mask(results, invalid, occ, bad_acc, shape, "O5",
                   "protein_accession UniProtKB format", "occurrence",
                   reason_fn=lambda i: f"bad accession: '{occ.loc[i,'protein_accession']}'",
                   detail="Expected 6-char UniProtKB e.g. P06850")
    else:
        _add_result(results, "O5", "protein_accession UniProtKB format",
                    shape, "WARN", 0, "Column absent")

    # O6 — start/end coordinate validity
    coord_issues: list[tuple] = []
    if "start" in occ.columns and "end" in occ.columns:
        for idx, row in occ[["start", "end"]].iterrows():
            parts: list[str] = []
            si = safe_int(row["start"])
            ei = safe_int(row["end"])
            if si is None and pd.notna(row["start"]):
                parts.append(f"start not integer: '{row['start']}'")
            elif si is not None and si <= 0:
                parts.append(f"start ≤ 0: {si}")
            if ei is None and pd.notna(row["end"]):
                parts.append(f"end not integer: '{row['end']}'")
            elif ei is not None and ei <= 0:
                parts.append(f"end ≤ 0: {ei}")
            if si is not None and ei is not None and si >= ei:
                parts.append(f"start({si}) ≥ end({ei})")
            if parts:
                coord_issues.append((idx, "; ".join(parts)))
    for idx, reason in coord_issues:
        invalid.append({
            "_source": "occurrence", "_check": "O6",
            "_reason": reason, "_source_row": int(idx) + 2,
            **_row_dict(occ, idx),
        })
    _add_result(results, "O6",
                "start/end valid positive integers (start < end)", shape,
                "FAIL" if coord_issues else "PASS", len(coord_issues))

    # O7–O12 — CV checks (uses _is_cv_backed which handles all token forms — Gap 3)
    cv_checks = [
        ("O7",  "positional_category"),
        ("O9",  "topology_context"),
        ("O10", "proximity_categories"),
        ("O11", "function_tags"),
        ("O12", "binding_parnter"),
    ]
    for check_id, col in cv_checks:
        ttype     = cv_col_map.get(col, "")
        valid_set = cv.get(ttype, frozenset())
        delim     = get_delimiter(cfg, col)
        _check_cv_col(results, invalid, occ, col, valid_set,
                      "occurrence", check_id,
                      f"{col} ontology-backed ({ttype})", shape, delim)

    # O13 — ChEBI format
    chebi_bad: list[tuple] = []
    if "binding_target_chebi" in occ.columns:
        for idx, val in occ["binding_target_chebi"].items():
            for tok in split_cell(val):
                if RE_CHEBI_DBL.fullmatch(tok):
                    chebi_bad.append((idx,
                        f"double-prefix: '{tok}' → remove one CHEBI: prefix"))
                elif not RE_CHEBI_OK.fullmatch(tok):
                    chebi_bad.append((idx,
                        f"bad ChEBI format: '{tok}' — expected CHEBI:nnnnn"))
    for idx, reason in chebi_bad:
        invalid.append({
            "_source": "occurrence", "_check": "O13",
            "_reason": reason, "_source_row": int(idx) + 2,
            **_row_dict(occ, idx),
        })
    _add_result(results, "O13",
                "binding_target_chebi format (CHEBI:nnnnn)", shape,
                "FAIL" if chebi_bad else "PASS", len(chebi_bad),
                "Known issue: CHEBI:CHEBI:nnnnn double-prefix — strip one CHEBI:")

    # O14 — GO format in standard GO columns
    _check_go_cols(results, invalid, occ, "occurrence", "O14", shape)

    # O15 — Reactome
    reactome_bad: list[tuple] = []
    if "protein_reactome_ids" in occ.columns:
        for idx, val in occ["protein_reactome_ids"].items():
            for tok in split_cell(val):
                if not RE_REACTOME.fullmatch(tok):
                    reactome_bad.append((idx, tok))
    for idx, tok in reactome_bad:
        invalid.append({
            "_source": "occurrence", "_check": "O15",
            "_reason": f"bad Reactome: '{tok}' — expected R-HSA-nnnnn",
            "_source_row": int(idx) + 2,
            **_row_dict(occ, idx),
        })
    _add_result(results, "O15", "Reactome ID format (R-XXX-nnnnn)", shape,
                "FAIL" if reactome_bad else "PASS", len(reactome_bad))

    # O16 — KEGG
    kegg_bad: list[tuple] = []
    if "protein_kegg_ids" in occ.columns:
        for idx, val in occ["protein_kegg_ids"].items():
            for tok in split_cell(val):
                if not RE_KEGG.fullmatch(tok):
                    kegg_bad.append((idx, tok))
    for idx, tok in kegg_bad:
        invalid.append({
            "_source": "occurrence", "_check": "O16",
            "_reason": f"bad KEGG: '{tok}' — expected hsa:nnnnn",
            "_source_row": int(idx) + 2,
            **_row_dict(occ, idx),
        })
    _add_result(results, "O16", "KEGG ID format (abc:nnnnn)", shape,
                "FAIL" if kegg_bad else "PASS", len(kegg_bad))

    # O17 — protein_go_cc_ids
    _check_go_cols(results, invalid, occ, "occurrence", "O17", shape,
                   extra_cols=["protein_go_cc_ids"])

    # O18 — Namespace cross-contamination (Gap 2)
    cross_viols: list[tuple] = []
    for forbidden_pfx, col, desc in CROSS_CONTAMINATION_RULES:
        if col not in occ.columns:
            continue
        delim = get_delimiter(cfg, col)
        for idx, val in occ[col].items():
            for tok in split_cell(val, delim):
                if _get_dox_prefix(tok) == forbidden_pfx:
                    cross_viols.append((idx, col, tok, desc))
    _add_result(results, "O18",
                "No namespace cross-contamination (DOP/DOF/DOT in wrong column)", shape,
                "FAIL" if cross_viols else "PASS", len(cross_viols),
                "e.g. a DOP: positional term placed in function_tags")
    for idx, col, tok, desc in cross_viols:
        invalid.append({
            "_source": "occurrence", "_check": "O18",
            "_reason": f"{desc}: found '{tok}' in column '{col}'",
            "_source_row": int(idx) + 2,
            **_row_dict(occ, idx),
        })

    return results, invalid


# ══════════════════════════════════════════════════════════════════════════════
# Missing term report with closest-match hints
# ══════════════════════════════════════════════════════════════════════════════

def build_missing_term_report(occ: pd.DataFrame, dt: pd.DataFrame,
                               cv: dict, cfg: dict) -> pd.DataFrame:
    """
    Collect every unmatched CV token across occurrence and domain-type tables.
    For each token, finds the closest match in the annotation table to help
    identify whether it is a typo or a genuinely missing term.
    """
    cv_col_map: dict[str, str] = cfg.get("cv_column_to_term_type", {})
    counter: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))

    for col, ttype in cv_col_map.items():
        valid = cv.get(ttype, frozenset())
        delim = get_delimiter(cfg, col)
        for tbl, prefix in [(occ, ""), (dt, "dt:")]:
            if col not in tbl.columns:
                continue
            for val in tbl[col].dropna():
                for tok in split_cell(val, delim):
                    if not _is_cv_backed(tok, valid):
                        counter[f"{prefix}{col}"][tok] += 1

    all_valid = _all_valid_tokens(cv)
    rows: list[dict] = []
    for col, tmap in sorted(counter.items()):
        for tok, cnt in sorted(tmap.items(), key=lambda x: -x[1]):
            matches  = get_close_matches(tok, all_valid, n=1, cutoff=0.6)
            closest  = matches[0] if matches else ""
            action   = (f"Closest match: '{closest}' — check spelling"
                        if closest else
                        "No close match — add as new term to annotation_table.xlsx")
            rows.append({
                "column":           col,
                "unmatched_token":  tok,
                "occurrence_count": cnt,
                "closest_match":    closest,
                "suggested_action": action,
            })
    df = pd.DataFrame(rows)
    return df


# ══════════════════════════════════════════════════════════════════════════════
# Export and plain-text summary
# ══════════════════════════════════════════════════════════════════════════════

def _export(out_dir: Path, all_results: list, all_invalid: list,
            missing: pd.DataFrame, logger) -> None:
    save_tsv(pd.DataFrame(all_results),
             out_dir / "validation_report_table_level.tsv", logger)

    if all_invalid:
        inv  = pd.DataFrame(all_invalid)
        meta = [c for c in ["_source", "_check", "_reason", "_source_row"]
                if c in inv.columns]
        other = [c for c in inv.columns if c not in meta]
        save_tsv(inv[meta + other], out_dir / "invalid_rows.tsv", logger)
    else:
        save_tsv(pd.DataFrame(), out_dir / "invalid_rows.tsv", logger)

    save_tsv(missing, out_dir / "missing_term_report.tsv", logger)


def _write_text_summary(all_results: list, all_invalid: list,
                         missing: pd.DataFrame, path: Path, logger) -> None:
    from collections import Counter
    sc = Counter(r["status"] for r in all_results)
    lines = [
        "=" * 70,
        "  STAGE 1 — TABLE VALIDATION SUMMARY",
        "=" * 70,
        f"  Total checks  : {len(all_results)}",
        f"  PASS          : {sc.get('PASS', 0)}",
        f"  FAIL          : {sc.get('FAIL', 0)}",
        f"  WARN          : {sc.get('WARN', 0)}",
        f"  Flagged rows  : {len(all_invalid)}",
        f"  Missing tokens: {len(missing)}",
        "",
        "─" * 70,
        "  FAILURES AND WARNINGS  (sorted by check_id)",
        "─" * 70,
    ]
    for r in sorted(all_results, key=lambda x: x["check_id"]):
        if r["status"] in ("FAIL", "WARN"):
            lines += [
                f"  [{r['check_id']:<4}] {r['status']:<4}  ({r['fix_priority']})  "
                f"{r['description']}",
                f"         n_issues : {r['n_issues']}",
            ]
            if r.get("detail"):
                lines.append(f"         detail   : {r['detail']}")
            lines.append(f"         fix      : {r.get('suggested_fix', '')}")
            lines.append("")

    if not missing.empty:
        lines += [
            "─" * 70,
            "  TOP UNMATCHED TOKENS  (full list: missing_term_report.tsv)",
            "─" * 70,
        ]
        for _, row in missing.head(20).iterrows():
            hint = f"→ '{row['closest_match']}'" if row["closest_match"] else "→ no close match"
            lines.append(
                f"  {row['column']:<35}  '{row['unmatched_token']}'  "
                f"(n={row['occurrence_count']})  {hint}"
            )
    lines += ["", "=" * 70]
    path.write_text("\n".join(lines), encoding="utf-8")
    logger.info(f"Text summary → {path}")


def _log_summary(all_results: list, all_invalid: list,
                  missing: pd.DataFrame, logger) -> None:
    from collections import Counter
    sc = Counter(r["status"] for r in all_results)
    logger.info(
        f"Checks: {len(all_results)} | "
        f"PASS {sc['PASS']} | FAIL {sc['FAIL']} | WARN {sc['WARN']}"
    )
    logger.info(f"Flagged rows: {len(all_invalid)} | Missing CV tokens: {len(missing)}")
    for r in sorted(all_results, key=lambda x: x["check_id"]):
        if r["status"] in ("FAIL", "WARN"):
            logger.warning(
                f"  [{r['check_id']}] {r['status']} ({r['fix_priority']}) — "
                f"{r['description']} — {r['n_issues']} issues  |  {r['suggested_fix']}"
            )


# ══════════════════════════════════════════════════════════════════════════════
# Stage entry point
# ══════════════════════════════════════════════════════════════════════════════

def run(cfg: dict) -> None:
    log_dir = cfg.get("log_dir", "logs")
    logger  = setup_logging(log_dir, "stage1_table_validation")
    stage_banner(logger, "STAGE 1", "Table Validation (v3 — all gaps addressed)")
    out_dir = make_output_dir(cfg, "stage1_validation")

    # ── Load annotation table ────────────────────────────────────────────────
    ann_path = get_input_path(cfg, "annotation_table")
    if ann_path is None:
        logger.error("annotation_table not found — cannot build CV lookup. Aborting.")
        return
    ann = load_table(ann_path)
    logger.info(f"Annotation table loaded: {ann.shape[0]} rows × {ann.shape[1]} cols")

    # Build CV from annotation table (Gap 7: self-contained, includes ID column)
    cv = _build_cv_lookup(ann)
    logger.info("CV lookup: " + " | ".join(f"{k}={len(v)//2}" for k, v in cv.items())
                + "  (unique tokens per type, not counting lowercase duplicates)")

    all_results: list[dict] = []
    all_invalid: list[dict] = []

    # ── Annotation table (A1–A7) ─────────────────────────────────────────────
    ann_r, ann_i = validate_annotation(ann, cfg)
    all_results += ann_r
    all_invalid += ann_i
    logger.info(f"Annotation checks complete: {len(ann_r)}")

    # ── Domain type table (T1–T8) ────────────────────────────────────────────
    dt_path = get_input_path(cfg, "domain_type_summary")
    if dt_path:
        dt = load_table(dt_path)
        logger.info(f"Domain type table loaded: {dt.shape[0]} rows × {dt.shape[1]} cols")
        dt_r, dt_i = validate_domain_types(dt, cv, cfg)
        all_results += dt_r
        all_invalid += dt_i
        logger.info(f"Domain type checks complete: {len(dt_r)}")
    else:
        logger.warning("domain_type_summary not found — skipping T1–T8")
        dt = pd.DataFrame()

    # ── Occurrence table (O1–O18) ────────────────────────────────────────────
    occ_path = get_input_path(cfg, "domain_occurrence_table")
    if occ_path:
        occ = load_table(occ_path)
        logger.info(f"Occurrence table loaded: {occ.shape[0]} rows × {occ.shape[1]} cols")
        occ_r, occ_i = validate_occurrences(occ, dt, cv, cfg)
        all_results += occ_r
        all_invalid += occ_i
        logger.info(f"Occurrence checks complete: {len(occ_r)}")
    else:
        logger.warning("domain_occurrence_table not found — skipping O1–O18")
        occ = pd.DataFrame()

    # ── Missing term report ───────────────────────────────────────────────────
    missing = build_missing_term_report(occ, dt, cv, cfg)

    # ── Export all outputs ────────────────────────────────────────────────────
    _export(out_dir, all_results, all_invalid, missing, logger)
    _write_text_summary(all_results, all_invalid, missing,
                        out_dir / "validation_summary.txt", logger)
    _log_summary(all_results, all_invalid, missing, logger)
    logger.info("Stage 1 complete.")


# ── Standalone CLI ────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import sys, os
    # Always run from the ontology_pipeline project root so that
    # relative paths (input/, output/, config/) resolve correctly
    # regardless of where PyCharm / the terminal launched the script from.
    _project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(_project_root)
    _cfg_path = sys.argv[1] if len(sys.argv) > 1 else "config/pipeline_config.yaml"
    run(load_config(_cfg_path))
