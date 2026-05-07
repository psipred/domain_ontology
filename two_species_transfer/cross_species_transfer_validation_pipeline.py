#!/usr/bin/env python3
"""
Cross-Species Ontology Transfer Validation Pipeline

Purpose
Validates two ontology-aligned protein domain datasets from the same pathway
but different species, measures how well the DOP/DOF/DOT ontology framework
transfers from the reference species (human) to a second species, and provides
a structured conservation analysis at the domain, function, and topology levels.

When an ortholog mapping file is provided the pipeline additionally computes
pairwise conservation scores across ortholog pairs (domain architecture,
function-tag, topology/positional annotation).

Designed to generalise to any number of species with minimal adaptation.

Usage
  # Without ortholog mapping
  python cross_species_transfer_validation_pipeline.py \
      --human   cAMP_human_mapped.tsv \
      --other   cAMP_mouse_mapped.tsv \
      --ann     annotation_table.xlsx \
      --out     output_cross_species \
      --species2 mouse

  # With ortholog mapping + plots
  python cross_species_transfer_validation_pipeline.py \
      --human   cAMP_human_mapped.tsv \
      --other   cAMP_mouse_mapped.tsv \
      --ann     annotation_table.xlsx \
      --ortho   ortholog_mapping.tsv \
      --out     output_cross_species \
      --species1 human \
      --species2 mouse \
      --plots

Outputs
  validation_report_human.tsv           per-check results for species 1
  validation_report_otherSpecies.tsv    per-check results for species 2
  invalid_rows.tsv                      all flagged rows from both species
  unknown_terms_report.tsv              tokens not found in any CV set
  missing_term_report.tsv               unmatched tokens + closest-match hints
  species_summary_metrics.tsv           per-species size and coverage stats
  transfer_metrics_cross_species.tsv    schema/vocab/rule portability scores
  overlap_metrics_cross_species.tsv     Jaccard and set overlap per category
  comparison_report_cross_species.tsv   side-by-side biological comparison
  ortholog_conservation_report.tsv      per-ortholog-pair conservation scores
                                        (only when --ortho is provided)
Optional plots (--plots)
  term_conservation_barplot.png
  jaccard_heatmap.png
  validation_summary_plot.png
  ortholog_conservation_plot.png        (only when --ortho is provided)
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from difflib import get_close_matches
from pathlib import Path
from typing import Any

import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

# ── Publication-quality global style (Nature/Cell) ────────────────────────────
plt.rcParams.update({
    "font.family":        "DejaVu Sans",
    "font.size":          8,
    "axes.titlesize":     9,
    "axes.labelsize":     8,
    "xtick.labelsize":    7,
    "ytick.labelsize":    7,
    "legend.fontsize":    7,
    "axes.linewidth":     0.7,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "xtick.major.width":  0.7,
    "ytick.major.width":  0.7,
    "xtick.major.size":   3,
    "ytick.major.size":   3,
    "xtick.direction":    "out",
    "ytick.direction":    "out",
    "legend.frameon":     False,
    "figure.dpi":         300,
    "savefig.dpi":        300,
    "savefig.bbox":       "tight",
    "axes.grid":          False,
    "figure.facecolor":   "white",
    "axes.facecolor":     "white",
    "pdf.fonttype":       42,
    "ps.fonttype":        42,
})

# ── Colorblind-safe palette (IBM / Okabe-Ito) ─────────────────────────────────
PAL = {
    "blue":      "#0077BB",
    "orange":    "#EE7733",
    "green":     "#009988",
    "red":       "#CC3311",
    "purple":    "#AA3377",
    "yellow":    "#CCBB44",
    "cyan":      "#33BBEE",
    "grey":      "#BBBBBB",
    "darkgrey":  "#555555",
    "lightgrey": "#EEEEEE",
    "pass":      "#44AA99",
    "warn":      "#DDAA33",
    "fail":      "#BB5566",
    "shared":    "#44AA99",
    "ref_only":  "#0077BB",
    "tgt_only":  "#EE7733",
    "unused":    "#BBBBBB",
}


def _panel_label(ax, letter: str, x: float = -0.10, y: float = 1.06) -> None:
    """Add a bold panel label (A, B, C…) in the upper-left corner of an axes."""
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="top", ha="left")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — Constants and configuration
# All regex patterns, column rules, and CV mappings are defined here so the
# script is easy to adapt for new species or extended schemas.
# ══════════════════════════════════════════════════════════════════════════════

# ── Identifier format patterns ────────────────────────────────────────────────
RE_ANN_ID    = re.compile(r"^(DOP|DOF|DOT):\d{6}$")
RE_GO        = re.compile(r"^GO:\d{7}$")
RE_REACTOME  = re.compile(r"^R-[A-Z]{3}-\d+$")
RE_KEGG      = re.compile(r"^[a-z]{2,5}:\d+$")
RE_CHEBI_OK  = re.compile(r"^CHEBI:\d+$")
RE_CHEBI_DBL = re.compile(r"^CHEBI:CHEBI:\d+$")
RE_UNIPROT   = re.compile(
    r"^[A-Z][0-9][A-Z0-9]{3}[0-9]([A-Z][0-9][A-Z0-9]{3}[0-9])?$"
)
RE_DT_ID = re.compile(
    r"^(CDD:cd\d+"
    r"|INTERPRO:IPR\d{6}"
    r"|NCBIFAM:TIGR\d{5}"
    r"|PFAM:PF\d{5}"
    r"|PRINTS:PR\d{5}"
    r"|PROFILE:PS\d{5}"
    r"|SMART:SM\d{5}"
    r"|SUPERFAMILY:SSF\d+"
    r"|PANTHER:PTHR\d+"
    r"|GENE3D:G3DSA:[\d.]+)$"
)
RE_ID_LABEL = re.compile(
    r"^((DOP|DOF|DOT|GO|IPR|CHEBI|PR|MOD|R-[A-Z]{3}):[^\s]+)\s+(.+)$"
)

# ── Relation-to-namespace rules ───────────────────────────────────────────────
# Each tuple: (column_name, required_namespace_prefix, rule_description)
# These rules are species-agnostic — the same DOP/DOF/DOT constraints apply
# regardless of which species the occurrence data comes from.
NAMESPACE_RULES: list[tuple[str, str, str]] = [
    ("positional_category",  "DOP",
     "positional_category must use DOP positional terms"),
    ("topology_context",     "DOP",
     "topology_context must use DOP topology terms"),
    ("proximity_categories", "DOP",
     "proximity_categories must use DOP proximity terms"),
    ("copy_number", "DOP",
     "copy_number must use DOP CopyNumber terms"),
    ("function_tags",        "DOF",
     "function_tags must use DOF functional terms"),
    ("binding_parnter",      "DOT",
     "binding_parnter must use DOT binding-target terms"),
    ("binding_partner",      "DOT",
     "binding_partner must use DOT binding-target terms"),
]

# ── CV column → annotation registry term_type ─────────────────────────────────
CV_COL_TO_TERM_TYPE: dict[str, str] = {
    "positional_category":  "positional_or_topology_term",
    "topology_context":     "positional_or_topology_term",
    "proximity_categories": "positional_or_topology_term",
    "copy_number": "CopyNumberCategory",
    "function_tags":        "functional_term",
    "binding_parnter":      "binding_target_term",
    "binding_partner":      "binding_target_term",
}

# ── Required columns — minimum schema for an occurrence table ─────────────────
# These columns must be present for the pipeline to run checks.
# Species-specific columns (e.g. kegg species prefix) are validated
# individually and do not block execution if absent.
REQUIRED_COLUMNS: list[str] = [
    "occurrence_id",
    "protein_accession",
    "domain_type_id",
    "start",
    "end",
    "positional_category",
    "copy_number",
    "topology_context",
    "function_tags",
]

# ── GO columns validated for ID format ────────────────────────────────────────
GO_COLS: list[str] = [
    "go_terms_all", "go_terms_mf", "go_terms_bp", "go_terms_cc",
    "protein_go_cc_ids",
]

# ── Conservation analysis column groups ───────────────────────────────────────
# Used in ortholog-pair comparison to measure per-category conservation.
CONSERVATION_GROUPS: dict[str, list[str]] = {
    "domain_types":        ["domain_type_id"],
    "function_tags":       ["function_tags"],
    "topology":            ["topology_context"],
    "positional":          ["positional_category"],
    "proximity":           ["proximity_categories"],
    "copy_number":         ["copy_number"],
    "binding_targets":     ["binding_parnter", "binding_partner"],
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — I/O helpers
# ══════════════════════════════════════════════════════════════════════════════

def load_table(path: Path) -> pd.DataFrame:
    """
    Load TSV or XLSX into a DataFrame with all columns typed as str.
    Strips whitespace from every string cell after loading.
    Raises FileNotFoundError with a clear message if the file is absent.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"Input file not found: {path}\n"
            f"Check that the path is correct and the file exists."
        )
    kw: dict[str, Any] = {"dtype": str}
    if path.suffix.lower() in (".xlsx", ".xls"):
        df = pd.read_excel(path, **kw)
    else:
        df = pd.read_csv(path, sep="\t", low_memory=False, **kw)
    for col in df.select_dtypes(include=["object", "string"]).columns:
        df[col] = df[col].str.strip()
    return df


def save_tsv(df: pd.DataFrame, path: Path, label: str = "") -> None:
    """Save DataFrame as TSV. Creates parent directories if needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    tag = f" [{label}]" if label else ""
    print(f"  ✓  Saved {len(df):>5} rows  →  {path}{tag}")


def print_banner(title: str) -> None:
    print(f"\n{'═' * 70}\n  {title}\n{'═' * 70}")


def print_section(title: str) -> None:
    print(f"\n{'─' * 70}\n  {title}\n{'─' * 70}")


def print_check(check_id: str, status: str, n: int, desc: str) -> None:
    icon   = {"PASS": "✓", "FAIL": "✗", "WARN": "⚠"}.get(status, "?")
    suffix = {"FAIL": "  ← ACTION REQUIRED",
              "WARN": "  ← review"}.get(status, "")
    print(f"  {icon} [{check_id}] {status:<4}  n={n:<5}  {desc}{suffix}")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — Token normalisation and CV utilities
# All token-handling logic is defined here so every section uses the same
# normalisation rules. This ensures consistency regardless of whether tokens
# appear as label-only, ID-only, combined "ID label", or "label:count" forms.
# ══════════════════════════════════════════════════════════════════════════════

def split_cell(val: Any, delimiter: str = ";") -> list[str]:
    """Split a multi-value cell string; return [] for null/empty values."""
    if pd.isna(val) or str(val).strip() == "":
        return []
    return [t.strip() for t in str(val).split(delimiter) if t.strip()]


def extract_id(token: str) -> str:
    """
    Extract only the ID portion from a combined 'ID label' token.
      'DOP:000003 Cterminal'  →  'DOP:000003'
      'Cterminal'             →  'Cterminal'
    """
    m = RE_ID_LABEL.match(token.strip())
    return m.group(1) if m else token.strip()


def extract_label(token: str) -> str:
    """
    Extract only the label portion from a combined 'ID label' token.
      'DOP:000003 Cterminal'  →  'Cterminal'
      'Cterminal'             →  'Cterminal'
    """
    m = RE_ID_LABEL.match(token.strip())
    return m.group(3) if m else token.strip()


def get_dox_namespace(token: str) -> str | None:
    """
    Return 'DOP', 'DOF', or 'DOT' if the ID part of a token carries one of
    those namespace prefixes. Returns None for tokens without such a prefix.
    """
    id_part = extract_id(token)
    for pfx in ("DOP", "DOF", "DOT"):
        if id_part.startswith(pfx + ":"):
            return pfx
    return None


def strip_count_suffix(token: str) -> str:
    """Strip ':N' count suffix: 'Cterminal:22' → 'Cterminal'."""
    return re.sub(r":\d+$", "", token.strip())


def is_cv_backed(token: str, valid_set: frozenset[str]) -> bool:
    """
    Return True if the token matches the valid set under any normalised form:
      - as-is                  (label-only: 'Cterminal')
      - lowercased             (case-insensitive)
      - ID part extracted      ('DOP:000003 Cterminal' → 'DOP:000003')
      - label part extracted   ('DOP:000003 Cterminal' → 'Cterminal')
      - count-suffix stripped  ('Cterminal:22' → 'Cterminal')
    All forms are also checked in lowercase.
    """
    t          = token.strip()
    id_part    = extract_id(t)
    label_part = extract_label(t)
    stripped   = strip_count_suffix(t)
    return any(v in valid_set for v in (
        t,          t.lower(),
        id_part,    id_part.lower(),
        label_part, label_part.lower(),
        stripped,   stripped.lower(),
    ))


def build_cv_lookup(ann: pd.DataFrame) -> dict[str, frozenset[str]]:
    """
    Build a controlled-vocabulary lookup from the annotation table.
    Returns {term_type: frozenset(valid_tokens)}.
    Valid tokens include: ID, Name, canonical_name, synonym — plus lowercase
    variants of each to enable case-insensitive matching.
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


def all_valid_label_tokens(cv: dict[str, frozenset[str]]) -> list[str]:
    """Flat list of non-lowercase valid tokens for difflib close-match hints."""
    return [v for s in cv.values() for v in s if not v.islower()]


def collect_column_tokens(df: pd.DataFrame, col: str) -> set[str]:
    """
    Collect all unique label-form tokens (lowercased) from a single column.
    Used for set-comparison and Jaccard calculations.
    """
    toks: set[str] = set()
    if col not in df.columns:
        return toks
    for val in df[col].dropna():
        for tok in split_cell(val):
            toks.add(extract_label(tok).lower())
    return toks


def collect_value_set(df: pd.DataFrame, col: str) -> set[str]:
    """
    Collect unique non-null entries from a single-value column (lowercased).
    Used for protein_accession and domain_type_id set comparisons.
    """
    if col not in df.columns:
        return set()
    return (
        set(df[col].dropna().str.strip().str.lower().unique()) - {"", "nan"}
    )


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — Validation engine (species-agnostic)
# Runs 18 checks (V1–V18) on one occurrence table.
# Identical logic is applied to both species tables so results are comparable.
# ══════════════════════════════════════════════════════════════════════════════

def _row_dict(df: pd.DataFrame, idx: Any) -> dict:
    """Serialise one DataFrame row to a flat dict, truncating long strings."""
    return {
        k: (str(v)[:300] if isinstance(v, str) else v)
        for k, v in df.loc[idx].to_dict().items()
    }


def validate_occurrence_table(
    df: pd.DataFrame,
    species_label: str,
    cv: dict[str, frozenset[str]],
) -> tuple[list[dict], list[dict], list[dict]]:
    """
    Run all validation checks on one species occurrence table.

    Parameters
    ----------
    df            : occurrence table loaded via load_table()
    species_label : short label for console output and report columns
    cv            : CV lookup dict from build_cv_lookup()

    Returns
    -------
    results  : list[dict] — one entry per check (the validation report)
    invalid  : list[dict] — one entry per flagged row
    unknown  : list[dict] — tokens not found in ANY CV set in the registry
    """
    results: list[dict] = []
    invalid: list[dict] = []
    unknown: list[dict] = []
    shape = df.shape

    # Flat set of ALL valid tokens across all term types, for unknown-term detection
    all_valid_flat: frozenset[str] = frozenset(
        tok for s in cv.values() for tok in s
    )

    def _add(check_id: str, desc: str, status: str, n: int,
             detail: str = "") -> None:
        results.append({
            "species":      species_label,
            "check_id":     check_id,
            "description":  desc,
            "table_rows":   shape[0],
            "table_cols":   shape[1],
            "status":       status,
            "n_issues":     n,
            "detail":       detail,
        })
        print_check(check_id, status, n, desc)

    def _flag(idx: Any, check_id: str, reason: str) -> None:
        invalid.append({
            "species":    species_label,
            "check_id":   check_id,
            "reason":     reason,
            "source_row": int(idx) + 2,  # 1-based + header → Excel row number
            **_row_dict(df, idx),
        })

    def _flag_unknown(idx: Any, col: str, tok: str, ttype: str) -> None:
        unknown.append({
            "species":    species_label,
            "column":     col,
            "token":      tok,
            "term_type":  ttype,
            "source_row": int(idx) + 2,
        })

    print(f"\n  Species: {species_label}  "
          f"({shape[0]:,} rows × {shape[1]} cols)")

    # ── Guard: cast expected string columns to str so .str accessor never
    #    receives float NaN (happens when an entire column is empty/NaN and
    #    pandas infers float64 dtype for that column).
    _STR_COLS = [
        "occurrence_id", "protein_accession", "domain_type_id",
        "positional_category", "copy_number", "topology_context",
        "proximity_categories", "function_tags", "binding_parnter",
        "binding_partner", "protein_go_bp_ids", "protein_go_mf_ids",
        "protein_go_cc_ids", "protein_reactome_ids", "protein_kegg_ids",
        "binding_target_chebi",
    ]
    for _c in _STR_COLS:
        if _c in df.columns:
            df[_c] = df[_c].fillna("").astype(str)

    # ── V1: Required columns ──────────────────────────────────────────────────
    missing_cols = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    _add("V1", "Required columns present",
         "FAIL" if missing_cols else "PASS",
         len(missing_cols), ", ".join(missing_cols))
    if missing_cols:
        # Cannot proceed with further checks if core columns are absent
        return results, invalid, unknown

    # Resolve occurrence_id column (handles both naming variants)
    oid_col = next(
        (c for c in ("occurrence_id", "domain_occurrence_id") if c in df.columns),
        "occurrence_id",
    )

    # ── V2: Blank occurrence_id ───────────────────────────────────────────────
    blank_oid = df[oid_col].isna() | (df[oid_col].str.strip() == "")
    n = int(blank_oid.sum())
    _add("V2", f"Blank {oid_col}", "FAIL" if n else "PASS", n)
    for idx in df[blank_oid].index:
        _flag(idx, "V2", f"blank {oid_col}")

    # ── V3a: Strict duplicate occurrence_id ──────────────────────────────────
    dup_strict = df.duplicated(subset=oid_col, keep=False) & ~blank_oid
    n = int(dup_strict.sum())
    _add("V3a", f"Duplicate {oid_col} (strict key)",
         "FAIL" if n else "PASS", n,
         "Genuine key violation — each occurrence_id must be unique")
    for idx in df[dup_strict].index:
        _flag(idx, "V3a",
              f"duplicate {oid_col}: '{df.loc[idx, oid_col]}'")

    # ── V3b: Full-row duplicates excluding batch_gene ────────────────────────
    key_cols = [c for c in df.columns if c != "batch_gene"]
    dup_row  = df.duplicated(subset=key_cols, keep=False) & ~blank_oid
    _add("V3b", "Full-row duplicate excl. batch_gene",
         "WARN" if dup_row.sum() else "PASS", int(dup_row.sum()),
         "Expected for batch_gene split rows; investigate if count is high")

    # ── V4: protein_accession UniProtKB format ────────────────────────────────
    if "protein_accession" in df.columns:
        filled = (df["protein_accession"].notna() &
                  (df["protein_accession"].str.strip() != ""))
        bad = filled & ~df["protein_accession"].str.strip().str.fullmatch(
            RE_UNIPROT.pattern).fillna(False)
        n = int(bad.sum())
        _add("V4", "protein_accession UniProtKB format",
             "FAIL" if n else "PASS", n,
             "Expected 6-char UniProtKB accession e.g. P06850")
        for idx in df[bad].index:
            _flag(idx, "V4",
                  f"bad accession: '{df.loc[idx,'protein_accession']}'")
    else:
        _add("V4", "protein_accession UniProtKB format",
             "WARN", 0, "Column absent")

    # ── V5: domain_type_id format ─────────────────────────────────────────────
    if "domain_type_id" in df.columns:
        filled = (df["domain_type_id"].notna() &
                  (df["domain_type_id"].str.strip() != ""))
        bad = filled & ~df["domain_type_id"].str.strip().str.fullmatch(
            RE_DT_ID.pattern).fillna(False)
        n = int(bad.sum())
        _add("V5", "domain_type_id format",
             "FAIL" if n else "PASS", n,
             "e.g. PFAM:PF00001 / INTERPRO:IPR000001 / CDD:cd15445")
        for idx in df[bad].index:
            _flag(idx, "V5",
                  f"bad domain_type_id: '{df.loc[idx,'domain_type_id']}'")

    # ── V6: start/end coordinate validity ────────────────────────────────────
    if "start" in df.columns and "end" in df.columns:
        coord_bad: list[tuple] = []
        for idx, row in df[["start", "end"]].iterrows():
            parts: list[str] = []
            try:
                si = int(float(row["start"]))
            except (ValueError, TypeError):
                si = None
                if pd.notna(row["start"]):
                    parts.append(f"start not integer: '{row['start']}'")
            try:
                ei = int(float(row["end"]))
            except (ValueError, TypeError):
                ei = None
                if pd.notna(row["end"]):
                    parts.append(f"end not integer: '{row['end']}'")
            if si is not None and si <= 0:
                parts.append(f"start ≤ 0: {si}")
            if ei is not None and ei <= 0:
                parts.append(f"end ≤ 0: {ei}")
            if si is not None and ei is not None and si >= ei:
                parts.append(f"start({si}) ≥ end({ei})")
            if parts:
                coord_bad.append((idx, "; ".join(parts)))
        n = len(coord_bad)
        _add("V6", "start/end valid positive integers (start < end)",
             "FAIL" if n else "PASS", n)
        for idx, reason in coord_bad:
            _flag(idx, "V6", reason)

    # ── V7: GO term format ────────────────────────────────────────────────────
    go_bad: list[tuple] = []
    for col in GO_COLS:
        if col not in df.columns:
            continue
        for idx, val in df[col].items():
            for tok in split_cell(val):
                if not RE_GO.fullmatch(tok):
                    go_bad.append((idx, f"{col}='{tok}'"))
    n = len(go_bad)
    _add("V7", "GO term format (GO:nnnnnnn)",
         "FAIL" if n else "PASS", n)
    for idx, reason in go_bad:
        _flag(idx, "V7", f"bad GO format: {reason}")

    # ── V8: binding_target_chebi format ──────────────────────────────────────
    chebi_bad: list[tuple] = []
    if "binding_target_chebi" in df.columns:
        for idx, val in df["binding_target_chebi"].items():
            for tok in split_cell(val):
                if RE_CHEBI_DBL.fullmatch(tok):
                    chebi_bad.append((idx,
                        f"double-prefix: '{tok}' → remove one 'CHEBI:'"))
                elif not RE_CHEBI_OK.fullmatch(tok):
                    chebi_bad.append((idx,
                        f"bad ChEBI: '{tok}' — expected CHEBI:nnnnn"))
    n = len(chebi_bad)
    _add("V8", "binding_target_chebi format (CHEBI:nnnnn)",
         "FAIL" if n else "PASS", n,
         "Common issue: CHEBI:CHEBI:nnnnn double-prefix — strip one CHEBI:")
    for idx, reason in chebi_bad:
        _flag(idx, "V8", reason)

    # ── V9: Reactome ID format ────────────────────────────────────────────────
    if "protein_reactome_ids" in df.columns:
        r_bad: list[tuple] = []
        for idx, val in df["protein_reactome_ids"].items():
            for tok in split_cell(val):
                if not RE_REACTOME.fullmatch(tok):
                    r_bad.append((idx, tok))
        n = len(r_bad)
        _add("V9", "Reactome ID format (R-XXX-nnnnn)",
             "FAIL" if n else "PASS", n)
        for idx, tok in r_bad:
            _flag(idx, "V9", f"bad Reactome: '{tok}'")
    else:
        _add("V9", "Reactome ID format", "WARN", 0, "Column absent")

    # ── V10: KEGG ID format ───────────────────────────────────────────────────
    if "protein_kegg_ids" in df.columns:
        k_bad: list[tuple] = []
        for idx, val in df["protein_kegg_ids"].items():
            for tok in split_cell(val):
                if not RE_KEGG.fullmatch(tok):
                    k_bad.append((idx, tok))
        n = len(k_bad)
        _add("V10", "KEGG ID format (abc:nnnnn)",
             "FAIL" if n else "PASS", n,
             "Note: species prefix changes across species e.g. hsa: vs mmu:")
        for idx, tok in k_bad:
            _flag(idx, "V10", f"bad KEGG: '{tok}'")
    else:
        _add("V10", "KEGG ID format", "WARN", 0, "Column absent")

    # ── V11–V16: CV-backed column checks ─────────────────────────────────────
    cv_checks: list[tuple[str, str, str]] = [
        ("V11", "positional_category",  "positional_or_topology_term"),
        ("V12", "copy_number", "CopyNumberCategory"),
        ("V13", "topology_context",     "positional_or_topology_term"),
        ("V14", "proximity_categories", "positional_or_topology_term"),
        ("V15", "function_tags",        "functional_term"),
        ("V16", "binding_parnter",      "binding_target_term"),
    ]
    for check_id, col, ttype in cv_checks:
        if col not in df.columns:
            _add(check_id, f"{col} CV-backed ({ttype})",
                 "WARN", 0, f"Column '{col}' absent")
            continue
        valid_set = cv.get(ttype, frozenset())
        bad_rows: dict[Any, list[str]] = defaultdict(list)
        for idx, val in df[col].items():
            for tok in split_cell(val):
                if not is_cv_backed(tok, valid_set):
                    bad_rows[idx].append(tok)
                    # Flag tokens absent from the entire registry
                    if not is_cv_backed(tok, all_valid_flat):
                        _flag_unknown(idx, col, tok, ttype)
        n_toks = sum(len(v) for v in bad_rows.values())
        _add(check_id, f"{col} CV-backed ({ttype})",
             "FAIL" if n_toks else "PASS", n_toks,
             f"{len(bad_rows)} rows contain unmatched tokens")
        for idx, toks in bad_rows.items():
            _flag(idx, check_id,
                  f"unmatched tokens in '{col}': " + "; ".join(toks))

    # ── V17: Invalid namespace prefix in wrong column ─────────────────────────
    ns_viols: list[tuple] = []
    for col, required_ns, rule_desc in NAMESPACE_RULES:
        if col not in df.columns:
            continue
        for idx, val in df[col].items():
            for tok in split_cell(val):
                actual_ns = get_dox_namespace(tok)
                if actual_ns is not None and actual_ns != required_ns:
                    ns_viols.append(
                        (idx, col, tok, actual_ns, required_ns, rule_desc)
                    )
    n = len(ns_viols)
    _add("V17",
         "Relation-to-namespace rules (DOP/DOF/DOT in correct column)",
         "FAIL" if n else "PASS", n,
         "e.g. a DOF functional term placed in positional_category column")
    for idx, col, tok, actual, req, desc in ns_viols:
        _flag(idx, "V17",
              f"{desc}: found '{tok}' (ns={actual}) in column '{col}'")

    # ── V18: Invalid namespace prefix format in annotation-linked columns ─────
    # Checks that any DOP/DOF/DOT IDs present conform to RE_ANN_ID format.
    format_viols: list[tuple] = []
    for col in CV_COL_TO_TERM_TYPE:
        if col not in df.columns:
            continue
        for idx, val in df[col].items():
            for tok in split_cell(val):
                id_part = extract_id(tok)
                if (id_part.startswith(("DOP:", "DOF:", "DOT:")) and
                        not RE_ANN_ID.fullmatch(id_part)):
                    format_viols.append((idx, col, id_part))
    n = len(format_viols)
    _add("V18",
         "DOP/DOF/DOT ID format (PREFIX:NNNNNN) in CV columns",
         "FAIL" if n else "PASS", n,
         "IDs must match DOP:000001 / DOF:000001 / DOT:000001 pattern")
    for idx, col, bad_id in format_viols:
        _flag(idx, "V18",
              f"malformed ontology ID '{bad_id}' in column '{col}'")

    return results, invalid, unknown

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — Missing term report
# ══════════════════════════════════════════════════════════════════════════════

def build_missing_term_report(
    dfs: dict[str, pd.DataFrame],
    cv: dict[str, frozenset[str]],
) -> pd.DataFrame:
    """
    Collect every CV token from all species tables that is not in the annotation
    registry. Provides closest-match hints to distinguish typos from genuinely
    missing terms that need to be added to annotation_table.xlsx.

    Parameters
    ----------
    dfs : {species_label: DataFrame}
    cv  : CV lookup from build_cv_lookup()

    Returns
    -------
    DataFrame with columns:
      species, column, unmatched_token, occurrence_count,
      closest_match, suggested_action
    """
    counter: dict[tuple[str, str, str], int] = defaultdict(int)
    for label, df in dfs.items():
        for col, ttype in CV_COL_TO_TERM_TYPE.items():
            if col not in df.columns:
                continue
            valid = cv.get(ttype, frozenset())
            for val in df[col].dropna():
                for tok in split_cell(val):
                    if not is_cv_backed(tok, valid):
                        counter[(label, col, tok)] += 1

    all_valid = all_valid_label_tokens(cv)
    rows: list[dict] = []
    for (species, col, tok), cnt in sorted(
            counter.items(), key=lambda x: (-x[1], x[0][0], x[0][1])):
        matches = get_close_matches(tok, all_valid, n=1, cutoff=0.6)
        closest = matches[0] if matches else ""
        action  = (
            f"Closest match: '{closest}' — check spelling"
            if closest else
            "No close match — add new term to annotation_table.xlsx"
        )
        rows.append({
            "species":          species,
            "column":           col,
            "unmatched_token":  tok,
            "occurrence_count": cnt,
            "closest_match":    closest,
            "suggested_action": action,
        })
    return pd.DataFrame(rows)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — Per-species summary metrics
# ══════════════════════════════════════════════════════════════════════════════

def compute_species_metrics(
    df: pd.DataFrame,
    species_label: str,
    cv: dict[str, frozenset[str]],
) -> dict:
    """
    Compute per-species biological and ontology coverage statistics.

    Returns a flat dict for one row of species_summary_metrics.tsv.
    """
    m: dict[str, Any] = {"species": species_label}

    # Basic counts
    m["n_rows"] = len(df)
    m["n_unique_proteins"] = (
        df["protein_accession"].nunique()
        if "protein_accession" in df.columns else 0
    )
    m["n_unique_domain_types"] = (
        df["domain_type_id"].nunique()
        if "domain_type_id" in df.columns else 0
    )

    # Schema coverage: fraction of required columns present
    present = sum(1 for c in REQUIRED_COLUMNS if c in df.columns)
    m["schema_coverage_pct"] = round(
        100 * present / len(REQUIRED_COLUMNS), 1
    )

    # CV column fill rate and unique term counts
    for col in CV_COL_TO_TERM_TYPE:
        if col in df.columns:
            filled = df[col].notna() & (df[col].str.strip() != "")
            m[f"pct_filled_{col}"] = round(
                100 * filled.sum() / max(len(df), 1), 1
            )
            m[f"n_unique_terms_{col}"] = len(collect_column_tokens(df, col))
        else:
            m[f"pct_filled_{col}"]     = None
            m[f"n_unique_terms_{col}"] = 0

    # CV-backed rate per column (fraction of tokens that match the registry)
    for col, ttype in CV_COL_TO_TERM_TYPE.items():
        if col not in df.columns:
            m[f"cv_backed_rate_{col}"] = None
            continue
        valid = cv.get(ttype, frozenset())
        total = backed = 0
        for val in df[col].dropna():
            for tok in split_cell(val):
                total += 1
                if is_cv_backed(tok, valid):
                    backed += 1
        m[f"cv_backed_rate_{col}"] = round(
            100 * backed / max(total, 1), 1
        )

    # Domain type source distribution
    if "domain_type_id" in df.columns:
        sources: dict[str, int] = defaultdict(int)
        for v in df["domain_type_id"].dropna():
            prefix = v.split(":")[0] if ":" in v else "UNKNOWN"
            sources[prefix] += 1
        for src, cnt in sorted(sources.items()):
            m[f"n_domain_type_{src}"] = cnt

    return m


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7 — Transfer metrics (cross-species)
# ══════════════════════════════════════════════════════════════════════════════

def compute_transfer_metrics(
    df_ref: pd.DataFrame,
    df_tgt: pd.DataFrame,
    cv: dict[str, frozenset[str]],
    label_ref: str = "human",
    label_tgt: str = "otherSpecies",
) -> dict:
    """
    Measure how well the human ontology framework transfers to the second species.

    Key metrics
    ───────────
    schema_reuse_rate          fraction of reference columns present in target
    vocab_reuse_rate           fraction of reference CV tokens re-used in target
    rule_portability_rate      fraction of namespace rules testable in target
    mapping_completeness       fraction of target CV tokens covered by registry
    ontology_term_conservation fraction of registry terms used by both species
    species_specific_rate      fraction of target tokens not found in reference
    """
    m: dict[str, Any] = {
        "reference_species": label_ref,
        "target_species":    label_tgt,
    }

    # ── Schema reuse ──────────────────────────────────────────────────────────
    ref_cols = set(df_ref.columns)
    tgt_cols = set(df_tgt.columns)
    shared   = ref_cols & tgt_cols
    m["schema_reuse_rate"]     = round(len(shared) / max(len(ref_cols), 1), 4)
    m["n_shared_columns"]      = len(shared)
    m["n_ref_only_columns"]    = len(ref_cols - tgt_cols)
    m["n_target_only_columns"] = len(tgt_cols - ref_cols)

    # ── CV token sets ─────────────────────────────────────────────────────────
    def _all_cv_tokens(df: pd.DataFrame) -> set[str]:
        toks: set[str] = set()
        for col in CV_COL_TO_TERM_TYPE:
            toks |= collect_column_tokens(df, col)
        return toks

    ref_toks = _all_cv_tokens(df_ref)
    tgt_toks = _all_cv_tokens(df_tgt)
    shared_toks = ref_toks & tgt_toks

    m["vocab_reuse_rate"]     = round(len(shared_toks) / max(len(ref_toks), 1), 4)
    m["n_ref_cv_tokens"]      = len(ref_toks)
    m["n_target_cv_tokens"]   = len(tgt_toks)
    m["n_shared_cv_tokens"]   = len(shared_toks)
    m["n_species_specific_tokens"] = len(tgt_toks - ref_toks)

    # ── Rule portability ──────────────────────────────────────────────────────
    testable = sum(
        1 for col, _, _ in NAMESPACE_RULES if col in df_tgt.columns
    )
    m["rule_portability_rate"] = round(
        testable / max(len(NAMESPACE_RULES), 1), 4
    )
    m["n_testable_rules"]  = testable
    m["n_total_rules"]     = len(NAMESPACE_RULES)

    # ── Mapping completeness ──────────────────────────────────────────────────
    total_tgt = covered_tgt = 0
    for col, ttype in CV_COL_TO_TERM_TYPE.items():
        if col not in df_tgt.columns:
            continue
        valid = cv.get(ttype, frozenset())
        for val in df_tgt[col].dropna():
            for tok in split_cell(val):
                total_tgt += 1
                if is_cv_backed(tok, valid):
                    covered_tgt += 1
    m["mapping_completeness"]     = round(covered_tgt / max(total_tgt, 1), 4)
    m["n_target_tokens_total"]    = total_tgt
    m["n_target_tokens_covered"]  = covered_tgt

    # ── Ontology term conservation ────────────────────────────────────────────
    # Registry terms (label form, lowercased, non-duplicate)
    registry: set[str] = set()
    for s in cv.values():
        for tok in s:
            if not tok.islower():
                registry.add(tok.lower())

    used_by_both   = ref_toks & tgt_toks & registry
    ref_only_reg   = (ref_toks & registry) - tgt_toks
    tgt_only_reg   = (tgt_toks & registry) - ref_toks
    new_required   = tgt_toks - registry
    unused_in_reg  = registry - ref_toks - tgt_toks

    m["ontology_term_conservation_rate"] = round(
        len(used_by_both) / max(len(ref_toks & registry), 1), 4
    )
    m["species_specific_annotation_rate"] = round(
        len(tgt_toks - ref_toks) / max(len(tgt_toks), 1), 4
    )
    m["n_terms_conserved"]           = len(used_by_both)
    m["n_terms_ref_only"]            = len(ref_only_reg)
    m["n_terms_target_only"]         = len(tgt_only_reg)
    m["n_new_terms_required"]        = len(new_required)
    m["n_registry_terms_unused"]     = len(unused_in_reg)

    return m


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 8 — Overlap and Jaccard metrics
# ══════════════════════════════════════════════════════════════════════════════

def jaccard(a: set, b: set) -> float:
    """Jaccard similarity: |A ∩ B| / |A ∪ B|. Returns 0.0 if both sets empty."""
    union = a | b
    return round(len(a & b) / max(len(union), 1), 4)


def compute_overlap_metrics(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    label1: str = "human",
    label2: str = "otherSpecies",
) -> pd.DataFrame:
    """
    Compute set overlap and Jaccard similarity across all biological and
    ontology categories.

    Returns a DataFrame with one row per category:
      category, <label1>_size, <label2>_size, shared, unique_to_<label1>,
      unique_to_<label2>, jaccard, shared_items (first 10), ...
    """
    categories: list[tuple[str, set, set]] = [
        ("domain_types",
         collect_value_set(df1, "domain_type_id"),
         collect_value_set(df2, "domain_type_id")),
        ("proteins",
         collect_value_set(df1, "protein_accession"),
         collect_value_set(df2, "protein_accession")),
        ("function_tags",
         collect_column_tokens(df1, "function_tags"),
         collect_column_tokens(df2, "function_tags")),
        ("binding_targets",
         collect_column_tokens(df1, "binding_partner") or collect_column_tokens(df1, "binding_parnter"),
         collect_column_tokens(df2, "binding_partner") or collect_column_tokens(df2, "binding_parnter")),
        ("topology_context",
         collect_column_tokens(df1, "topology_context"),
         collect_column_tokens(df2, "topology_context")),
        ("positional_category",
         collect_column_tokens(df1, "positional_category"),
         collect_column_tokens(df2, "positional_category")),
        ("proximity_categories",
         collect_column_tokens(df1, "proximity_categories"),
         collect_column_tokens(df2, "proximity_categories")),
        ("copy_number",
         collect_column_tokens(df1, "copy_number"),
         collect_column_tokens(df2, "copy_number")),
        ("go_terms_mf",
         collect_column_tokens(df1, "go_terms_mf"),
         collect_column_tokens(df2, "go_terms_mf")),
        ("go_terms_bp",
         collect_column_tokens(df1, "go_terms_bp"),
         collect_column_tokens(df2, "go_terms_bp")),
        ("go_terms_cc",
         collect_column_tokens(df1, "go_terms_cc"),
         collect_column_tokens(df2, "go_terms_cc")),
    ]

    rows: list[dict] = []
    for cat, s1, s2 in categories:
        rows.append({
            "category":              cat,
            f"{label1}_size":        len(s1),
            f"{label2}_size":        len(s2),
            "shared":                len(s1 & s2),
            f"unique_to_{label1}":   len(s1 - s2),
            f"unique_to_{label2}":   len(s2 - s1),
            "jaccard":               jaccard(s1, s2),
            "shared_items":          "; ".join(sorted(s1 & s2)[:10]),
            f"{label1}_unique_items":"; ".join(sorted(s1 - s2)[:10]),
            f"{label2}_unique_items":"; ".join(sorted(s2 - s1)[:10]),
        })
    return pd.DataFrame(rows)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 9 — Ortholog conservation analysis
# Runs only when an ortholog_mapping.tsv file is provided.
# For each ortholog pair the pipeline computes per-category conservation:
# whether the two proteins carry the same domain types, function tags,
# topology context, positional category, proximity categories, etc.
# ══════════════════════════════════════════════════════════════════════════════

def load_ortholog_mapping(path: Path) -> pd.DataFrame:
    """
    Load ortholog mapping TSV.

    Expected columns (tab-separated, header row required):
      human_accession       UniProtKB accession of the human protein
      other_accession       UniProtKB accession of the orthologue
      gene_name             (optional) gene symbol
      ortholog_type         (optional) e.g. one2one, one2many, many2many

    Returns a clean DataFrame with columns normalised to lowercase.
    """
    df = load_table(path)
    # Normalise column names
    df.columns = df.columns.str.strip().str.lower().str.replace(" ", "_")
    # Accept common column name variants
    rename_map = {
        "human_protein": "human_accession",
        "human_protein_accession": "human_accession",
        "species1_accession": "human_accession",
        "other_protein": "other_accession",
        "other_protein_accession": "other_accession",
        "species2_accession": "other_accession",
        "ortholog_accession": "other_accession",
    }
    df = df.rename(columns={k: v for k, v in rename_map.items()
                             if k in df.columns})
    required = {"human_accession", "other_accession"}
    missing  = required - set(df.columns)
    if missing:
        raise ValueError(
            f"ortholog_mapping.tsv is missing required columns: {missing}\n"
            f"Found columns: {list(df.columns)}\n"
            f"Required: human_accession, other_accession"
        )
    return df


def _protein_profile(
    df: pd.DataFrame,
    protein_accession: str,
) -> dict[str, set[str]]:
    """
    Build a per-protein annotation profile: one set of tokens per category.
    Used to compare annotation between orthologous protein pairs.

    Parameters
    ----------
    df                : species occurrence table
    protein_accession : UniProtKB accession to extract

    Returns
    -------
    dict mapping category_name → set of unique label tokens for that protein
    """
    if "protein_accession" not in df.columns:
        return {cat: set() for cat in CONSERVATION_GROUPS}
    mask = df["protein_accession"].str.upper() == protein_accession.upper()
    sub  = df[mask]
    profile: dict[str, set[str]] = {}
    for cat, cols in CONSERVATION_GROUPS.items():
        toks: set[str] = set()
        for col in cols:
            if col in sub.columns:
                for val in sub[col].dropna():
                    for tok in split_cell(val):
                        toks.add(extract_label(tok).lower())
        profile[cat] = toks
    return profile


def compute_ortholog_conservation(
    df_human: pd.DataFrame,
    df_other: pd.DataFrame,
    ortho_map: pd.DataFrame,
    label_ref: str = "human",
    label_tgt: str = "otherSpecies",
) -> pd.DataFrame:
    """
    Compute per-ortholog-pair conservation scores for all annotation categories.

    For each pair (human_accession, other_accession):
      - Retrieve the annotation profile for each protein
      - Compute Jaccard similarity per category
      - Compute an overall mean conservation score
      - Flag pairs where one or both proteins have no occurrences

    Also computes aggregate summary statistics:
      - ortholog_coverage: fraction of ortholog pairs where both proteins
        have at least one occurrence row
      - per-category mean conservation across all complete pairs

    Returns a DataFrame with one row per ortholog pair.
    """
    rows: list[dict] = []

    for _, pair in ortho_map.iterrows():
        h_acc = str(pair["human_accession"]).strip()
        o_acc = str(pair["other_accession"]).strip()
        gene  = str(pair.get("gene_name", "")).strip()
        otype = str(pair.get("ortholog_type", "")).strip()

        h_profile = _protein_profile(df_human, h_acc)
        o_profile = _protein_profile(df_other, o_acc)

        h_has_data = any(len(v) > 0 for v in h_profile.values())
        o_has_data = any(len(v) > 0 for v in o_profile.values())

        row: dict[str, Any] = {
            "human_accession":   h_acc,
            "other_accession":   o_acc,
            "gene_name":         gene,
            "ortholog_type":     otype,
            f"{label_ref}_has_data": h_has_data,
            f"{label_tgt}_has_data": o_has_data,
            "pair_complete":     h_has_data and o_has_data,
        }

        # Per-category Jaccard scores
        cat_scores: list[float] = []
        for cat in CONSERVATION_GROUPS:
            h_set = h_profile.get(cat, set())
            o_set = o_profile.get(cat, set())
            j     = jaccard(h_set, o_set)
            row[f"jaccard_{cat}"]       = j
            row[f"{label_ref}_{cat}"]   = "; ".join(sorted(h_set)[:8])
            row[f"{label_tgt}_{cat}"]   = "; ".join(sorted(o_set)[:8])
            row[f"shared_{cat}"]        = "; ".join(
                sorted(h_set & o_set)[:8])
            if h_has_data and o_has_data:
                cat_scores.append(j)

        row["mean_conservation_score"] = (
            round(sum(cat_scores) / len(cat_scores), 4)
            if cat_scores else None
        )
        rows.append(row)

    df_report = pd.DataFrame(rows)

    # ── Print aggregate summary ───────────────────────────────────────────────
    total_pairs    = len(df_report)
    complete_pairs = int(df_report["pair_complete"].sum())
    coverage       = round(complete_pairs / max(total_pairs, 1), 4)

    print(f"  Ortholog pairs total     : {total_pairs:,}")
    print(f"  Pairs with data in both  : {complete_pairs:,}")
    print(f"  Ortholog coverage        : {coverage:.1%}")

    complete = df_report[df_report["pair_complete"]]
    if len(complete) > 0:
        print(f"  Mean conservation scores (complete pairs):")
        for cat in CONSERVATION_GROUPS:
            col = f"jaccard_{cat}"
            if col in complete.columns:
                mean_j = complete[col].mean()
                print(f"    {cat:<35} {mean_j:.3f}")

    return df_report


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 10 — Comparison report
# ══════════════════════════════════════════════════════════════════════════════

def build_comparison_report(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    label1: str,
    label2: str,
    cv: dict[str, frozenset[str]],
) -> pd.DataFrame:
    """
    Build a side-by-side biological and ontology comparison of two species tables.
    Each row is a named metric with values for both species and a delta column.
    """

    def _cv_rate(df: pd.DataFrame, col: str, ttype: str) -> float:
        if col not in df.columns:
            return 0.0
        valid = cv.get(ttype, frozenset())
        total = backed = 0
        for val in df[col].dropna():
            for tok in split_cell(val):
                total += 1
                if is_cv_backed(tok, valid):
                    backed += 1
        return round(100 * backed / max(total, 1), 1)

    def _n_unique_toks(df: pd.DataFrame, col: str) -> int:
        return len(collect_column_tokens(df, col))

    rows: list[tuple[str, Any, Any]] = [
        ("Total occurrence rows",        len(df1),                    len(df2)),
        ("Unique proteins",
         df1["protein_accession"].nunique() if "protein_accession" in df1.columns else 0,
         df2["protein_accession"].nunique() if "protein_accession" in df2.columns else 0),
        ("Unique domain types",
         df1["domain_type_id"].nunique() if "domain_type_id" in df1.columns else 0,
         df2["domain_type_id"].nunique() if "domain_type_id" in df2.columns else 0),
        ("Schema columns present",       len(df1.columns),            len(df2.columns)),
        ("Unique function tags",
         _n_unique_toks(df1, "function_tags"),
         _n_unique_toks(df2, "function_tags")),
        ("Unique topology terms",
         _n_unique_toks(df1, "topology_context"),
         _n_unique_toks(df2, "topology_context")),
        ("Unique positional terms",
         _n_unique_toks(df1, "positional_category"),
         _n_unique_toks(df2, "positional_category")),
        ("Unique proximity terms",
         _n_unique_toks(df1, "proximity_categories"),
         _n_unique_toks(df2, "proximity_categories")),
        ("Unique binding targets",
         _n_unique_toks(df1, "binding_parnter"),
         _n_unique_toks(df2, "binding_parnter")),
        ("Unique copy-number terms",
         _n_unique_toks(df1, "copy_number"),
         _n_unique_toks(df2, "copy_number")),
        ("CV-backed rate: function_tags (%)",
         _cv_rate(df1, "function_tags", "functional_term"),
         _cv_rate(df2, "function_tags", "functional_term")),
        ("CV-backed rate: topology_context (%)",
         _cv_rate(df1, "topology_context", "positional_or_topology_term"),
         _cv_rate(df2, "topology_context", "positional_or_topology_term")),
        ("CV-backed rate: positional_category (%)",
         _cv_rate(df1, "positional_category", "positional_or_topology_term"),
         _cv_rate(df2, "positional_category", "positional_or_topology_term")),
        ("CV-backed rate: binding_parnter (%)",
         _cv_rate(df1, "binding_parnter", "binding_target_term"),
         _cv_rate(df2, "binding_parnter", "binding_target_term")),
        ("CV-backed rate: copy_number (%)",
         _cv_rate(df1, "copy_number", "CopyNumberCategory"),
         _cv_rate(df2, "copy_number", "CopyNumberCategory")),
    ]

    out: list[dict] = []
    for name, v1, v2 in rows:
        try:
            diff = round(float(v2) - float(v1), 2)
        except (TypeError, ValueError):
            diff = None
        out.append({
            "metric":          name,
            label1:            v1,
            label2:            v2,
            f"delta_{label2}_minus_{label1}": diff,
        })
    return pd.DataFrame(out)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 11 — Publication-quality figures  (7 figures, 300 DPI, Nature style)
#
# Figure inventory
# ────────────────
#   Fig 1  validation_summary_plot.png
#          A: PASS/WARN/FAIL stacked bars per species
#          B: Per-check status grid heatmap
#          C: FAIL/WARN detail table
#
#   Fig 2  species_comparison_plot.png
#          A: Dataset sizes  B: Function-tag composition
#          C: CV coverage rates  D: Topology context distributions
#          E: Transfer scores  F: Summary statistics table
#
#   Fig 3  jaccard_overlap_plot.png
#          A: Lollipop Jaccard per category
#          B: Set-size breakdown (shared / unique per species)
#
#   Fig 4  transfer_metrics_plot.png
#          A: Spider chart of 4 key portability metrics
#          B: Ontology term conservation breakdown bars
#
#   Fig 5  domain_conservation_plot.png
#          A: Top conserved and species-specific domain types
#          B: Occurrence-count comparison scatter (shared vs unique domains)
#
#   Fig 6  function_profile_plot.png
#          A: Stacked function-tag composition bars per species
#          B: Per-tag conservation lollipop chart
#
#   Fig 7  ortholog_conservation_plot.png  (only when --ortho provided)
#          A: Mean Jaccard per annotation category
#          B: Per-pair conservation score distribution
# ══════════════════════════════════════════════════════════════════════════════


# ─────────────────────────────────────────────────────────────────────────────
# Fig 1 — Validation summary
# ─────────────────────────────────────────────────────────────────────────────

def plot_validation_summary(
    results1: list[dict],
    results2: list[dict],
    out_path: Path,
    label1: str = "human",
    label2: str = "otherSpecies",
) -> None:
    """
    Three-panel validation summary.
      A — PASS/WARN/FAIL stacked bars per species
      B — Per-check status grid (colour-coded)
      C — FAIL/WARN detail table
    """
    from collections import Counter
    sp = {"PASS": PAL["pass"], "WARN": PAL["warn"], "FAIL": PAL["fail"]}

    fig = plt.figure(figsize=(12, 5.5))
    gs  = gridspec.GridSpec(1, 3, figure=fig, wspace=0.38,
                            left=0.06, right=0.97, top=0.88, bottom=0.08,
                            width_ratios=[2, 3, 3])
    ax_b = fig.add_subplot(gs[0])
    ax_g = fig.add_subplot(gs[1])
    ax_t = fig.add_subplot(gs[2])

    # A — stacked bars
    for yi, (lbl, res) in enumerate([(label2, results2), (label1, results1)]):
        cnt = Counter(r["status"] for r in res)
        tot = sum(cnt.values()); left = 0
        for st in ["PASS", "WARN", "FAIL"]:
            v = cnt.get(st, 0)
            if v == 0: continue
            ax_b.barh(yi, v, height=0.5, left=left,
                      color=sp[st], edgecolor="white", lw=0.7)
            col = "white" if st != "WARN" else "#2c3e50"
            ax_b.text(left + v/2, yi,
                      f"{v}\n({100*v/max(tot,1):.0f}%)",
                      ha="center", va="center", fontsize=6.5,
                      fontweight="bold", color=col)
            left += v
    ax_b.set_yticks([0, 1])
    ax_b.set_yticklabels([label2, label1], fontsize=8)
    ax_b.set_xlabel("Number of checks")
    ax_b.set_title("Validation Summary", fontweight="bold", pad=6)
    ax_b.legend(handles=[mpatches.Patch(facecolor=sp[s], label=s)
                          for s in ["PASS", "WARN", "FAIL"]],
                fontsize=6.5, loc="lower right", handlelength=1.2)
    _panel_label(ax_b, "A")

    # B — per-check status grid
    all_cids = sorted({r["check_id"] for r in results1 + results2})
    snum     = {"PASS": 2, "WARN": 1, "FAIL": 0}
    gd       = np.full((2, len(all_cids)), np.nan)
    for pi, res in enumerate([results1, results2]):
        lk = {r["check_id"]: r["status"] for r in res}
        for ci, cid in enumerate(all_cids):
            if cid in lk: gd[pi, ci] = snum.get(lk[cid], -1)

    cmap3 = LinearSegmentedColormap.from_list(
        "pf", [PAL["fail"], PAL["warn"], PAL["pass"]], N=3)
    ax_g.imshow(gd, aspect="auto", cmap=cmap3,
                vmin=-0.5, vmax=2.5, interpolation="nearest")
    ax_g.set_xticks(range(len(all_cids)))
    ax_g.set_xticklabels(all_cids, rotation=60, ha="right", fontsize=6.5)
    ax_g.set_yticks([0, 1])
    ax_g.set_yticklabels([label1, label2], fontsize=8)
    ax_g.set_title("Per-Check Status Grid", fontweight="bold", pad=6)
    for pi, res in enumerate([results1, results2]):
        lk = {r["check_id"]: r["status"] for r in res}
        for ci, cid in enumerate(all_cids):
            st = lk.get(cid, "")
            ax_g.text(ci, pi, st[:1], ha="center", va="center",
                      fontsize=5.5,
                      color="white" if st == "FAIL" else "#2c3e50",
                      fontweight="bold")
    _panel_label(ax_g, "B")

    # C — FAIL/WARN detail table
    ax_t.axis("off")
    rows_t = []
    for lbl, res in [(label1, results1), (label2, results2)]:
        for r in res:
            if r["status"] in ("FAIL", "WARN"):
                desc = (r["description"][:32] + "…"
                        if len(r["description"]) > 32 else r["description"])
                rows_t.append([r["check_id"], r["status"], lbl[:10], desc])
    if rows_t:
        tbl = ax_t.table(
            cellText=rows_t,
            colLabels=["Check", "Status", "Species", "Description"],
            loc="center", cellLoc="left", bbox=[0, 0, 1, 1],
        )
        tbl.auto_set_font_size(False); tbl.set_fontsize(6.5)
        tbl.auto_set_column_width([0, 1, 2, 3])
        for (r, c), cell in tbl.get_celld().items():
            cell.set_edgecolor("#dddddd"); cell.set_linewidth(0.4)
            if r == 0:
                cell.set_facecolor("#2c3e50")
                cell.set_text_props(color="white", fontweight="bold")
            elif rows_t[r-1][1] == "FAIL": cell.set_facecolor("#ffeaea")
            elif rows_t[r-1][1] == "WARN": cell.set_facecolor("#fff8e1")
            else: cell.set_facecolor("white")
    else:
        ax_t.text(0.5, 0.5, "All checks PASS",
                  ha="center", va="center",
                  color=PAL["pass"], fontsize=10, fontweight="bold")
    ax_t.set_title("FAIL / WARN Detail", fontweight="bold", pad=6)
    _panel_label(ax_t, "C")

    fig.suptitle("Cross-Species Ontology Validation Report",
                 fontsize=10, fontweight="bold", y=0.97)
    fig.savefig(out_path); plt.close(fig)
    print(f"  ✓  Fig 1 saved → {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Fig 2 — Species comparison overview (6-panel)
# ─────────────────────────────────────────────────────────────────────────────

def plot_species_comparison(
    df1: "pd.DataFrame",
    df2: "pd.DataFrame",
    comp_df: "pd.DataFrame",
    transfer: dict,
    label1: str = "human",
    label2: str = "otherSpecies",
    out_path: Path = Path("species_comparison_plot.png"),
) -> None:
    """
    Six-panel comparison overview:
      A — Dataset sizes (occurrences, proteins, domain types)
      B — Function-tag composition (stacked bars)
      C — CV-backed rates per column
      D — Topology context distribution
      E — Cross-species transfer scores
      F — Summary statistics table
    """
    from collections import Counter

    c1, c2 = PAL["blue"], PAL["orange"]

    def _cv(name, lbl):
        r = comp_df[comp_df["metric"] == name]
        return float(r[lbl].values[0]) if not r.empty else 0.0

    def _tag_counts(df):
        ctr: Counter = Counter()
        if "function_tags" not in df.columns: return ctr
        for val in df["function_tags"].dropna():
            for tok in str(val).split(";"):
                tok = tok.strip()
                # use label part after space if present
                lbl = tok.split(" ", 1)[1] if " " in tok else tok
                if lbl: ctr[lbl.lower()] += 1
        return ctr

    fig = plt.figure(figsize=(11, 11))
    gs  = fig.add_gridspec(3, 2, hspace=0.55, wspace=0.40,
                           left=0.09, right=0.97, top=0.93, bottom=0.04)
    ax_A = fig.add_subplot(gs[0, 0])
    ax_B = fig.add_subplot(gs[0, 1])
    ax_C = fig.add_subplot(gs[1, 0])
    ax_D = fig.add_subplot(gs[1, 1])
    ax_E = fig.add_subplot(gs[2, 0])
    ax_F = fig.add_subplot(gs[2, 1])

    # A — Dataset sizes
    labels_A = ["Occurrences", "Unique\nproteins", "Domain\ntypes"]
    names_A  = ["Total occurrence rows", "Unique proteins", "Unique domain types"]
    x = np.arange(3); w = 0.33
    v1 = [_cv(n, label1) for n in names_A]
    v2 = [_cv(n, label2) for n in names_A]
    ax_A.bar(x-w/2, v1, w, color=c1, edgecolor="white", lw=0.5, label=label1)
    ax_A.bar(x+w/2, v2, w, color=c2, edgecolor="white", lw=0.5, label=label2)
    m_ = max(v1 + v2)
    for xi, (a, b) in enumerate(zip(v1, v2)):
        ax_A.text(xi-w/2, a+m_*0.02, f"{int(a)}", ha="center", va="bottom", fontsize=6.5)
        ax_A.text(xi+w/2, b+m_*0.02, f"{int(b)}", ha="center", va="bottom", fontsize=6.5)
    ax_A.set_xticks(x); ax_A.set_xticklabels(labels_A)
    ax_A.set_ylabel("Count"); ax_A.set_title("Dataset Size", fontweight="bold", pad=6)
    ax_A.legend(fontsize=6.5, handlelength=1.2, loc="upper right")
    _panel_label(ax_A, "A")

    # B — Function-tag composition stacked bars
    ct1 = _tag_counts(df1); ct2 = _tag_counts(df2)
    all_tags = sorted((ct1 | ct2).keys(), key=lambda t: -(ct1.get(t,0)+ct2.get(t,0)))
    top_tags = all_tags[:14]
    tpal = [PAL["blue"],PAL["orange"],PAL["green"],PAL["red"],PAL["purple"],
            PAL["yellow"],PAL["cyan"],"#997700","#004488","#882255",
            "#BBBBBB","#44AA99","#EE7733","#009988"][:len(top_tags)]
    n1 = max(sum(ct1.values()), 1); n2 = max(sum(ct2.values()), 1)
    for yi, (ctr, total) in enumerate([(ct1, n1), (ct2, n2)]):
        left = 0.0
        for tag, col_ in zip(top_tags, tpal):
            frac = ctr.get(tag, 0) / total
            if frac == 0: continue
            ax_B.barh(1-yi, frac*100, height=0.38, left=left,
                      color=col_, edgecolor="white", lw=0.4)
            if frac*100 > 4:
                ax_B.text(left+frac*100/2, 1-yi, f"{frac*100:.0f}%",
                          ha="center", va="center",
                          fontsize=5.5, color="white", fontweight="bold")
            left += frac*100
    ax_B.set_yticks([0, 1]); ax_B.set_yticklabels([label2, label1], fontsize=8)
    ax_B.set_xlabel("% of domain occurrences")
    ax_B.set_title("Function-Tag Composition\n(top 14 tags)", fontweight="bold", pad=6)
    ax_B.set_xlim(0, 105); ax_B.spines["left"].set_visible(False)
    ax_B.tick_params(axis="y", length=0)
    ax_B.legend(handles=[mpatches.Patch(color=c, label=t[:20])
                          for t, c in zip(top_tags, tpal)],
                fontsize=5, ncol=2, loc="lower right",
                handlelength=1.0, handletextpad=0.4, columnspacing=0.8,
                bbox_to_anchor=(1.0, -0.22))
    _panel_label(ax_B, "B")

    # C — CV-backed rates (dot-lollipop)
    rate_rows = [
        ("CV-backed rate: function_tags (%)",     "Function tags"),
        ("CV-backed rate: topology_context (%)",  "Topology"),
        ("CV-backed rate: positional_category (%)", "Positional"),
        ("CV-backed rate: binding_partner (%)",   "Binding partner"),
    ]
    yC = np.arange(len(rate_rows))
    r1C = [_cv(m, label1) for m, _ in rate_rows]
    r2C = [_cv(m, label2) for m, _ in rate_rows]
    for i, (r1, r2) in enumerate(zip(r1C, r2C)):
        ax_C.plot([r1], [i-0.15], "o", color=c1, ms=7,
                  markeredgecolor="white", markeredgewidth=0.6, zorder=3)
        ax_C.plot([r2], [i+0.15], "D", color=c2, ms=7,
                  markeredgecolor="white", markeredgewidth=0.6, zorder=3)
        ax_C.plot([r1, r2], [i-0.15, i+0.15], color="#cccccc", lw=0.8, zorder=1)
        ax_C.text(r1+1, i-0.15, f"{r1:.0f}%", va="center", fontsize=6, color=c1)
        ax_C.text(r2+1, i+0.15, f"{r2:.0f}%", va="center", fontsize=6, color=c2)
    ax_C.set_yticks(yC); ax_C.set_yticklabels([s for _, s in rate_rows])
    ax_C.set_xlim(-5, 115); ax_C.axvline(100, color="#cccccc", lw=0.7, ls="--")
    ax_C.set_xlabel("CV-backed rate (%)")
    ax_C.set_title("Ontology Coverage per Column", fontweight="bold", pad=6)
    ax_C.legend(handles=[mpatches.Patch(color=c1, label=label1),
                          mpatches.Patch(color=c2, label=label2)],
                fontsize=6.5, loc="lower right", handlelength=1.2)
    _panel_label(ax_C, "C")

    # D — Topology context distribution comparison
    def _topo_counts(df):
        c: Counter = Counter()
        if "topology_context" not in df.columns: return c
        for val in df["topology_context"].dropna():
            for tok in str(val).split(";"):
                tok = tok.strip().split(" ", 1)
                label = tok[1] if len(tok) > 1 else tok[0]
                if label: c[label.lower()] += 1
        return c
    tc1 = _topo_counts(df1); tc2 = _topo_counts(df2)
    all_topo = sorted((tc1 | tc2).keys(),
                       key=lambda t: -(tc1.get(t,0)+tc2.get(t,0)))[:8]
    xD = np.arange(len(all_topo)); wD = 0.35
    v1D = [tc1.get(t,0) for t in all_topo]
    v2D = [tc2.get(t,0) for t in all_topo]
    ax_D.bar(xD-wD/2, v1D, wD, color=c1, edgecolor="white", lw=0.4)
    ax_D.bar(xD+wD/2, v2D, wD, color=c2, edgecolor="white", lw=0.4)
    ax_D.set_xticks(xD)
    ax_D.set_xticklabels([t[:12] for t in all_topo], rotation=40, ha="right")
    ax_D.set_ylabel("Count")
    ax_D.set_title("Topology Context Distribution", fontweight="bold", pad=6)
    ax_D.legend(handles=[mpatches.Patch(color=c1, label=label1),
                          mpatches.Patch(color=c2, label=label2)],
                fontsize=6.5, handlelength=1.2, loc="upper right")
    _panel_label(ax_D, "D")

    # E — Transfer scores (horizontal bars)
    tm_names = ["Schema reuse", "Vocabulary reuse",
                 "Rule portability", "Mapping completeness"]
    tm_keys  = ["schema_reuse_rate", "vocab_reuse_rate",
                 "rule_portability_rate", "mapping_completeness"]
    tm_vals  = [transfer.get(k, 0.0) for k in tm_keys]
    yE       = np.arange(4); cmap_E = plt.cm.RdYlGn
    bE = ax_E.barh(yE, tm_vals, height=0.55,
                   color=[cmap_E(v) for v in tm_vals],
                   edgecolor="white", lw=0.5)
    for bar, v in zip(bE, tm_vals):
        fc = "white" if v < 0.35 or v > 0.75 else "#2c3e50"
        ax_E.text(min(v-0.01, 1.08),
                  bar.get_y()+bar.get_height()/2,
                  f"{v:.1%}", ha="right", va="center",
                  fontsize=7, fontweight="bold", color=fc, clip_on=True)
    ax_E.axvline(0.5, color="#aaaaaa", lw=0.7, ls="--", alpha=0.7)
    ax_E.set_xlim(0, 1.12); ax_E.set_yticks(yE)
    ax_E.set_yticklabels(tm_names)
    ax_E.set_xlabel("Score (0 → 1)")
    ax_E.set_title(f"Transfer Performance\n{label1} → {label2}",
                   fontweight="bold", pad=6)
    sm_E = plt.cm.ScalarMappable(cmap=cmap_E, norm=plt.Normalize(0, 1))
    sm_E.set_array([])
    fig.colorbar(sm_E, ax=ax_E, shrink=0.65, pad=0.02, aspect=18).ax.tick_params(labelsize=5.5)
    _panel_label(ax_E, "E")

    # F — Summary statistics table
    ax_F.axis("off")
    rows_F = [
        ["Metric",               label1[:12],           label2[:12]],
        ["Occurrence rows",      str(int(_cv("Total occurrence rows", label1))),
                                 str(int(_cv("Total occurrence rows", label2)))],
        ["Unique proteins",      str(int(_cv("Unique proteins", label1))),
                                 str(int(_cv("Unique proteins", label2)))],
        ["Domain types",         str(int(_cv("Unique domain types", label1))),
                                 str(int(_cv("Unique domain types", label2)))],
        ["Schema reuse",         f"{transfer.get('schema_reuse_rate',0):.1%}", "←"],
        ["Vocab reuse",          f"{transfer.get('vocab_reuse_rate',0):.1%}", "←"],
        ["Shared CV tokens",     str(transfer.get("n_shared_cv_tokens","—")), "↑↑"],
        ["New terms needed",     "—", str(transfer.get("n_new_terms_required","—"))],
        ["Conserved ont. terms", str(transfer.get("n_terms_conserved","—")), "↑↑"],
        ["Mapping completeness", f"{transfer.get('mapping_completeness',0):.1%}", "same"],
    ]
    tbl = ax_F.table(cellText=rows_F[1:], colLabels=rows_F[0],
                     loc="center", cellLoc="center", bbox=[0, 0, 1, 1])
    tbl.auto_set_font_size(False); tbl.set_fontsize(6.5)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("#dddddd"); cell.set_linewidth(0.4)
        if r == 0:
            cell.set_facecolor("#2c3e50")
            cell.set_text_props(color="white", fontweight="bold")
        elif r % 2 == 0: cell.set_facecolor("#f7f7f7")
        else: cell.set_facecolor("white")
    ax_F.set_title("Summary Statistics", fontweight="bold", pad=6)
    _panel_label(ax_F, "F", x=-0.04)

    fig.suptitle(f"Cross-Species Comparison: {label1} vs {label2}",
                 fontsize=10, fontweight="bold")
    fig.savefig(out_path); plt.close(fig)
    print(f"  ✓  Fig 2 saved → {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Fig 3 — Jaccard overlap (lollipop + set-size bars)
# ─────────────────────────────────────────────────────────────────────────────

def plot_jaccard_heatmap(
    overlap_df: "pd.DataFrame",
    out_path: Path,
    label1: str = "human",
    label2: str = "otherSpecies",
) -> None:
    """
    Two-panel biological overlap figure.
      A — Lollipop chart of Jaccard similarity per category
      B — Stacked set-size bars (shared / unique per species)
    """
    cats    = overlap_df["category"].tolist()
    jacs    = overlap_df["jaccard"].astype(float).tolist()
    shared  = overlap_df["shared"].astype(int).tolist()
    u1c, u2c = f"unique_to_{label1}", f"unique_to_{label2}"
    s1c, s2c = f"{label1}_size", f"{label2}_size"
    uniq1 = (overlap_df[u1c].astype(int).tolist() if u1c in overlap_df.columns
             else (overlap_df[s1c].astype(int)-overlap_df["shared"].astype(int)).tolist())
    uniq2 = (overlap_df[u2c].astype(int).tolist() if u2c in overlap_df.columns
             else (overlap_df[s2c].astype(int)-overlap_df["shared"].astype(int)).tolist())

    n  = len(cats); y = np.arange(n)
    cmap_j = plt.cm.RdYlGn

    fig, (ax_A, ax_B) = plt.subplots(
        1, 2, figsize=(9, max(3.5, n*0.52+1.2)),
        gridspec_kw={"width_ratios": [2.2, 3]},
    )

    # A — lollipop
    for i, (j, cat) in enumerate(zip(jacs, cats)):
        ax_A.plot([0, j], [i, i], color="#d5d5d5", lw=1.2, zorder=1)
        ax_A.scatter(j, i, s=65, color=cmap_j(j), zorder=3,
                     linewidths=0.6, edgecolors="white")
        ax_A.text(j+0.025, i, f"{j:.3f}", va="center",
                  fontsize=6.5, color="#2c3e50")
    ax_A.axvline(0.5, color="#aaaaaa", lw=0.8, ls="--", alpha=0.7)
    ax_A.set_xlim(-0.04, 1.22); ax_A.set_yticks(y)
    ax_A.set_yticklabels(cats, fontsize=7.5)
    ax_A.set_xlabel("Jaccard similarity index")
    ax_A.set_title(f"Cross-Species Category Overlap\n{label1}  vs  {label2}",
                   fontweight="bold", pad=6)
    sm = plt.cm.ScalarMappable(cmap=cmap_j, norm=plt.Normalize(0, 1))
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax_A, shrink=0.5, pad=0.03, aspect=16)
    cb.set_label("Jaccard", fontsize=6.5); cb.ax.tick_params(labelsize=6)
    _panel_label(ax_A, "A")

    # B — stacked set-size bars
    bh = 0.22
    ax_B.barh(y+bh, shared, height=bh, color=PAL["shared"], label="Shared")
    ax_B.barh(y,    uniq1,  height=bh, color=PAL["ref_only"], label=f"Unique {label1}")
    ax_B.barh(y-bh, uniq2,  height=bh, color=PAL["tgt_only"], label=f"Unique {label2}")
    for i, (sv, u1v, u2v) in enumerate(zip(shared, uniq1, uniq2)):
        for val, dy in [(sv, bh), (u1v, 0), (u2v, -bh)]:
            if val > 0:
                ax_B.text(val+0.3, i+dy, str(val), va="center",
                          fontsize=6, color="#2c3e50")
    ax_B.set_yticks(y); ax_B.set_yticklabels([])
    ax_B.set_xlabel("Item count")
    ax_B.set_title("Set Breakdown\n(shared / unique per species)",
                   fontweight="bold", pad=6)
    ax_B.legend(fontsize=6.5, loc="lower right", handlelength=1.2)
    _panel_label(ax_B, "B")

    fig.tight_layout()
    fig.savefig(out_path); plt.close(fig)
    print(f"  ✓  Fig 3 saved → {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Fig 4 — Transfer metrics (radar + term conservation bars)
# ─────────────────────────────────────────────────────────────────────────────

def plot_transfer_metrics(
    transfer: dict,
    out_path: Path,
    label1: str = "human",
    label2: str = "otherSpecies",
) -> None:
    """
    Two-panel transfer metrics figure.
      A — Spider/radar chart of 4 key portability metrics
      B — Ontology term conservation breakdown bars
    """
    fig = plt.figure(figsize=(9.5, 4.5))
    gs  = fig.add_gridspec(1, 2, wspace=0.45,
                           left=0.04, right=0.97, top=0.88, bottom=0.10)
    ax_B = fig.add_subplot(gs[1])    # bars first (easier)
    ax_A = fig.add_subplot(gs[0], polar=True)

    # A — radar
    mls   = ["Schema\nreuse", "Vocabulary\nreuse",
              "Rule\nportability", "Mapping\ncompleteness"]
    vals  = [transfer.get(k, 0.0) for k in
             ["schema_reuse_rate", "vocab_reuse_rate",
              "rule_portability_rate", "mapping_completeness"]]
    N     = len(mls)
    angs  = np.linspace(0, 2*np.pi, N, endpoint=False).tolist()
    angs += angs[:1]; vs = vals + vals[:1]

    for lv, ls in [(0.25,":"),(0.5,"--"),(0.75,":"),(1.0,"-")]:
        ax_A.plot(angs, [lv]*(N+1), color="#cccccc", lw=0.6, ls=ls, zorder=1)
        ax_A.text(np.pi/2, lv+0.04, f"{lv:.0%}",
                  ha="center", va="bottom", fontsize=5.5, color="#999999")
    ax_A.fill(angs, vs, color=PAL["blue"], alpha=0.18, zorder=2)
    ax_A.plot(angs, vs, color=PAL["blue"], lw=2.0,
              marker="o", ms=6, markerfacecolor=PAL["blue"],
              markeredgecolor="white", markeredgewidth=1.0, zorder=3)
    for angle, val in zip(angs[:-1], vals):
        ax_A.text(angle, val+0.12, f"{val:.1%}",
                  ha="center", va="center",
                  fontsize=7.5, fontweight="bold", color=PAL["blue"])
    ax_A.set_xticks(angs[:-1]); ax_A.set_xticklabels(mls, fontsize=8)
    ax_A.set_ylim(0, 1.18); ax_A.set_yticks([])
    ax_A.spines["polar"].set_color("#cccccc")
    ax_A.set_title(f"Ontology Transfer Metrics\n{label1} → {label2}",
                   fontsize=9, fontweight="bold", pad=18)
    _panel_label(ax_A, "A", x=-0.08, y=1.08)

    # B — term conservation breakdown
    cats_B   = ["Conserved\n(both)", f"{label1[:8]}\nonly",
                 f"{label2[:8]}\nonly", "New terms\nneeded", "Registry\nunused"]
    keys_B   = ["n_terms_conserved", "n_terms_ref_only", "n_terms_target_only",
                 "n_new_terms_required", "n_registry_terms_unused"]
    vals_B   = [transfer.get(k, 0) for k in keys_B]
    cols_B   = [PAL["pass"], PAL["ref_only"], PAL["tgt_only"], PAL["warn"], PAL["grey"]]
    yB       = np.arange(len(cats_B))
    barsB    = ax_B.barh(yB, vals_B, height=0.55, color=cols_B,
                          edgecolor="white", lw=0.5)
    mB = max(vals_B + [1])
    for bar, v in zip(barsB, vals_B):
        if v > 0:
            ax_B.text(v + mB*0.02, bar.get_y()+bar.get_height()/2,
                      str(v), va="center", fontsize=7, color="#2c3e50")
    ax_B.set_yticks(yB); ax_B.set_yticklabels(cats_B)
    ax_B.set_xlabel("Number of CV terms")
    ax_B.set_title(f"Ontology Term Conservation\n{label1} → {label2}",
                   fontweight="bold", pad=6)
    _panel_label(ax_B, "B")

    fig.suptitle(f"Ontology Transfer Analysis: {label1} → {label2}",
                 fontsize=9.5, fontweight="bold", y=1.01)
    fig.savefig(out_path); plt.close(fig)
    print(f"  ✓  Fig 4 saved → {out_path}")


# Keep legacy name so existing callers still work
def plot_term_conservation(transfer, out_path, label_ref="human",
                            label_tgt="otherSpecies"):
    plot_transfer_metrics(transfer, out_path, label_ref, label_tgt)


# ─────────────────────────────────────────────────────────────────────────────
# Fig 5 — Domain conservation profile
# ─────────────────────────────────────────────────────────────────────────────

def plot_domain_conservation(
    df1: "pd.DataFrame",
    df2: "pd.DataFrame",
    label1: str = "human",
    label2: str = "otherSpecies",
    out_path: Path = Path("domain_conservation_plot.png"),
) -> None:
    """
    Two-panel domain-level conservation figure.
      A — Top shared and species-specific domain types by occurrence count
      B — Occurrence-count scatter: shared vs species-specific domains
    """
    if "domain_type_id" not in df1.columns or "domain_type_id" not in df2.columns:
        print("  ⚠  domain_type_id absent — skipping domain conservation plot")
        return

    s1  = set(df1["domain_type_id"].dropna().unique())
    s2  = set(df2["domain_type_id"].dropna().unique())
    cnt1 = df1["domain_type_id"].value_counts()
    cnt2 = df2["domain_type_id"].value_counts()

    def _short(d):
        return d.split(":")[-1] if ":" in d else d

    sh_d   = sorted(s1 & s2,   key=lambda x: -(cnt1.get(x,0)+cnt2.get(x,0)))[:8]
    only1  = sorted(s1 - s2,   key=lambda x: -cnt1.get(x,0))[:6]
    only2  = sorted(s2 - s1,   key=lambda x: -cnt2.get(x,0))[:6]
    doms   = [_short(d) for d in only1] + [_short(d) for d in sh_d] + [_short(d) for d in only2]
    v1s    = [cnt1.get(d,0) for d in only1] + [cnt1.get(d,0) for d in sh_d] + [0]*len(only2)
    v2s    = [0]*len(only1) + [cnt2.get(d,0) for d in sh_d] + [cnt2.get(d,0) for d in only2]
    cc     = ([PAL["ref_only"]]*len(only1) + [PAL["shared"]]*len(sh_d)
              + [PAL["tgt_only"]]*len(only2))

    fig, (ax_A, ax_B) = plt.subplots(1, 2, figsize=(11, 5.5),
                                      gridspec_kw={"width_ratios": [3, 2.5]},
                                      constrained_layout=True)

    # A — top domain types
    yA = np.arange(len(doms)); wA = 0.33
    ax_A.barh(yA+wA/2, v1s, height=wA, color=cc, edgecolor="white", lw=0.4, alpha=0.85)
    ax_A.barh(yA-wA/2, v2s, height=wA,
              color=[PAL["tgt_only"] if c==PAL["ref_only"] else
                     PAL["ref_only"] if c==PAL["tgt_only"] else c for c in cc],
              edgecolor="white", lw=0.4, alpha=0.65)
    ax_A.set_yticks(yA); ax_A.set_yticklabels(doms, fontsize=6.5)
    ax_A.set_xlabel("Domain occurrence count")
    ax_A.set_title(
        f"Top Domain Types by Category\n"
        f"Green=shared  Blue={label1}-only  Orange={label2}-only",
        fontweight="bold", pad=6, fontsize=7.5)
    mA = max(v1s+v2s+[1])
    for lbt, mid_y, col_ in [
        (f"{label1}-only", len(only1)/2-0.5,                PAL["ref_only"]),
        ("Shared",         len(only1)+len(sh_d)/2-0.5,      PAL["shared"]),
        (f"{label2}-only", len(only1)+len(sh_d)+len(only2)/2-0.5, PAL["tgt_only"]),
    ]:
        ax_A.text(-mA*0.18, mid_y, lbt, ha="right", va="center",
                  fontsize=6, color=col_, fontweight="bold")
    ax_A.legend(handles=[mpatches.Patch(color=PAL["blue"],   label=label1),
                          mpatches.Patch(color=PAL["orange"], label=label2)],
                fontsize=6.5, loc="lower right", handlelength=1.2)
    _panel_label(ax_A, "A")

    # B — scatter: shared vs species-specific
    all_dom  = list(s1 | s2)
    x_vals   = np.array([cnt1.get(d, 0) for d in all_dom], dtype=float)
    y_vals   = np.array([cnt2.get(d, 0) for d in all_dom], dtype=float)
    is_shared = np.array([d in s1 and d in s2 for d in all_dom])

    ax_B.scatter(x_vals[is_shared],  y_vals[is_shared],
                 s=22, color=PAL["shared"],  alpha=0.75, lw=0.3,
                 edgecolors="#2c3e50", label="Shared", zorder=3)
    ax_B.scatter(x_vals[~is_shared & (x_vals>0)], y_vals[~is_shared & (x_vals>0)],
                 s=22, color=PAL["ref_only"], alpha=0.75, lw=0.3,
                 edgecolors="#2c3e50", label=f"{label1}-only", zorder=3)
    ax_B.scatter(x_vals[~is_shared & (y_vals>0)], y_vals[~is_shared & (y_vals>0)],
                 s=22, color=PAL["tgt_only"], alpha=0.75, lw=0.3,
                 edgecolors="#2c3e50", label=f"{label2}-only", zorder=3)

    lim = max(x_vals.max(), y_vals.max()) * 1.08
    ax_B.plot([0, lim], [0, lim], color="#aaaaaa", lw=0.8, ls="--",
              alpha=0.7, zorder=1, label="1:1 line")
    ax_B.set_xlim(-0.5, lim); ax_B.set_ylim(-0.5, lim)
    ax_B.set_xlabel(f"{label1} occurrence count")
    ax_B.set_ylabel(f"{label2} occurrence count")
    ax_B.set_title("Domain Occurrence Counts\n(shared vs species-specific)",
                   fontweight="bold", pad=6)
    ax_B.legend(fontsize=6.5, handlelength=1.2)
    _panel_label(ax_B, "B")

    fig.suptitle(f"Domain Conservation: {label1} vs {label2}",
                 fontsize=9.5, fontweight="bold")
    fig.savefig(out_path); plt.close(fig)
    print(f"  ✓  Fig 5 saved → {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Fig 6 — Function-tag cross-species profile
# ─────────────────────────────────────────────────────────────────────────────

def plot_function_profile(
    df1: "pd.DataFrame",
    df2: "pd.DataFrame",
    label1: str = "human",
    label2: str = "otherSpecies",
    out_path: Path = Path("function_profile_plot.png"),
) -> None:
    """
    Two-panel function-tag cross-species profile.
      A — Stacked horizontal bars: function-tag fraction per species
      B — Per-tag conservation lollipop (Jaccard per DOF category)
    """
    from collections import Counter

    def _tag_ctr(df):
        c: Counter = Counter()
        if "function_tags" not in df.columns: return c
        for val in df["function_tags"].dropna():
            for tok in str(val).split(";"):
                tok = tok.strip()
                lbl = tok.split(" ", 1)[1] if " " in tok else tok
                if lbl: c[lbl.lower()] += 1
        return c

    def _tag_set(df):
        s: dict[str, set] = {}
        if "function_tags" not in df.columns: return s
        if "protein_accession" not in df.columns: return s
        for prot, grp in df.groupby("protein_accession"):
            tags: set = set()
            for val in grp["function_tags"].dropna():
                for tok in str(val).split(";"):
                    tok = tok.strip()
                    lbl = tok.split(" ", 1)[1] if " " in tok else tok
                    if lbl: tags.add(lbl.lower())
            s[prot] = tags
        return s

    ct1 = _tag_ctr(df1); ct2 = _tag_ctr(df2)
    all_tags = sorted((ct1|ct2).keys(), key=lambda t: -(ct1.get(t,0)+ct2.get(t,0)))
    top_tags = all_tags[:14]
    tpal     = [PAL["blue"],PAL["orange"],PAL["green"],PAL["red"],PAL["purple"],
                PAL["yellow"],PAL["cyan"],"#997700","#004488","#882255",
                "#BBBBBB","#44AA99","#EE7733","#009988"][:len(top_tags)]
    n1 = max(sum(ct1.values()), 1); n2 = max(sum(ct2.values()), 1)

    fig = plt.figure(figsize=(12, 5.5))
    gs  = fig.add_gridspec(1, 2, wspace=0.42,
                           left=0.07, right=0.97, top=0.88, bottom=0.15,
                           width_ratios=[2, 3])
    ax_A = fig.add_subplot(gs[0]); ax_B = fig.add_subplot(gs[1])

    # A — stacked composition bars
    for yi, (ctr, total) in enumerate([(ct1, n1), (ct2, n2)]):
        left = 0.0
        for tag, col_ in zip(top_tags, tpal):
            frac = ctr.get(tag, 0) / total
            if frac == 0: continue
            ax_A.barh(1-yi, frac*100, height=0.38, left=left,
                      color=col_, edgecolor="white", lw=0.4)
            if frac*100 > 4:
                ax_A.text(left+frac*100/2, 1-yi, f"{frac*100:.0f}%",
                          ha="center", va="center",
                          fontsize=5.5, color="white", fontweight="bold")
            left += frac*100
    ax_A.set_yticks([0, 1]); ax_A.set_yticklabels([label2, label1], fontsize=8)
    ax_A.set_xlabel("% of domain occurrences")
    ax_A.set_title("Function-Tag Composition\n(top 14 tags)",
                   fontweight="bold", pad=6)
    ax_A.set_xlim(0, 105); ax_A.spines["left"].set_visible(False)
    ax_A.tick_params(axis="y", length=0)
    ax_A.legend(handles=[mpatches.Patch(color=c, label=t[:22])
                          for t, c in zip(top_tags, tpal)],
                fontsize=5, ncol=2, loc="lower right",
                handlelength=1.0, handletextpad=0.4, columnspacing=0.8,
                bbox_to_anchor=(1.0, -0.25))
    _panel_label(ax_A, "A")

    # B — per-tag Jaccard lollipop
    def _jaccard_sets(a: set, b: set) -> float:
        u = a | b; return len(a&b)/max(len(u),1)

    # per-tag presence as sets of proteins that carry each tag
    pts1: dict[str, set] = {}
    pts2: dict[str, set] = {}
    for tag in top_tags:
        if "protein_accession" in df1.columns:
            pts1[tag] = set(
                df1[df1["function_tags"].fillna("").str.contains(
                    tag.split()[0] if " " in tag else tag, regex=False)
                    ]["protein_accession"])
        else: pts1[tag] = set()
        if "protein_accession" in df2.columns:
            pts2[tag] = set(
                df2[df2["function_tags"].fillna("").str.contains(
                    tag.split()[0] if " " in tag else tag, regex=False)
                    ]["protein_accession"])
        else: pts2[tag] = set()

    jacs_B = [_jaccard_sets(pts1.get(t,set()), pts2.get(t,set())) for t in top_tags]
    # sort by Jaccard descending for readability
    order  = sorted(range(len(top_tags)), key=lambda i: -jacs_B[i])
    s_tags = [top_tags[i] for i in order]
    s_jacs = [jacs_B[i]   for i in order]
    s_cols = [tpal[i]      for i in order]

    yB = np.arange(len(s_tags)); cmap_j = plt.cm.RdYlGn
    for i, (j, col_) in enumerate(zip(s_jacs, s_cols)):
        ax_B.plot([0, j], [i, i], color="#d5d5d5", lw=1.0, zorder=1)
        ax_B.scatter(j, i, s=55, color=cmap_j(j), zorder=3,
                     linewidths=0.6, edgecolors=col_)
        ax_B.text(j+0.025, i, f"{j:.3f}", va="center",
                  fontsize=6.5, color="#2c3e50")
    ax_B.axvline(0.5, color="#aaaaaa", lw=0.8, ls="--", alpha=0.7)
    ax_B.axvline(0.7, color=PAL["green"], lw=0.7, ls=":", alpha=0.6,
                 label="High (≥0.7)")
    ax_B.set_yticks(yB)
    ax_B.set_yticklabels([t[:25] for t in s_tags], fontsize=7)
    ax_B.set_xlim(-0.04, 1.22)
    ax_B.set_xlabel("Jaccard similarity (protein-set overlap per tag)")
    ax_B.set_title("Per-Tag Protein-Set Conservation\n(Jaccard per function category)",
                   fontweight="bold", pad=6)
    sm_B = plt.cm.ScalarMappable(cmap=cmap_j, norm=plt.Normalize(0, 1))
    sm_B.set_array([])
    cb_B = fig.colorbar(sm_B, ax=ax_B, shrink=0.55, pad=0.02, aspect=16)
    cb_B.set_label("Jaccard", fontsize=6.5); cb_B.ax.tick_params(labelsize=6)
    _panel_label(ax_B, "B")

    fig.suptitle(f"Functional Profile Conservation: {label1} vs {label2}",
                 fontsize=9.5, fontweight="bold", y=0.97)
    fig.savefig(out_path); plt.close(fig)
    print(f"  ✓  Fig 6 saved → {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Fig 7 — Ortholog conservation (only when ortholog file provided)
# ─────────────────────────────────────────────────────────────────────────────

def plot_ortholog_conservation(
    ortho_df: "pd.DataFrame",
    out_path: Path,
    label1: str = "human",
    label2: str = "otherSpecies",
) -> None:
    """
    Two-panel ortholog conservation figure.
      A — Mean Jaccard per annotation category (horizontal lollipop)
      B — Distribution of overall conservation scores across ortholog pairs
    """
    complete = ortho_df[ortho_df["pair_complete"]] if "pair_complete" in ortho_df.columns else ortho_df
    if len(complete) == 0:
        print("  ⚠  No complete ortholog pairs — skipping ortholog conservation plot")
        return

    cats  = list(CONSERVATION_GROUPS.keys())
    means = []
    for cat in cats:
        col = f"jaccard_{cat}"
        means.append(round(complete[col].mean(), 3)
                     if col in complete.columns else 0.0)

    fig, (ax_A, ax_B) = plt.subplots(1, 2, figsize=(9.5, max(3.5, len(cats)*0.55+1.2)),
                                      constrained_layout=True)
    cmap_j = plt.cm.RdYlGn; y = np.arange(len(cats))

    # A — lollipop
    for i, (cat, v) in enumerate(zip(cats, means)):
        ax_A.plot([0, v], [i, i], color="#d5d5d5", lw=1.2, zorder=1)
        ax_A.scatter(v, i, s=65, color=cmap_j(v), zorder=3,
                     linewidths=0.6, edgecolors="white")
        ax_A.text(v+0.025, i, f"{v:.3f}", va="center",
                  fontsize=6.5, color="#2c3e50")
    ax_A.axvline(0.7, color=PAL["green"],  ls="--", lw=0.8, alpha=0.7, label="High (≥0.7)")
    ax_A.axvline(0.4, color=PAL["orange"], ls="--", lw=0.8, alpha=0.7, label="Medium (≥0.4)")
    ax_A.set_xlim(-0.04, 1.22); ax_A.set_yticks(y)
    ax_A.set_yticklabels(cats, fontsize=8)
    ax_A.set_xlabel("Mean Jaccard conservation score")
    ax_A.set_title(f"Ortholog Conservation per Category\n"
                   f"{label1} → {label2}  ({len(complete):,} pairs)",
                   fontweight="bold", pad=6)
    ax_A.legend(fontsize=7, loc="lower right", handlelength=1.2)
    sm = plt.cm.ScalarMappable(cmap=cmap_j, norm=plt.Normalize(0, 1)); sm.set_array([])
    fig.colorbar(sm, ax=ax_A, shrink=0.55, pad=0.02, aspect=16,
                 label="Conservation").ax.tick_params(labelsize=6)
    _panel_label(ax_A, "A")

    # B — distribution of mean conservation score across pairs
    if "mean_conservation" in complete.columns:
        vals = complete["mean_conservation"].dropna().astype(float).values
        ax_B.hist(vals, bins=20, color=PAL["blue"], edgecolor="white", lw=0.5, alpha=0.85)
        mu = float(np.mean(vals)); med = float(np.median(vals))
        ax_B.axvline(mu,  color=PAL["orange"], lw=1.5, ls="--",
                     label=f"Mean={mu:.3f}")
        ax_B.axvline(med, color=PAL["green"],  lw=1.5, ls=":",
                     label=f"Median={med:.3f}")
        ax_B.axvline(0.7, color="#aaaaaa", lw=0.8, ls="-.",
                     label="High threshold (0.7)", alpha=0.6)
        ax_B.set_xlabel("Overall conservation score")
        ax_B.set_ylabel("Number of ortholog pairs")
        ax_B.set_title("Distribution of Pair Conservation Scores",
                        fontweight="bold", pad=6)
        ax_B.legend(fontsize=7, handlelength=1.2)
    else:
        ax_B.text(0.5, 0.5, "mean_conservation\ncolumn absent",
                  ha="center", va="center", transform=ax_B.transAxes,
                  fontsize=9, color="#888888")
        ax_B.set_title("Distribution of Pair Conservation Scores",
                        fontweight="bold", pad=6)
    _panel_label(ax_B, "B")

    fig.suptitle(f"Ortholog-Level Conservation: {label1} vs {label2}",
                 fontsize=9.5, fontweight="bold")
    fig.savefig(out_path); plt.close(fig)
    print(f"  ✓  Fig 7 saved → {out_path}")


def _check_mpl() -> bool:
    """Kept for backward compatibility — always True when matplotlib is installed."""
    try:
        import matplotlib  # noqa: F401
        import numpy       # noqa: F401
        return True
    except ImportError:
        print("  ⚠  matplotlib/numpy not installed — skipping plots")
        return False


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 12 — Main pipeline orchestrator
# ══════════════════════════════════════════════════════════════════════════════


def plot_term_conservation(
    transfer: dict,
    out_path: Path,
    label_ref: str = "human",
    label_tgt: str = "otherSpecies",
) -> None:
    """
    Grouped bar chart showing conserved, species-specific, and unused terms.
    """
    if not _check_mpl():
        return
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.use("Agg")

    categories = [
        "Conserved\n(both species)",
        f"Ref only\n({label_ref})",
        f"Target only\n({label_tgt})",
        "New terms\nrequired",
        "Registry\nunused",
    ]
    counts = [
        transfer.get("n_terms_conserved", 0),
        transfer.get("n_terms_ref_only", 0),
        transfer.get("n_terms_target_only", 0),
        transfer.get("n_new_terms_required", 0),
        transfer.get("n_registry_terms_unused", 0),
    ]
    colours = ["#2ecc71", "#3498db", "#e67e22", "#e74c3c", "#95a5a6"]

    fig, ax = plt.subplots(figsize=(10, 5))
    bars = ax.bar(categories, counts, color=colours,
                  edgecolor="white", linewidth=0.8, width=0.55)
    for bar, val in zip(bars, counts):
        if val > 0:
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + max(counts) * 0.015,
                    str(val), ha="center", va="bottom",
                    fontsize=10, fontweight="bold")
    ax.set_ylabel("Number of CV terms", fontsize=12)
    ax.set_title(
        f"Ontology Term Conservation: {label_ref} → {label_tgt}",
        fontsize=13, fontweight="bold"
    )
    ax.set_ylim(0, max(counts) * 1.2 if max(counts) > 0 else 10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  ✓  Plot saved → {out_path}")


def plot_jaccard_heatmap(
    overlap_df: pd.DataFrame,
    out_path: Path,
    label1: str = "human",
    label2: str = "otherSpecies",
) -> None:
    """
    Horizontal heatmap of Jaccard similarity values per biological category.
    """
    if not _check_mpl():
        return
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    matplotlib.use("Agg")

    cats     = overlap_df["category"].tolist()
    jaccards = overlap_df["jaccard"].tolist()
    data     = [[j] for j in jaccards]

    fig, ax = plt.subplots(figsize=(4.5, max(4, len(cats) * 0.58)))
    im = ax.imshow(data, aspect="auto",
                   cmap="YlGnBu", vmin=0, vmax=1)
    ax.set_xticks([0])
    ax.set_xticklabels([f"Jaccard\n{label1} vs {label2}"], fontsize=10)
    ax.set_yticks(range(len(cats)))
    ax.set_yticklabels(cats, fontsize=9)
    for i, j in enumerate(jaccards):
        ax.text(0, i, f"{j:.3f}", ha="center", va="center",
                fontsize=9,
                color="black" if j < 0.6 else "white",
                fontweight="bold")
    plt.colorbar(im, ax=ax, label="Jaccard similarity", shrink=0.7)
    ax.set_title("Cross-Species Category Overlap",
                 fontsize=12, fontweight="bold", pad=10)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  ✓  Plot saved → {out_path}")


def plot_validation_summary(
    results1: list[dict],
    results2: list[dict],
    out_path: Path,
    label1: str = "human",
    label2: str = "otherSpecies",
) -> None:
    """
    Grouped bar chart of PASS / FAIL / WARN check counts per species.
    """
    if not _check_mpl():
        return
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    matplotlib.use("Agg")

    from collections import Counter
    c1 = Counter(r["status"] for r in results1)
    c2 = Counter(r["status"] for r in results2)
    statuses   = ["PASS", "FAIL", "WARN"]
    colour_map = {"PASS": "#2ecc71", "FAIL": "#e74c3c", "WARN": "#f39c12"}
    colours    = [colour_map[s] for s in statuses]
    v1 = [c1.get(s, 0) for s in statuses]
    v2 = [c2.get(s, 0) for s in statuses]

    x     = np.arange(len(statuses))
    width = 0.35
    fig, ax = plt.subplots(figsize=(7, 5))
    b1 = ax.bar(x - width / 2, v1, width, label=label1,
                color=colours, alpha=0.9, edgecolor="white")
    b2 = ax.bar(x + width / 2, v2, width, label=label2,
                color=colours, alpha=0.55, edgecolor="white", hatch="//")
    for bar in list(b1) + list(b2):
        h = bar.get_height()
        if h > 0:
            ax.text(bar.get_x() + bar.get_width() / 2,
                    h + 0.1, str(int(h)),
                    ha="center", va="bottom", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(statuses, fontsize=12)
    ax.set_ylabel("Number of checks", fontsize=11)
    ax.set_title("Validation Check Results by Species",
                 fontsize=13, fontweight="bold")
    ax.legend()
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  ✓  Plot saved → {out_path}")


def plot_ortholog_conservation(
    ortho_df: pd.DataFrame,
    out_path: Path,
    label1: str = "human",
    label2: str = "otherSpecies",
) -> None:
    """
    Horizontal bar chart of mean Jaccard conservation score per category,
    computed across all complete ortholog pairs.
    """
    if not _check_mpl():
        return
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.use("Agg")

    complete = ortho_df[ortho_df["pair_complete"]]
    if len(complete) == 0:
        print("  ⚠  No complete ortholog pairs — skipping conservation plot")
        return

    cats   = list(CONSERVATION_GROUPS.keys())
    means  = []
    for cat in cats:
        col = f"jaccard_{cat}"
        means.append(
            round(complete[col].mean(), 3) if col in complete.columns else 0.0
        )

    colours = ["#2ecc71" if v >= 0.7 else
               "#f39c12" if v >= 0.4 else
               "#e74c3c" for v in means]

    fig, ax = plt.subplots(figsize=(8, max(4, len(cats) * 0.6)))
    bars = ax.barh(cats, means, color=colours, edgecolor="white", height=0.55)
    for bar, val in zip(bars, means):
        ax.text(val + 0.01, bar.get_y() + bar.get_height() / 2,
                f"{val:.3f}", va="center", fontsize=9, fontweight="bold")
    ax.set_xlim(0, 1.15)
    ax.set_xlabel("Mean Jaccard conservation score", fontsize=11)
    ax.set_title(
        f"Ortholog Conservation: {label1} → {label2}\n"
        f"({len(complete):,} complete pairs)",
        fontsize=12, fontweight="bold"
    )
    ax.axvline(0.7, color="#2ecc71", linestyle="--", linewidth=0.8,
               label="high (≥0.7)")
    ax.axvline(0.4, color="#f39c12", linestyle="--", linewidth=0.8,
               label="medium (≥0.4)")
    ax.legend(fontsize=8, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  ✓  Plot saved → {out_path}")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 12 — Main pipeline orchestrator
# ══════════════════════════════════════════════════════════════════════════════

def run_pipeline(
    path_human:  Path,
    path_other:  Path,
    path_ann:    Path,
    out_dir:     Path,
    path_ortho:  Path | None = None,
    label1:      str = "human",
    label2:      str = "otherSpecies",
    make_plots:  bool = False,
) -> None:
    """
    Orchestrate the full cross-species transfer validation pipeline.

    Steps
    ─────
    1  Load all input files and build CV lookup
    2  Validate each species occurrence table independently (18 checks each)
    3  Build missing term report with closest-match hints
    4  Compute per-species summary metrics
    5  Compute transfer metrics (schema/vocab/rule portability, conservation)
    6  Compute overlap and Jaccard metrics per category
    7  Build side-by-side comparison report
    8  Run ortholog conservation analysis (if ortholog file provided)
    9  Generate optional plots
    10 Print final summary and save all outputs
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Banner ────────────────────────────────────────────────────────────────
    print_banner("CROSS-SPECIES ONTOLOGY TRANSFER VALIDATION PIPELINE")
    print(f"  Reference species  : {path_human.name}  [{label1}]")
    print(f"  Target species     : {path_other.name}  [{label2}]")
    print(f"  Annotation registry: {path_ann.name}")
    if path_ortho:
        print(f"  Ortholog mapping   : {path_ortho.name}")
    print(f"  Output directory   : {out_dir}")

    # ── Step 1: Load inputs ───────────────────────────────────────────────────
    print_section("1 / 10   Loading inputs")
    df1  = load_table(path_human)
    print(f"  ✓  {label1}: {df1.shape[0]:,} rows × {df1.shape[1]} cols")
    df2  = load_table(path_other)
    print(f"  ✓  {label2}: {df2.shape[0]:,} rows × {df2.shape[1]} cols")
    ann  = load_table(path_ann)
    print(f"  ✓  Annotation registry: {len(ann):,} terms")
    cv   = build_cv_lookup(ann)
    print("  ✓  CV lookup: " +
          " | ".join(f"{k}={len(v)//2}" for k, v in cv.items()))

    ortho_df: pd.DataFrame | None = None
    if path_ortho is not None:
        ortho_df = load_ortholog_mapping(path_ortho)
        print(f"  ✓  Ortholog mapping: {len(ortho_df):,} pairs")

    # ── Step 2: Validate both species tables ─────────────────────────────────
    print_section("2 / 10   Validating species occurrence tables")
    print(f"\n── {label1} ──")
    res1, inv1, unk1 = validate_occurrence_table(df1, label1, cv)
    print(f"\n── {label2} ──")
    res2, inv2, unk2 = validate_occurrence_table(df2, label2, cv)

    save_tsv(pd.DataFrame(res1),
             out_dir / "validation_report_human.tsv",
             "validation checks")
    save_tsv(pd.DataFrame(res2),
             out_dir / "validation_report_otherSpecies.tsv",
             "validation checks")

    all_invalid = pd.DataFrame(inv1 + inv2)
    if not all_invalid.empty:
        meta = [c for c in ["species", "check_id", "reason", "source_row"]
                if c in all_invalid.columns]
        other_cols = [c for c in all_invalid.columns if c not in meta]
        save_tsv(all_invalid[meta + other_cols],
                 out_dir / "invalid_rows.tsv", "flagged rows")
    else:
        save_tsv(pd.DataFrame(), out_dir / "invalid_rows.tsv")

    all_unknown = pd.DataFrame(unk1 + unk2)
    save_tsv(
        all_unknown if not all_unknown.empty else pd.DataFrame(),
        out_dir / "unknown_terms_report.tsv",
        "tokens not in registry",
    )

    # ── Step 3: Missing term report ───────────────────────────────────────────
    print_section("3 / 10   Building missing term report")
    missing_df = build_missing_term_report(
        {label1: df1, label2: df2}, cv
    )
    save_tsv(missing_df, out_dir / "missing_term_report.tsv",
             "unmatched CV tokens with hints")

    # ── Step 4: Per-species summary metrics ───────────────────────────────────
    print_section("4 / 10   Computing per-species summary metrics")
    m1  = compute_species_metrics(df1, label1, cv)
    m2  = compute_species_metrics(df2, label2, cv)
    save_tsv(pd.DataFrame([m1, m2]),
             out_dir / "species_summary_metrics.tsv", "summary")

    # ── Step 5: Transfer metrics ──────────────────────────────────────────────
    print_section("5 / 10   Computing cross-species transfer metrics")
    transfer = compute_transfer_metrics(df1, df2, cv, label1, label2)
    save_tsv(pd.DataFrame([transfer]),
             out_dir / "transfer_metrics_cross_species.tsv", "transfer")

    print(f"\n  Schema reuse rate              : {transfer['schema_reuse_rate']:.1%}")
    print(f"  Vocabulary reuse rate          : {transfer['vocab_reuse_rate']:.1%}")
    print(f"  Rule portability rate          : {transfer['rule_portability_rate']:.1%}")
    print(f"  Mapping completeness           : {transfer['mapping_completeness']:.1%}")
    print(f"  Ontology term conservation     : {transfer['ontology_term_conservation_rate']:.1%}")
    print(f"  Species-specific annotation    : {transfer['species_specific_annotation_rate']:.1%}")
    print(f"  Terms conserved (both species) : {transfer['n_terms_conserved']}")
    print(f"  New terms required             : {transfer['n_new_terms_required']}")
    print(f"  Registry terms unused          : {transfer['n_registry_terms_unused']}")

    # ── Step 6: Overlap and Jaccard metrics ───────────────────────────────────
    print_section("6 / 10   Computing overlap and Jaccard metrics")
    overlap_df = compute_overlap_metrics(df1, df2, label1, label2)
    save_tsv(overlap_df,
             out_dir / "overlap_metrics_cross_species.tsv",
             "Jaccard per category")
    print()
    for _, row in overlap_df.iterrows():
        print(f"  {row['category']:<35}  "
              f"Jaccard = {row['jaccard']:.3f}  "
              f"(shared {row['shared']})")

    # ── Step 7: Comparison report ─────────────────────────────────────────────
    print_section("7 / 10   Building comparison report")
    comp_df = build_comparison_report(df1, df2, label1, label2, cv)
    save_tsv(comp_df,
             out_dir / "comparison_report_cross_species.tsv",
             "side-by-side comparison")

    # ── Step 8: Ortholog conservation analysis ────────────────────────────────
    if ortho_df is not None:
        print_section("8 / 10   Running ortholog conservation analysis")
        ortho_report = compute_ortholog_conservation(
            df1, df2, ortho_df, label1, label2
        )
        save_tsv(ortho_report,
                 out_dir / "ortholog_conservation_report.tsv",
                 "per-pair conservation scores")

        # Proteins with no ortholog
        no_ortho_human = (
            set(df1["protein_accession"].dropna().str.upper())
            - set(ortho_df["human_accession"].str.upper())
            if "protein_accession" in df1.columns else set()
        )
        no_ortho_other = (
            set(df2["protein_accession"].dropna().str.upper())
            - set(ortho_df["other_accession"].str.upper())
            if "protein_accession" in df2.columns else set()
        )
        if no_ortho_human or no_ortho_other:
            unmapped_rows: list[dict] = []
            for acc in sorted(no_ortho_human):
                unmapped_rows.append({
                    "species": label1,
                    "protein_accession": acc,
                    "note": "No ortholog in mapping file",
                })
            for acc in sorted(no_ortho_other):
                unmapped_rows.append({
                    "species": label2,
                    "protein_accession": acc,
                    "note": "No ortholog in mapping file",
                })
            save_tsv(
                pd.DataFrame(unmapped_rows),
                out_dir / "proteins_without_ortholog.tsv",
                "unmapped proteins",
            )
            print(f"  {label1} proteins without ortholog : {len(no_ortho_human)}")
            print(f"  {label2} proteins without ortholog : {len(no_ortho_other)}")
    else:
        print_section("8 / 10   Ortholog analysis skipped (no --ortho file)")
        save_tsv(
            pd.DataFrame(columns=[
                "human_accession", "other_accession",
                "gene_name", "ortholog_type", "pair_complete",
                "mean_conservation_score",
            ]),
            out_dir / "ortholog_conservation_report.tsv",
            "empty — no ortholog file provided",
        )
        ortho_report = None

    # ── Step 9: Optional plots ────────────────────────────────────────────────
    if make_plots:
        print_section("9 / 10   Generating publication-quality figures (7 figures)")

        # Fig 1 — Validation summary (3-panel)
        plot_validation_summary(
            res1, res2,
            out_dir / "validation_summary_plot.png",
            label1, label2,
        )
        # Fig 2 — Species comparison overview (6-panel)
        plot_species_comparison(
            df1, df2, comp_df, transfer, label1, label2,
            out_dir / "species_comparison_plot.png",
        )
        # Fig 3 — Jaccard overlap (lollipop + set-size bars)
        plot_jaccard_heatmap(
            overlap_df, out_dir / "jaccard_overlap_plot.png",
            label1, label2,
        )
        # Fig 4 — Transfer metrics (radar + term conservation)
        plot_transfer_metrics(
            transfer,
            out_dir / "transfer_metrics_plot.png",
            label1, label2,
        )
        # Fig 5 — Domain conservation profile
        plot_domain_conservation(
            df1, df2, label1, label2,
            out_dir / "domain_conservation_plot.png",
        )
        # Fig 6 — Function-tag cross-species profile
        plot_function_profile(
            df1, df2, label1, label2,
            out_dir / "function_profile_plot.png",
        )
        # Fig 7 — Ortholog conservation (only when --ortho provided)
        if ortho_report is not None:
            plot_ortholog_conservation(
                ortho_report,
                out_dir / "ortholog_conservation_plot.png",
                label1, label2,
            )
    else:
        print_section("9 / 10   Plots skipped (pass --plots to enable)")

    # ── Step 10: Final summary ────────────────────────────────────────────────
    print_section("10 / 10   Pipeline complete")
    total1 = len(res1)
    pass1  = sum(1 for r in res1 if r["status"] == "PASS")
    fail1  = sum(1 for r in res1 if r["status"] == "FAIL")
    total2 = len(res2)
    pass2  = sum(1 for r in res2 if r["status"] == "PASS")
    fail2  = sum(1 for r in res2 if r["status"] == "FAIL")

    print(f"\n  {label1:<25}  {total1} checks | {pass1} PASS | {fail1} FAIL")
    print(f"  {label2:<25}  {total2} checks | {pass2} PASS | {fail2} FAIL")
    print(f"  Flagged rows          : {len(inv1) + len(inv2):,}")
    print(f"  Missing terms         : {len(missing_df):,}")
    print(f"  All outputs           → {out_dir}/")
    if fail1 + fail2 > 0:
        print(
            "\n  ⚠  Some checks FAILED — review:\n"
            f"     • {out_dir}/invalid_rows.tsv\n"
            f"     • {out_dir}/missing_term_report.tsv"
        )
    print()


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 13 — CLI entry point
# ══════════════════════════════════════════════════════════════════════════════

def main() -> None:
    """
    Entry point — supports both CLI arguments and hardcoded paths (PyCharm mode).

    CLI usage:
      python cross_species_transfer_validation_pipeline.py \\
          --human  path/to/human.tsv \\
          --other  path/to/mouse.tsv \\
          --ann    path/to/Annotation_table.xlsx \\
          --out    output_dir \\
          --species1 human --species2 mouse \\
          --plots

    PyCharm mode (no args):
      Edit the PATH_* variables below and run directly.
    """
    import sys

    # ── PyCharm / direct-run configuration ────────────────────────────────────
    # Edit these paths when running without CLI arguments.
    PATH_HUMAN  = Path(r"C:\Users\32045\PycharmProjects\PythonProject\cross_testing\c_AMP_pathway_MERGED_domain_occurrences_WITH_IDS.tsv")
    PATH_OTHER  = Path(r"C:\Users\32045\PycharmProjects\PythonProject\cross_testing\c_AMP_pathway_mouse_MERGED_domain_occurrences_WITH_IDS.tsv")
    PATH_ANN    = Path(r"C:\Users\32045\PycharmProjects\PythonProject\cross_testing\Annotation_table_cAMP_mouse_improved.xlsx")
    PATH_ORTHO  = None
    OUT_DIR     = Path(r"output_cross_species_cAMP_human_vs_mouse")
    SPECIES1    = "cAMP_human"
    SPECIES2    = "cAMP_mouse"
    MAKE_PLOTS  = True

    # ── CLI argument parsing ───────────────────────────────────────────────────
    if len(sys.argv) > 1:
        parser = argparse.ArgumentParser(
            description="Cross-Species Ontology Transfer Validation Pipeline",
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )
        parser.add_argument("--human",    required=True, help="Reference species TSV")
        parser.add_argument("--other",    required=True, help="Target species TSV")
        parser.add_argument("--ann",      required=True, help="Annotation table XLSX")
        parser.add_argument("--out",      default="output_cross_species",
                            help="Output directory")
        parser.add_argument("--ortho",    default=None,
                            help="Ortholog mapping TSV (optional)")
        parser.add_argument("--species1", default="human",
                            help="Label for reference species (default: human)")
        parser.add_argument("--species2", default="otherSpecies",
                            help="Label for target species (default: otherSpecies)")
        parser.add_argument("--plots",    action="store_true",
                            help="Generate publication-quality figures")
        args = parser.parse_args()
        PATH_HUMAN = Path(args.human)
        PATH_OTHER = Path(args.other)
        PATH_ANN   = Path(args.ann)
        PATH_ORTHO = Path(args.ortho) if args.ortho else None
        OUT_DIR    = Path(args.out)
        SPECIES1   = args.species1
        SPECIES2   = args.species2
        MAKE_PLOTS = args.plots

    run_pipeline(
        path_human  = PATH_HUMAN,
        path_other  = PATH_OTHER,
        path_ann    = PATH_ANN,
        out_dir     = OUT_DIR,
        path_ortho  = PATH_ORTHO,
        label1      = SPECIES1,
        label2      = SPECIES2,
        make_plots  = MAKE_PLOTS,
    )


if __name__ == "__main__":
    main()
