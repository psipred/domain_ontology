"""
02_id_injection.py
────────────────────────────────────────────────────────────────────────────────
Stage 3 of the domain-ontology pipeline.

Reads the DOP/DOF/DOT annotation table and the merged domain-occurrence TSV
produced by 01_data_and_merge.py, then replaces every plain-text vocabulary
token in the annotated columns with a stable "ID:XXXXXX label" string, e.g.:

    "Cterminal"       →  "DOP:000003 Cterminal"
    "kinase"          →  "DOF:000039 kinase"
    "cAMP"            →  "DOT:000002 cAMP"

Columns that receive ID injection
    positional_category       – DOP terms  (positional_or_topology_term)
    copy_number_property      – DOP terms  (CopyNumberCategory)
    topology_context          – DOP terms  (positional_or_topology_term)
    proximity_categories      – DOP terms  (positional_or_topology_term)
    function_tags             – DOF terms  (functional_term)
    binding_partner           – DOT terms  (binding_target_term)
    binding_partner_specific  – DOT terms  (binding_target_term)

Columns passed through unchanged (free-text or computed — no ID mapping)
    note, reference, domain_sequence, cooccurring_domain_occurrences,
    start_residue, end_residue, and all other columns

Optional: carry manual annotations forward from a previously-corrected file
    If you have already filled in function_tags / binding_partner in an
    older version of the occurrence table, pass it via --carry-over so those
    values are merged into the new pipeline output before injection.
    Matching is done on protein_accession + domain_type_id + start + end.

Each cell may contain multiple tokens separated by ";".
Tokens not found in the annotation table are left unchanged and written to
    unmapped_tokens.txt  grouped by column, so you know exactly which
    annotation-table row to add next.

Usage
    # Minimum required arguments (output name derived automatically)
    py 02_id_injection.py ^
        --annotation Annotation_table_human_cAMP.xlsx ^
        --input      cAMP_human_MERGED_domain_occurrences.tsv

    # Specify output path explicitly
    py 02_id_injection.py ^
        --annotation Annotation_table_human_cAMP.xlsx ^
        --input      cAMP_human_MERGED_domain_occurrences.tsv ^
        --output     cAMP_human_MERGED_domain_occurrences_WITH_IDS.tsv

    # Carry manual annotations over from a previously-corrected file
    py 02_id_injection.py ^
        --annotation Annotation_table_human_cAMP.xlsx ^
        --input      cAMP_human_MERGED_domain_occurrences.tsv ^
        --carry-over c_AMP_pathway_MERGED_domain_occurrences.tsv

    # Strict mode: exit with error if any token is unmapped
    py 02_id_injection.py ^
        --annotation Annotation_table_human_cAMP.xlsx ^
        --input      cAMP_human_MERGED_domain_occurrences.tsv ^
        --strict

Pipeline position
    01_data_and_merge.py   → produces *_MERGED_domain_occurrences.tsv
    02_id_injection.py     → produces *_MERGED_domain_occurrences_WITH_IDS.tsv  ← YOU ARE HERE
    03_owl_translation.py  → produces OWL / Turtle population file
────────────────────────────────────────────────────────────────────────────────
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

# ── Output naming convention ─────────────────────────────────────────────────
# --output is optional: if omitted, the output file is derived automatically
# from --input by inserting "_WITH_IDS" before the .tsv extension, e.g.:
#   cAMP_human_MERGED_domain_occurrences.tsv
#     → cAMP_human_MERGED_domain_occurrences_WITH_IDS.tsv
# --unmapped defaults to unmapped_tokens.txt in the same directory as --input.

# ── Columns that receive ID injection ─────────────────────────────────────────
# Note: 'note', 'reference', 'domain_sequence', 'cooccurring_domain_occurrences',
# 'start_residue', 'end_residue' are intentionally excluded — free-text or
# computed fields that must not be transformed.
TARGET_COLUMNS: List[str] = [
    "positional_category",
    "copy_number_property",
    "topology_context",
    "proximity_categories",
    "function_tags",
    "binding_partner",
    "binding_partner_specific",
]

# ── Which annotation-table term_types serve each column ───────────────────────
TERM_TYPE_FOR_COLUMN: Dict[str, List[str]] = {
    "positional_category"      : ["positional_or_topology_term"],
    "copy_number_property"     : ["CopyNumberCategory"],
    "topology_context"         : ["positional_or_topology_term"],
    "proximity_categories"     : ["positional_or_topology_term"],
    "function_tags"            : ["functional_term"],
    "binding_partner"          : ["binding_target_term"],
    "binding_partner_specific" : ["binding_target_term"],
}

# ── Columns whose manual content can be carried over from an older file ───────
CARRY_OVER_COLUMNS: List[str] = [
    "function_tags",
    "binding_partner",
    "binding_partner_specific",
    "note",
    "reference",
]

# ── Legacy column-name aliases (old typo → correct name) ─────────────────────
LEGACY_ALIASES: Dict[str, str] = {
    "binding_parnter"          : "binding_partner",
    "binding_parnter_specific" : "binding_partner_specific",
}

# ── Match key for carry-over merge ────────────────────────────────────────────
MERGE_KEY: List[str] = ["protein_accession", "domain_type_id", "start", "end"]

TOKEN_DELIMITER = ";"


# ══════════════════════════════════════════════════════════════════════════════
# Helpers
# ══════════════════════════════════════════════════════════════════════════════

def _find_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    """Return the first matching column from *candidates* (case-insensitive).
    Returns None if none found."""
    lower_map = {c.lower(): c for c in df.columns}
    for name in candidates:
        if name.lower() in lower_map:
            return lower_map[name.lower()]
    return None


def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, low_memory=False, encoding="latin-1").fillna("")


def _normalise_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Rename any legacy-typo column names to their corrected equivalents."""
    rename = {old: new for old, new in LEGACY_ALIASES.items() if old in df.columns}
    if rename:
        df = df.rename(columns=rename)
        for old, new in rename.items():
            print(f"  [RENAME] '{old}' → '{new}'  (legacy typo correction)")
    return df


# ══════════════════════════════════════════════════════════════════════════════
# Step 1 — Build lookup dictionary from the annotation table
# ══════════════════════════════════════════════════════════════════════════════

def build_lookup(annotation_path: str) -> Dict[str, Dict[str, str]]:
    """
    Return  lookup[column_name][token_lowercase] = "ID:XXXXXX PrimaryName"

    Every alias (Name, canonical_name, synonym columns) is indexed so that
    alternate spellings all resolve to the same ID string.
    Lookup is always case-insensitive.
    """
    ann = pd.read_excel(annotation_path, dtype=str).fillna("")

    # Normalise column names to lowercase for robust access
    ann.columns = [c.strip().lower() for c in ann.columns]

    for required in ("id", "term_type"):
        if required not in ann.columns:
            raise ValueError(
                f"Annotation table is missing required column '{required}'.\n"
                f"  Found columns: {list(ann.columns)}\n"
                f"  Check that the file is the correct Annotation_table_human_cAMP.xlsx."
            )

    ann = ann[ann["id"].str.strip() != ""]  # skip rows with no ID

    lookup: Dict[str, Dict[str, str]] = {col: {} for col in TARGET_COLUMNS}

    for _, row in ann.iterrows():
        term_id   = row["id"].strip()
        term_type = row.get("term_type", "").strip()

        aliases: List[str] = []
        for field in ("name", "canonical_name", "synonym"):
            val = row.get(field, "")
            if val.strip():
                aliases.append(val.strip())

        if not aliases:
            continue

        primary  = aliases[0]
        id_label = f"{term_id} {primary}"

        for col, allowed_types in TERM_TYPE_FOR_COLUMN.items():
            if term_type in allowed_types:
                for alias in aliases:
                    lookup[col][alias.lower()] = id_label

    return lookup


# ══════════════════════════════════════════════════════════════════════════════
# Step 2 — Carry over manual annotations from an older corrected file
# ══════════════════════════════════════════════════════════════════════════════

def carry_over_annotations(
        target: pd.DataFrame,
        source_path: Path,
) -> pd.DataFrame:
    """
    Merge manual annotation columns from *source_path* (an older, hand-corrected
    occurrence file) into *target* (the new pipeline output).

    Matching is done on MERGE_KEY (protein_accession + domain_type_id + start + end).
    Only rows whose function_tags / binding_partner are currently empty in *target*
    are updated — existing target values are never overwritten.
    """
    print(f"\nCarrying over manual annotations from {source_path.name} …")
    src = _read_tsv(source_path)
    src = _normalise_columns(src)

    # Check that source has at least one carry-over column worth using
    src_useful = [c for c in CARRY_OVER_COLUMNS if c in src.columns
                  and (src[c].str.strip() != "").any()]
    if not src_useful:
        print(f"  [SKIP] Source file has no filled carry-over columns — nothing to merge.")
        return target

    # Check merge keys exist in both files
    missing_key_src = [k for k in MERGE_KEY if k not in src.columns]
    missing_key_tgt = [k for k in MERGE_KEY if k not in target.columns]
    if missing_key_src or missing_key_tgt:
        print(f"  [WARN] Merge key columns missing — cannot carry over.")
        if missing_key_src:
            print(f"         Missing in source: {missing_key_src}")
        if missing_key_tgt:
            print(f"         Missing in target: {missing_key_tgt}")
        return target

    # Build source lookup: key_tuple → {col: value}
    src_lookup: Dict[Tuple, Dict[str, str]] = {}
    for _, row in src.iterrows():
        key = tuple(str(row.get(k, "")).strip() for k in MERGE_KEY)
        src_lookup[key] = {c: str(row.get(c, "")).strip() for c in src_useful}

    merged_counts: Dict[str, int] = {c: 0 for c in src_useful}

    for idx, row in target.iterrows():
        key = tuple(str(row.get(k, "")).strip() for k in MERGE_KEY)
        src_vals = src_lookup.get(key)
        if src_vals is None:
            continue
        for col, val in src_vals.items():
            if col not in target.columns:
                target[col] = ""
            # Only fill if the target cell is currently empty
            if val and str(target.at[idx, col]).strip() == "":
                target.at[idx, col] = val
                merged_counts[col] += 1

    print(f"  Matched {len(src_lookup)} source rows against {len(target)} target rows.")
    for col, count in merged_counts.items():
        if count:
            print(f"  ✓ {col:<32} {count:>4} cells filled from source")
    if not any(merged_counts.values()):
        print(f"  [INFO] No new values added (target cells were already filled or no key matches).")

    return target


# ══════════════════════════════════════════════════════════════════════════════
# Step 3 — Transform a single cell
# ══════════════════════════════════════════════════════════════════════════════

def transform_cell(
        cell_value,
        mapping: Dict[str, str],
        unmapped: Set[Tuple[str, str]],
        column_name: str,
) -> str:
    """
    Split the cell on TOKEN_DELIMITER, map each token, and rejoin.

    Idempotent: tokens already in "XXX:NNNNNN label" format are passed through
    unchanged — safe to run the script twice on the same file.
    Unrecognised tokens are passed through unchanged AND logged with their
    column name so the unmapped-token report is immediately actionable.
    """
    if pd.isna(cell_value) or str(cell_value).strip() == "":
        return cell_value

    tokens = [t.strip() for t in str(cell_value).split(TOKEN_DELIMITER) if t.strip()]
    result_tokens = []

    for token in tokens:
        # Idempotent check: already looks like "DOP:000003 Cterminal"
        if len(token) >= 4 and ":" in token[:7] and token.split(":")[0].isupper():
            result_tokens.append(token)
            continue

        mapped = mapping.get(token.lower())
        if mapped:
            result_tokens.append(mapped)
        else:
            result_tokens.append(token)
            unmapped.add((column_name, token))

    return TOKEN_DELIMITER.join(result_tokens)


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="02_id_injection.py",
        description=(
            "Stage 3: inject DOP/DOF/DOT IDs into the domain-occurrence TSV.\n"
            "Reads the annotation table, maps vocabulary tokens to stable IDs,\n"
            "and writes *_WITH_IDS.tsv ready for OWL translation."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "EXAMPLES\n"
            "  # Minimum — output name is derived automatically\n"
            "  py 02_id_injection.py \\\n"
            "      --annotation Annotation_table_human_cAMP.xlsx \\\n"
            "      --input  cAMP_human_MERGED_domain_occurrences.tsv\n\n"
            "  # With explicit output path\n"
            "  py 02_id_injection.py \\\n"
            "      --annotation Annotation_table_human_cAMP.xlsx \\\n"
            "      --input  cAMP_human_MERGED_domain_occurrences.tsv \\\n"
            "      --output cAMP_human_MERGED_domain_occurrences_WITH_IDS.tsv\n\n"
            "  # Carry over manual annotations from an older corrected file\n"
            "  py 02_id_injection.py \\\n"
            "      --annotation Annotation_table_human_cAMP.xlsx \\\n"
            "      --input  cAMP_human_MERGED_domain_occurrences.tsv \\\n"
            "      --carry-over c_AMP_pathway_MERGED_domain_occurrences.tsv\n\n"
            "  # Strict mode (fail if any tokens unmapped)\n"
            "  py 02_id_injection.py \\\n"
            "      --annotation Annotation_table_human_cAMP.xlsx \\\n"
            "      --input  cAMP_human_MERGED_domain_occurrences.tsv \\\n"
            "      --strict\n"
        ),
    )
    ap.add_argument("--annotation", required=True, metavar="XLSX",
                    help="Annotation table Excel file  "
                         "(e.g. Annotation_table_human_cAMP.xlsx)")
    ap.add_argument("--input",      required=True, metavar="TSV",
                    help="Domain-occurrence TSV produced by 01_data_and_merge.py  "
                         "(e.g. cAMP_human_MERGED_domain_occurrences.tsv)")
    ap.add_argument("--output",     default=None,  metavar="TSV",
                    help="Output TSV path  "
                         "(default: <input>_WITH_IDS.tsv, derived automatically)")
    ap.add_argument("--unmapped",   default=None,  metavar="TXT",
                    help="Unmapped-token log  "
                         "(default: unmapped_tokens.txt beside the output file)")
    ap.add_argument("--carry-over",  default=None,               metavar="TSV",
                    help="Older corrected occurrence file — manual annotations "
                         "(function_tags, binding_partner, note, reference) are "
                         "merged into the new file before injection")
    ap.add_argument("--strict", action="store_true",
                    help="Exit with error code 1 if any tokens remain unmapped")
    return ap


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main(argv: Optional[List[str]] = None) -> None:
    args = build_arg_parser().parse_args(argv)

    ann_path  = Path(args.annotation)
    in_path   = Path(args.input)

    # Derive output path: insert "_WITH_IDS" before the file extension
    if args.output:
        out_path = Path(args.output)
    else:
        stem = in_path.stem          # e.g. "cAMP_human_MERGED_domain_occurrences"
        out_path = in_path.parent / f"{stem}_WITH_IDS.tsv"

    # Derive unmapped log path: same directory as output
    if args.unmapped:
        unmap_path = Path(args.unmapped)
    else:
        unmap_path = out_path.parent / "unmapped_tokens.txt"

    carry_path = Path(args.carry_over) if args.carry_over else None

    # ── Validate inputs ───────────────────────────────────────────────────
    for p in (ann_path, in_path):
        if not p.exists():
            print(f"[ERROR] File not found: {p}", file=sys.stderr)
            print(f"        Pass the correct path with --annotation / --input,\n"
                  f"        or place the file in the current working directory.",
                  file=sys.stderr)
            sys.exit(1)

    if carry_path and not carry_path.exists():
        print(f"[WARN] --carry-over file not found: {carry_path} — skipping carry-over step.")
        carry_path = None

    # ── Step 1: Build lookup ──────────────────────────────────────────────
    print("=" * 60)
    print("02_id_injection.py — DOP/DOF/DOT ID injection")
    print("=" * 60)
    print(f"\n[1] Reading annotation table: {ann_path.name}")
    try:
        lookup = build_lookup(str(ann_path))
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)

    print(f"  {'Column':<32} {'Aliases indexed':>15}")
    print(f"  {'-'*32} {'-'*15}")
    for col, mapping in lookup.items():
        print(f"  {col:<32} {len(mapping):>15,}")

    # ── Step 2: Read and normalise data ───────────────────────────────────
    print(f"\n[2] Reading occurrence table: {in_path.name}")
    data = _read_tsv(in_path)
    data = _normalise_columns(data)
    print(f"  {len(data):,} rows × {len(data.columns)} columns")

    # ── Step 3: Carry over manual annotations (optional) ─────────────────
    if carry_path:
        data = carry_over_annotations(data, carry_path)

    # ── Step 4: Identify columns to transform ────────────────────────────
    print(f"\n[3] Identifying target columns …")
    present_cols: List[str] = []
    for col in TARGET_COLUMNS:
        if col in data.columns:
            present_cols.append(col)
            print(f"  ✓ {col}")
        else:
            print(f"  [SKIP] '{col}' not in data — skipping "
                  f"(normal if column is new/empty in this pipeline run)")

    if not present_cols:
        print("[ERROR] None of the target columns were found.", file=sys.stderr)
        sys.exit(1)

    # ── Step 5: Inject IDs ────────────────────────────────────────────────
    print(f"\n[4] Injecting IDs …")
    print(f"  {'Column':<32} {'Cells mapped':>12}  {'Non-empty':>10}")
    print(f"  {'-'*32} {'-'*12}  {'-'*10}")

    unmapped: Set[Tuple[str, str]] = set()  # (column, token)

    for col in present_cols:
        col_mapping = lookup[col]
        orig = data[col].copy()

        data[col] = data[col].apply(
            lambda v, _m=col_mapping, _c=col: transform_cell(v, _m, unmapped, _c)
        )

        changed      = int((data[col] != orig).sum())
        non_empty    = int((orig.str.strip() != "").sum())
        print(f"  {col:<32} {changed:>12,}  {non_empty:>10,}")

    # ── Step 6: Save output ───────────────────────────────────────────────
    out_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(out_path, sep="\t", index=False, encoding="utf-8")
    print(f"\n[5] Saved → {out_path}")
    print(f"    {len(data):,} rows × {len(data.columns)} columns")

    # ── Step 7: Report unmapped tokens ────────────────────────────────────
    if unmapped:
        by_col: Dict[str, List[str]] = {}
        for col, token in sorted(unmapped):
            by_col.setdefault(col, []).append(token)

        lines = [
            "Unmapped tokens — not found in annotation table",
            "=" * 52,
            f"Total: {len(unmapped)} unique token(s)",
            "",
            "To fix: add a row to the annotation table for each",
            "token, then re-run 02_id_injection.py.",
            "",
        ]
        for col, tokens in sorted(by_col.items()):
            lines.append(f"[{col}]  ({len(tokens)} token(s))")
            for t in sorted(tokens):
                lines.append(f"    {t}")
            lines.append("")

        unmap_path.write_text("\n".join(lines), encoding="utf-8")
        print(f"\n⚠  {len(unmapped)} unmapped token(s) — see {unmap_path}")

        for col, tokens in sorted(by_col.items()):
            print(f"   [{col}]")
            for t in sorted(tokens):
                print(f"     {t!r}")

        if args.strict:
            print("\n[STRICT] Exiting with code 1 — unmapped tokens present.",
                  file=sys.stderr)
            sys.exit(1)
    else:
        print("\n✓  All tokens mapped — no unmapped tokens.")

    print("\n[DONE] Next step: run 03_owl_translation.py on the WITH_IDS file.")


if __name__ == "__main__":
    main()
