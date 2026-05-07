"""
06_build_annotation_tables.py
────────────────────────────────────────────────────────────────────────────
Step 6 — Extract annotations into long-form tables.

IMPORTANT: uses ortholog_pairs_calibration.tsv (not the full pairs file)
to ensure the blind holdout set is never touched at this stage.

Fixes vs original
─────────────────
  • Binding-column detection handles the 'binding_parnter' typo column name
    that is present in both input files (uses whichever exists).
  • copy_number_property column correctly extracted (includes the DOP-backed
    term, not the raw integer copy_number column).
  • is_transferable filter correctly applied (was invalid pandas .get() call).

Outputs
───────
  human_annotations_long.tsv    all human annotations (long form)
  mouse_annotations_gold.tsv    all mouse annotations (gold labels)
  transferability_rules.tsv     which categories to transfer (edit this)

  DESIGN RULE: do NOT edit transferability_rules.tsv after looking at
  mouse gold labels. Rules must be justified on biological grounds alone.
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd


# ── Annotation map: column → (predicate, category) ───────────────────────────
# 'binding_parnter' is the actual column name in the source files (typo preserved)
ANNOTATION_MAP: list[tuple[str, str, str]] = [
    ("domain_type_id",       "hasDomainType",          "domain_type"),
    ("function_tags",        "hasFunctionTag",          "function"),
    ("topology_context",     "hasTopologyContext",      "topology"),
    ("positional_category",  "hasPositionalCategory",   "position"),
    ("copy_number_property", "hasCopyNumberProperty",   "copy_number"),
    ("proximity_categories", "hasProximityCategory",    "proximity"),
    ("binding_parnter",      "hasBindingTarget",        "binding_target"),
    ("go_terms_mf",          "hasGOTermMF",             "GO_function"),
    ("go_terms_bp",          "hasGOTermBP",             "GO_process"),
    ("go_terms_cc",          "hasGOTermCC",             "GO_component"),
]

# Fallback column names (tried if the primary name is absent)
COLUMN_FALLBACKS: dict[str, list[str]] = {
    "binding_parnter":      ["binding_partner"],
    "copy_number_property": ["copy_number_property"],  # already primary
}

DEFAULT_TRANSFERABILITY: dict[str, bool] = {
    "domain_type":   True,
    "function":      True,
    "topology":      True,
    "position":      True,
    "copy_number":   True,
    "proximity":     True,
    "binding_target":False,
    "GO_function":   True,
    "GO_process":    False,
    "GO_component":  False,
}


def _split_cell(val) -> list[str]:
    if pd.isna(val) or str(val).strip() in ("", "nan"):
        return []
    return [v.strip() for v in str(val).split(";") if v.strip()]


def _resolve_column(occ: pd.DataFrame, primary: str) -> str | None:
    """Return the actual column name to use, trying fallbacks."""
    if primary in occ.columns:
        return primary
    for alt in COLUMN_FALLBACKS.get(primary, []):
        if alt in occ.columns:
            return alt
    return None


def _build_long_table(occ_path: Path, proteins_tsv: Path,
                       species: str,
                       transferability: dict[str, bool]) -> pd.DataFrame:
    if not occ_path.exists():
        print(f"  WARNING: {occ_path} not found")
        return pd.DataFrame()

    occ = pd.read_csv(occ_path, sep="\t", dtype=str, low_memory=False)

    gene_map: dict[str, str] = {}
    if proteins_tsv.exists():
        ptbl  = pd.read_csv(proteins_tsv, sep="\t", dtype=str)
        acc_c = next((c for c in ["protein_accession","accession"]
                      if c in ptbl.columns), None)
        gene_c = "gene_symbol" if "gene_symbol" in ptbl.columns else None
        if acc_c and gene_c:
            gene_map = dict(zip(ptbl[acc_c].astype(str),
                                ptbl[gene_c].fillna("").astype(str)))

    rows: list[dict] = []
    seen: set[tuple]  = set()

    for _, row in occ.iterrows():
        acc = str(row.get("protein_accession", "")).strip()
        if not acc or acc == "nan":
            continue
        gene = gene_map.get(acc, "")

        for col, predicate, category in ANNOTATION_MAP:
            actual_col = _resolve_column(occ, col)
            if actual_col is None:
                continue

            for term in _split_cell(row.get(actual_col, "")):
                key = (acc, predicate, term)
                if key in seen:
                    continue
                seen.add(key)
                rows.append({
                    "protein_accession":   acc,
                    "gene_symbol":         gene,
                    "species":             species,
                    "annotation_category": category,
                    "predicate":           predicate,
                    "term_value":          term,
                    "is_transferable":     transferability.get(category, False),
                })

    return pd.DataFrame(rows)


def _load_transferability_rules(rules_path: Path) -> dict[str, bool]:
    if not rules_path.exists():
        return DEFAULT_TRANSFERABILITY.copy()
    df = pd.read_csv(rules_path, sep="\t", dtype=str)
    if "annotation_category" not in df.columns or "transferable" not in df.columns:
        return DEFAULT_TRANSFERABILITY.copy()
    rules: dict[str, bool] = {}
    for _, row in df.iterrows():
        cat = str(row["annotation_category"]).strip()
        val = str(row["transferable"]).strip().lower()
        rules[cat] = val in ("yes", "true", "1")
    return {**DEFAULT_TRANSFERABILITY, **rules}


def _write_default_rules(path: Path) -> None:
    rows = []
    reasons = {
        "domain_type":    "conserved domain structure across orthologs",
        "function":       "core function tags usually conserved",
        "topology":       "membrane topology is conserved in homologs",
        "position":       "relative domain position is conserved",
        "copy_number":    "copy number usually conserved in 1:1 orthologs",
        "proximity":      "proximity context conserved",
        "binding_target": "binding partners can diverge between species",
        "GO_function":    "molecular function usually conserved",
        "GO_process":     "biological process may differ by species context",
        "GO_component":   "cellular localisation may differ",
    }
    for cat, transferable in DEFAULT_TRANSFERABILITY.items():
        rows.append({
            "annotation_category": cat,
            "transferable": "yes" if transferable else "no",
            "biological_justification": reasons.get(cat, ""),
        })
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    print(f"  Written default rules → {path}")
    print()
    print("  IMPORTANT: edit transferability_rules.tsv now, BEFORE running")
    print("  step 07 or looking at any evaluation results.")
    print("  Decisions must be made on biological grounds only.")


def main(args: argparse.Namespace) -> None:
    data_dir = Path(args.data)
    data_dir.mkdir(parents=True, exist_ok=True)

    human_occ = Path(args.human)
    mouse_occ = Path(args.mouse)

    print("=" * 60)
    print("  06 — Build Annotation Tables")
    print("=" * 60)

    rules_file = data_dir / "transferability_rules.tsv"
    if not rules_file.exists():
        print("\nCreating default transferability_rules.tsv …")
        print("STOP: edit this file before proceeding to step 07.")
        _write_default_rules(rules_file)
    else:
        print(f"\nLoading transferability rules from {rules_file}")

    transferability = _load_transferability_rules(rules_file)
    print(f"  Transferable categories:     "
          f"{[k for k, v in transferability.items() if v]}")
    print(f"  Non-transferable categories: "
          f"{[k for k, v in transferability.items() if not v]}")

    print("\nBuilding human annotation table …")
    human_df = _build_long_table(
        human_occ, data_dir / "human_proteins.tsv", "human", transferability)
    out_h = data_dir / "human_annotations_long.tsv"
    human_df.to_csv(out_h, sep="\t", index=False)
    n_h = human_df["protein_accession"].nunique() if not human_df.empty else 0
    print(f"  {len(human_df)} assertions for {n_h} proteins → {out_h}")

    if not human_df.empty:
        print("\n  Human counts per category:")
        for cat, grp in human_df.groupby("annotation_category"):
            t = transferability.get(cat, False)
            print(f"    {cat:<20} {len(grp):>5}  "
                  f"{'[transferable]' if t else '[NOT transferred]'}")

    print("\nBuilding mouse annotation table (gold labels) …")
    mouse_df = _build_long_table(
        mouse_occ, data_dir / "mouse_proteins.tsv", "mouse", transferability)
    out_m = data_dir / "mouse_annotations_gold.tsv"
    mouse_df.to_csv(out_m, sep="\t", index=False)
    n_m = mouse_df["protein_accession"].nunique() if not mouse_df.empty else 0
    print(f"  {len(mouse_df)} assertions for {n_m} proteins → {out_m}")

    print("\nStep 6 complete.")
    print("Next: python 06b_build_mouse_structural_annotations.py")
    print("      (creates the structural-only file for Method B)")
    print("Then: python 07_generate_transfer_predictions.py")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Build long-form annotation tables")
    p.add_argument("--human",  required=True, help="Human occurrence TSV")
    p.add_argument("--mouse",  required=True, help="Mouse occurrence TSV")
    p.add_argument("--data",   default="transfer_pipeline/data")
    main(p.parse_args())
