"""
06b_build_mouse_structural_annotations.py
────────────────────────────────────────────────────────────────────────────
Step 6b — Extract ONLY structural annotations from the mouse table.

Purpose
───────
This file is what Method B in step 07 is allowed to see.
It contains only annotations that exist INDEPENDENTLY of what we are
trying to predict (function tags, binding targets, GO terms).

  ALLOWED (prediction-time structural features):
    domain_type       — which domains are present (domain_type_id)
    topology          — membrane topology context
    position          — N-terminal / C-terminal / Internal / etc.
    proximity         — proximity to signal peptide, TM helix, active site

  NOT INCLUDED (these are what we want to predict):
    function_tags     — DOF functional categories
    binding_target    — DOT binding partner categories
    go_terms_mf/bp/cc — Gene Ontology terms

By using this file in Method B, the predictor cannot peek at any of the
annotations it is trying to predict — enforcing the blind test protocol.

Output
──────
  mouse_annotations_structural.tsv
    protein_accession | gene_symbol | species |
    annotation_category | predicate | term_value | source_column
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd


STRUCTURAL_SPECS: dict[str, dict] = {
    "domain_type": {
        "predicate":  "hasDomainType",
        "candidates": ["domain_type_id", "domain_types", "interpro_id"],
    },
    "topology": {
        "predicate":  "hasTopologyContext",
        "candidates": ["topology_context", "topology"],
    },
    "position": {
        "predicate":  "hasPositionalCategory",
        "candidates": ["positional_category", "positional_categories",
                       "position_category"],
    },
    "proximity": {
        "predicate":  "hasProximityCategory",
        "candidates": ["proximity_categories", "proximity_category"],
    },
}


def _split_cell(val) -> list[str]:
    if pd.isna(val) or str(val).strip() in ("", "nan"):
        return []
    s = str(val).replace("||", ";").replace("|", ";").replace("//", ";")
    return [p.strip() for p in s.split(";")
            if p.strip() and p.strip().lower() != "nan"]


def _first_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def main(args: argparse.Namespace) -> None:
    out_dir  = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    occ_path = Path(args.mouse)
    out_file = out_dir / "mouse_annotations_structural.tsv"

    print("=" * 60)
    print("  06b — Build Mouse Structural Annotations (blind-safe)")
    print("=" * 60)

    if not occ_path.exists():
        raise FileNotFoundError(f"Input not found: {occ_path}")

    df = pd.read_csv(occ_path, sep="\t", dtype=str, low_memory=False)
    print(f"\nLoaded: {occ_path.name}  ({df.shape[0]} rows × {df.shape[1]} cols)")

    if "protein_accession" not in df.columns:
        raise ValueError("Input must contain 'protein_accession'")

    gene_col = next((c for c in ["batch_gene","gene_symbol","gene"]
                     if c in df.columns), None)

    used:    dict[str, str] = {}
    skipped: list[str]      = []

    for cat, spec in STRUCTURAL_SPECS.items():
        col = _first_col(df, spec["candidates"])
        if col is None:
            skipped.append(cat)
        else:
            used[cat] = col

    if "domain_type" not in used:
        raise ValueError(
            "No domain-type column found. "
            f"Tried: {STRUCTURAL_SPECS['domain_type']['candidates']}"
        )

    rows: list[dict] = []
    for _, row in df.iterrows():
        acc = str(row.get("protein_accession", "")).strip()
        if not acc or acc.lower() == "nan":
            continue
        gene = ""
        if gene_col:
            gene = str(row.get(gene_col, "")).strip()
            if gene.lower() == "nan":
                gene = ""
        for cat, src_col in used.items():
            pred   = STRUCTURAL_SPECS[cat]["predicate"]
            values = _split_cell(row.get(src_col, ""))
            for term in values:
                rows.append({
                    "protein_accession":   acc,
                    "gene_symbol":         gene,
                    "species":             "mouse",
                    "annotation_category": cat,
                    "predicate":           pred,
                    "term_value":          term,
                    "source_column":       src_col,
                })

    out_df = (pd.DataFrame(rows)
                .drop_duplicates()
                .sort_values(["protein_accession","annotation_category","term_value"])
                .reset_index(drop=True))

    if out_df.empty:
        raise ValueError("No structural annotations extracted. Check input.")

    out_df.to_csv(out_file, sep="\t", index=False)

    print("\nUsed source columns:")
    for cat, col in used.items():
        print(f"  {cat:<12} ← {col}")
    if skipped:
        print(f"\nSkipped (column absent): {skipped}")

    print(f"\nUnique proteins:   {out_df['protein_accession'].nunique()}")
    print(f"Total rows:        {len(out_df)}")
    print("\nRows per category:")
    for cat, n in out_df["annotation_category"].value_counts().sort_index().items():
        print(f"  {cat:<12} {n}")
    print(f"\nSaved → {out_file}")
    print("\nStep 6b complete.  Next: python 07_generate_transfer_predictions.py")


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Extract structural-only mouse annotations for blind Method B")
    p.add_argument("--mouse", required=True, help="Mouse occurrence TSV")
    p.add_argument("--out",   default="transfer_pipeline/data")
    main(p.parse_args())
