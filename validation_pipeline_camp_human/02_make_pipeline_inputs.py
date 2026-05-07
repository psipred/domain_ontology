#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
make_pipeline_inputs.py

Create pipeline input TSV files from your merged domain-occurrence table.

This version uses hardcoded Windows paths, so you can just edit the
PATH SETTINGS section and run the script directly in PyCharm.

Outputs:
    foreground_protein_domains.tsv
    protein_group_labels.tsv
    protein_group_labels_helper.tsv
    gold_standard_relationships.tsv
    ontology_vs_external_mapping.tsv

Optional:
    background_protein_domains.tsv
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, List, Tuple

import pandas as pd


# =============================================================================
# PATH SETTINGS
# =============================================================================
DOMAIN_TABLE_PATH = r"C:\Users\32045\Desktop\测试\integrated_pipeline _calcium_human\input\calcium_signaling_pathway_human_MERGED_domain_occurrences_WITH_IDS.tsv"

OUTPUT_DIR = r"C:\Users\32045\Desktop\测试\integrated_pipeline _calcium_human\input"

BACKGROUND_SOURCE_PATH = r"C:\Users\32045\Desktop\测试\integrated_pipeline _calcium_human\input\background_protein_domains.tsv"


# =============================================================================
# REGEX PATTERNS
# =============================================================================
GO_RE = re.compile(r"GO:\d{7}")
REACTOME_RE = re.compile(r"R-HSA-\d+")
KEGG_RE = re.compile(r"(?:hsa:\d+|hsa\d+)")
DOF_RE = re.compile(r"DOF:\d+")
DOT_RE = re.compile(r"DOT:\d+")
DOP_RE = re.compile(r"DOP:\d+")
IPR_RE = re.compile(r"IPR\d+")


# =============================================================================
# HELPERS
# =============================================================================
def ordered_unique(values: Iterable[str]) -> List[str]:
    seen = set()
    out = []
    for v in values:
        if v is None:
            continue
        s = str(v).strip()
        if not s or s.lower() == "nan":
            continue
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


def join_semicolon(values: Iterable[str]) -> str:
    return ";".join(ordered_unique(values))


def split_semicolon(cell: str) -> List[str]:
    if cell is None:
        return []
    s = str(cell).strip()
    if not s or s.lower() == "nan":
        return []
    return [x.strip() for x in s.split(";") if x and str(x).strip()]


def extract_regex_ids(text: str, regex: re.Pattern) -> List[str]:
    if text is None:
        return []
    return ordered_unique(regex.findall(str(text)))


def extract_ontology_ids_from_semicolon_field(
    cell: str,
    regex: re.Pattern,
    skip_unknown: bool = False
) -> List[str]:
    out = []
    for part in split_semicolon(cell):
        low = part.lower()
        if skip_unknown and "unknown" in low:
            continue
        found = regex.findall(part)
        out.extend(found)
    return ordered_unique(out)


def require_columns(df: pd.DataFrame, required: List[str]) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            "Missing required columns:\n"
            f"{missing}\n\n"
            f"Found columns are:\n{list(df.columns)}"
        )


def normalize_interpro_id(value: str) -> str:
    """
    Normalize:
      IPR017452 -> INTERPRO:IPR017452
      INTERPRO:IPR017452 -> INTERPRO:IPR017452
    """
    if value is None:
        return ""
    s = str(value).strip()
    if not s or s.lower() == "nan":
        return ""
    if s.startswith("INTERPRO:"):
        return s
    if re.fullmatch(r"IPR\d+", s):
        return f"INTERPRO:{s}"
    return s


def detect_background_columns(df: pd.DataFrame) -> Tuple[str, str]:
    protein_candidates = [
        "protein_accession", "accession", "Entry", "entry",
        "uniprot_accession", "UniProtKB-AC", "Protein"
    ]
    domain_candidates = [
        "domain_type_id", "interpro_id", "InterPro", "interpro",
        "interpro_accession", "InterPro accession", "ipr", "IPR"
    ]

    protein_col = next((c for c in protein_candidates if c in df.columns), None)
    domain_col = next((c for c in domain_candidates if c in df.columns), None)

    if protein_col is None or domain_col is None:
        raise ValueError(
            "Could not auto-detect background columns.\n"
            f"Found columns: {list(df.columns)}\n"
            "Need one protein column and one domain column."
        )
    return protein_col, domain_col


# =============================================================================
# BUILDERS
# =============================================================================
def build_foreground(df: pd.DataFrame, outdir: Path) -> None:
    fg = df[["protein_accession", "domain_type_id"]].copy()

    fg["protein_accession"] = fg["protein_accession"].astype(str).str.strip()
    fg["domain_type_id"] = fg["domain_type_id"].astype(str).str.strip()

    fg = fg[
        (fg["protein_accession"] != "") &
        (fg["domain_type_id"] != "") &
        (fg["protein_accession"].str.lower() != "nan") &
        (fg["domain_type_id"].str.lower() != "nan")
    ].drop_duplicates()

    fg = fg.sort_values(["protein_accession", "domain_type_id"])
    fg.to_csv(outdir / "foreground_protein_domains.tsv", sep="\t", index=False)


def build_protein_group_labels(df: pd.DataFrame, outdir: Path) -> None:
    proteins = (
        df[["protein_accession"]]
        .dropna()
        .drop_duplicates()
        .sort_values("protein_accession")
        .copy()
    )
    proteins["protein_accession"] = proteins["protein_accession"].astype(str).str.strip()
    proteins = proteins[
        (proteins["protein_accession"] != "") &
        (proteins["protein_accession"].str.lower() != "nan")
    ]

    proteins["protein_group"] = ""
    proteins.to_csv(outdir / "protein_group_labels.tsv", sep="\t", index=False)

    # Helper file to make manual filling easier
    if "batch_gene" in df.columns:
        helper = (
            df[["protein_accession", "batch_gene"]]
            .drop_duplicates()
            .sort_values(["protein_accession", "batch_gene"])
            .copy()
        )
        helper["protein_accession"] = helper["protein_accession"].fillna("").astype(str).str.strip()
        helper["batch_gene"] = helper["batch_gene"].fillna("").astype(str).str.strip()
    else:
        helper = proteins.copy()
        helper["batch_gene"] = ""

    helper.to_csv(outdir / "protein_group_labels_helper.tsv", sep="\t", index=False)


def build_external_mapping(df: pd.DataFrame, outdir: Path) -> None:
    records = []

    for protein, sub in df.groupby("protein_accession", dropna=True):
        protein = str(protein).strip()
        if not protein or protein.lower() == "nan":
            continue

        go_ids = []
        reactome_ids = []
        kegg_ids = []

        for _, row in sub.iterrows():
            if "go_terms_all" in row.index:
                go_ids.extend(extract_regex_ids(row.get("go_terms_all", ""), GO_RE))
            if "protein_go_cc_ids" in row.index:
                go_ids.extend(extract_regex_ids(row.get("protein_go_cc_ids", ""), GO_RE))
            if "protein_reactome_ids" in row.index:
                reactome_ids.extend(extract_regex_ids(row.get("protein_reactome_ids", ""), REACTOME_RE))
            if "protein_kegg_ids" in row.index:
                kegg_ids.extend(extract_regex_ids(row.get("protein_kegg_ids", ""), KEGG_RE))

        records.append({
            "protein_accession": protein,
            "go_ids": join_semicolon(go_ids),
            "reactome_ids": join_semicolon(reactome_ids),
            "kegg_ids": join_semicolon(kegg_ids),
        })

    out = pd.DataFrame(records).sort_values("protein_accession")
    out.to_csv(outdir / "ontology_vs_external_mapping.tsv", sep="\t", index=False)


def build_gold_standard(df: pd.DataFrame, outdir: Path) -> None:
    """
    Auto-seed a benchmark file from the ontology IDs already present in the table.
    Review this file manually before using it as a final gold-standard benchmark.
    """
    rows = []

    for domain_type, sub in df.groupby("domain_type_id", dropna=True):
        domain_type = str(domain_type).strip()
        if not domain_type or domain_type.lower() == "nan":
            continue

        function_ids = []
        binding_ids = []
        topology_ids = []
        positional_ids = []
        proximity_ids = []
        copy_number_ids = []

        for _, row in sub.iterrows():
            function_ids.extend(
                extract_ontology_ids_from_semicolon_field(
                    row.get("function_tags", ""), DOF_RE, skip_unknown=True
                )
            )

            # Keep your current spelling exactly: binding_parnter
            binding_ids.extend(
                extract_ontology_ids_from_semicolon_field(
                    row.get("binding_parnter", row.get("binding_partner", "")), DOT_RE, skip_unknown=True
                )
            )

            topology_ids.extend(
                extract_ontology_ids_from_semicolon_field(
                    row.get("topology_context", ""), DOP_RE, skip_unknown=True
                )
            )

            positional_ids.extend(
                extract_ontology_ids_from_semicolon_field(
                    row.get("positional_category", ""), DOP_RE, skip_unknown=True
                )
            )

            proximity_ids.extend(
                extract_ontology_ids_from_semicolon_field(
                    row.get("proximity_categories", ""), DOP_RE, skip_unknown=True
                )
            )

            copy_number_ids.extend(
                extract_ontology_ids_from_semicolon_field(
                    row.get("copy_number_property", ""), DOP_RE, skip_unknown=True
                )
            )

        for val in ordered_unique(function_ids):
            rows.append({
                "domain_type": domain_type,
                "relationship_type": "hasFunctionTag",
                "expected_value": val,
            })

        for val in ordered_unique(binding_ids):
            rows.append({
                "domain_type": domain_type,
                "relationship_type": "hasBindingTargetCategory",
                "expected_value": val,
            })

        for val in ordered_unique(topology_ids):
            rows.append({
                "domain_type": domain_type,
                "relationship_type": "hasTopologyContext",
                "expected_value": val,
            })

        for val in ordered_unique(positional_ids):
            rows.append({
                "domain_type": domain_type,
                "relationship_type": "hasPositionalCategory",
                "expected_value": val,
            })

        for val in ordered_unique(proximity_ids):
            rows.append({
                "domain_type": domain_type,
                "relationship_type": "hasProximityCategory",
                "expected_value": val,
            })

        for val in ordered_unique(copy_number_ids):
            rows.append({
                "domain_type": domain_type,
                "relationship_type": "hasCopyNumberProperty",
                "expected_value": val,
            })

    out = pd.DataFrame(rows).drop_duplicates()

    if len(out) == 0:
        out = pd.DataFrame(columns=["domain_type", "relationship_type", "expected_value"])
    else:
        out = out.sort_values(["domain_type", "relationship_type", "expected_value"])

    out.to_csv(outdir / "gold_standard_relationships.tsv", sep="\t", index=False)


def build_background(background_source: Path, outdir: Path) -> None:
    bg = pd.read_csv(background_source, sep="\t", dtype=str).fillna("")

    protein_col, domain_col = detect_background_columns(bg)

    out = bg[[protein_col, domain_col]].copy()
    out.columns = ["protein_accession", "domain_type_id"]

    out["protein_accession"] = out["protein_accession"].astype(str).str.strip()
    out["domain_type_id"] = out["domain_type_id"].astype(str).map(normalize_interpro_id)

    out = out[
        (out["protein_accession"] != "") &
        (out["domain_type_id"] != "") &
        (out["protein_accession"].str.lower() != "nan") &
        (out["domain_type_id"].str.lower() != "nan")
    ].drop_duplicates()

    out = out.sort_values(["protein_accession", "domain_type_id"])
    out.to_csv(outdir / "background_protein_domains.tsv", sep="\t", index=False)


# =============================================================================
# MAIN
# =============================================================================
def main() -> None:
    domain_path = Path(DOMAIN_TABLE_PATH)
    outdir = Path(OUTPUT_DIR)
    background_path = Path(BACKGROUND_SOURCE_PATH) if str(BACKGROUND_SOURCE_PATH).strip() else None

    outdir.mkdir(parents=True, exist_ok=True)

    if not domain_path.exists():
        raise FileNotFoundError(f"Domain table not found:\n{domain_path}")

    print(f"[INFO] Reading domain table:\n{domain_path}")
    # Try utf-8-sig first (handles BOM), fall back to latin-1 for Chinese Windows files
    try:
        df = pd.read_csv(domain_path, sep="\t", dtype=str, encoding="utf-8-sig").fillna("")
    except UnicodeDecodeError:
        df = pd.read_csv(domain_path, sep="\t", dtype=str, encoding="latin-1").fillna("")

    # Auto-detect column name variants so the script works with both
    # old format (binding_parnter, copy_number_property) and
    # new format (binding_partner, copy_number)
    col_aliases = {
        "binding_parnter":      ["binding_parnter", "binding_partner"],
        "copy_number_property": ["copy_number_property", "copy_number"],
        "protein_reactome_ids": ["protein_reactome_ids", "reactome_ids"],
        "protein_kegg_ids":     ["protein_kegg_ids", "kegg_ids"],
    }
    for canonical, candidates in col_aliases.items():
        if canonical not in df.columns:
            for alt in candidates:
                if alt in df.columns:
                    df[canonical] = df[alt]
                    break
            else:
                df[canonical] = ""   # column absent — fill with empty string

    required_columns = [
        "protein_accession",
        "domain_type_id",
        "function_tags",
        "topology_context",
        "positional_category",
        "proximity_categories",
        "go_terms_all",
        "protein_go_cc_ids",
    ]
    require_columns(df, required_columns)

    build_foreground(df, outdir)
    print("[OK] foreground_protein_domains.tsv created")

    build_protein_group_labels(df, outdir)
    print("[OK] protein_group_labels.tsv created")
    print("[OK] protein_group_labels_helper.tsv created")

    build_external_mapping(df, outdir)
    print("[OK] ontology_vs_external_mapping.tsv created")

    build_gold_standard(df, outdir)
    print("[OK] gold_standard_relationships.tsv created")

    if background_path:
        if not background_path.exists():
            raise FileNotFoundError(f"Background source not found:\n{background_path}")
        build_background(background_path, outdir)
        print("[OK] background_protein_domains.tsv created")
    else:
        print("[INFO] No background source path set, so background_protein_domains.tsv was not created.")

    print(f"\n[OK] Files written to:\n{outdir.resolve()}")
    print("Generated:")
    print("  - foreground_protein_domains.tsv")
    print("  - protein_group_labels.tsv")
    print("  - protein_group_labels_helper.tsv")
    print("  - gold_standard_relationships.tsv")
    print("  - ontology_vs_external_mapping.tsv")
    if background_path:
        print("  - background_protein_domains.tsv")


if __name__ == "__main__":
    main()