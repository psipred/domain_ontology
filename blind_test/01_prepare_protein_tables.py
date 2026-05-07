"""
01_prepare_protein_tables.py
────────────────────────────────────────────────────────────────────────────
Step 1 — Build clean, normalised protein tables for human and mouse.

Inputs
──────
  --human   human cAMP merged domain occurrence TSV  (WITH_IDS version)
  --mouse   mouse cAMP merged domain occurrence TSV  (WITH_IDS version)
  --out     output directory (default: transfer_pipeline/data)

  Optional per-species group-label files (protein_accession | protein_group):
  --human_groups  (skip if absent — proteins assigned 'other')
  --mouse_groups

Outputs
───────
  <out>/human_proteins.tsv   protein_accession | gene_symbol | species |
  <out>/mouse_proteins.tsv     protein_group | family_module | n_domains
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd


FAMILY_TO_MODULE: dict[str, str] = {
    "GPCR_class_A": "receptor",       "GPCR_class_B": "receptor",
    "GPCR_class_C": "receptor",       "peptide_ligand": "receptor",
    "secreted_binding_protein": "receptor",
    "adenylyl_cyclase": "effector",   "phosphodiesterase": "effector",
    "cAMP_effector_POPDC": "effector",
    "AGC_kinase_PKA_like": "kinase",  "AGC_kinase_other": "kinase",
    "AKT_kinase": "kinase",           "CaMK": "kinase",
    "MAP2K": "kinase",                "MAPK": "kinase",
    "MAPK_related_kinase_or_review": "kinase",
    "MAPK_scaffold": "kinase",        "RAF_kinase": "kinase",
    "PAK_kinase": "kinase",           "ROCK_kinase": "kinase",
    "PI3K_catalytic": "kinase",       "PI3K_regulatory": "kinase",
    "heterotrimeric_G_alpha": "G_protein",
    "small_GTPase_Rac_family": "G_protein",
    "small_GTPase_other": "G_protein",
    "RapGEF_EPAC": "G_protein",       "Rho_GEF": "G_protein",
    "calcium_pump": "calcium",        "calmodulin": "calcium",
    "voltage_gated_calcium_channel": "calcium",
    "cyclic_nucleotide_gated_channel": "calcium",
    "ryanodine_receptor": "calcium",
    "ABC_channel": "channel",         "ABC_transporter": "channel",
    "HCN_channel": "channel",
    "protein_phosphatase_PP1": "phosphatase",
    "phospholipase_C": "lipid",       "phospholipase_D": "lipid",
    "transcription_factor": "transcription",
    "transcriptional_coactivator": "transcription",
}


def _build_table(occ_path: Path, groups_path: Path | None,
                 species: str) -> pd.DataFrame:
    # ── Hard-fail if input file is missing ───────────────────────────────────
    if not occ_path.exists():
        raise FileNotFoundError(
            f"\nInput file not found: {occ_path}\n"
            f"  Absolute path tried: {occ_path.resolve()}\n"
            f"  Tip: make sure you pass the correct path to --human / --mouse.\n"
            f"  Example (from the folder containing the TSV files):\n"
            f"    python 01_prepare_protein_tables.py \\\n"
            f"        --human c_AMP_pathway_human_MERGED_domain_occurrences_WITH_IDS.tsv \\\n"
            f"        --mouse c_AMP_pathway_mouse_MERGED_domain_occurrences_WITH_IDS.tsv\n"
            f"  Or use an absolute path:\n"
            f"    --human C:/path/to/c_AMP_pathway_human_MERGED_domain_occurrences_WITH_IDS.tsv"
        )

    occ = pd.read_csv(occ_path, sep="\t", dtype=str, low_memory=False)
    print(f"  Loaded {species}: {occ.shape[0]} rows × {occ.shape[1]} cols")
    print(f"  Columns present: {sorted(occ.columns.tolist())[:12]} ...")

    # ── Hard-fail if protein_accession column is missing ─────────────────────
    if "protein_accession" not in occ.columns:
        raise ValueError(
            f"\n'protein_accession' column not found in {occ_path.name}\n"
            f"  Columns available: {sorted(occ.columns.tolist())}\n"
            f"  Make sure you are pointing --{species} at the occurrence table\n"
            f"  (the one ending in _WITH_IDS.tsv), not a protein list file."
        )

    # Gene symbols — try common column names in priority order
    gene_map: dict[str, str] = {}
    for gcol in ["batch_gene", "gene_symbol", "gene"]:
        if gcol in occ.columns and "protein_accession" in occ.columns:
            for _, row in occ.iterrows():
                acc  = str(row.get("protein_accession", "")).strip()
                gene = str(row.get(gcol, "")).strip()
                if acc and gene and gene.lower() not in ("nan", ""):
                    gene_map.setdefault(acc, gene)
            break

    # Domain counts
    n_domains: dict[str, int] = {}
    if "protein_accession" in occ.columns:
        n_domains = occ.groupby("protein_accession").size().to_dict()

    proteins = sorted(occ["protein_accession"].dropna().unique()
                      if "protein_accession" in occ.columns else [])

    # Group labels
    group_map: dict[str, str] = {}
    if groups_path and groups_path.exists():
        gdf  = pd.read_csv(groups_path, sep="\t", dtype=str)
        acc_c = next((c for c in ["protein_accession","accession"] if c in gdf.columns), None)
        grp_c = next((c for c in ["protein_group","group"] if c in gdf.columns), None)
        if acc_c and grp_c:
            group_map = dict(zip(gdf[acc_c].astype(str), gdf[grp_c].astype(str)))

    rows = []
    for acc in proteins:
        grp    = group_map.get(acc, "other")
        module = FAMILY_TO_MODULE.get(grp, "other")
        rows.append({
            "protein_accession": acc,
            "gene_symbol":       gene_map.get(acc, ""),
            "species":           species,
            "protein_group":     grp,
            "family_module":     module,
            "n_domains":         n_domains.get(acc, 0),
        })

    df = pd.DataFrame(rows)
    df["protein_accession"] = df["protein_accession"].str.strip()
    df["species"]           = df["species"].str.lower()
    df = df.drop_duplicates(subset="protein_accession").reset_index(drop=True)
    return df


def main(args: argparse.Namespace) -> None:
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    print("=" * 60)
    print("  01 — Prepare Protein Tables")
    print("=" * 60)

    print("\nBuilding human protein table …")
    human_df = _build_table(
        Path(args.human),
        Path(args.human_groups) if args.human_groups else None,
        "human",
    )
    print(f"  {len(human_df)} proteins  |  "
          f"{human_df['protein_group'].nunique()} groups  |  "
          f"{human_df['family_module'].nunique()} modules")
    human_df.to_csv(out_dir / "human_proteins.tsv", sep="\t", index=False)
    print(f"  Saved → {out_dir}/human_proteins.tsv")

    print("\nBuilding mouse protein table …")
    mouse_df = _build_table(
        Path(args.mouse),
        Path(args.mouse_groups) if args.mouse_groups else None,
        "mouse",
    )
    print(f"  {len(mouse_df)} proteins  |  "
          f"{mouse_df['protein_group'].nunique()} groups  |  "
          f"{mouse_df['family_module'].nunique()} modules")
    mouse_df.to_csv(out_dir / "mouse_proteins.tsv", sep="\t", index=False)
    print(f"  Saved → {out_dir}/mouse_proteins.tsv")

    print(f"\nFamily module distribution:")
    combined = pd.concat([human_df, mouse_df])
    for (sp, mod), grp in combined.groupby(["species", "family_module"]):
        print(f"  {sp:<8} {mod:<18} {len(grp)} proteins")

    # ── Validate outputs are non-empty ────────────────────────────────────────
    if human_df.empty or mouse_df.empty:
        raise RuntimeError(
            "\nERROR: one or both protein tables are empty after processing.\n"
            f"  Human rows: {len(human_df)}  Mouse rows: {len(mouse_df)}\n"
            "  This usually means the path passed to --human / --mouse was wrong.\n"
            "  Check:\n"
            "    1. The paths you passed to --human and --mouse exist\n"
            "    2. The files contain a 'protein_accession' column\n"
            "    3. You are pointing at the occurrence table (_WITH_IDS.tsv),\n"
            "       not a domain-types or protein-list file."
        )

    print(f"\n{'='*50}")
    print(f"  Step 1 SUCCESS")
    print(f"  Human: {len(human_df)} proteins → {out_dir}/human_proteins.tsv")
    print(f"  Mouse: {len(mouse_df)} proteins → {out_dir}/mouse_proteins.tsv")
    print(f"{'='*50}")
    print("  Next: python 02_fetch_sequences.py")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Build protein tables")
    p.add_argument("--human",        required=True, help="Human occurrence TSV")
    p.add_argument("--mouse",        required=True, help="Mouse occurrence TSV")
    p.add_argument("--out",          default="transfer_pipeline/data")
    p.add_argument("--human_groups", default=None)
    p.add_argument("--mouse_groups", default=None)
    main(p.parse_args())
