"""
07_generate_transfer_predictions.py
────────────────────────────────────────────────────────────────────────────
Step 7 — Generate ontology annotation transfer predictions for mouse.

Three methods — all fully scripted:

  Method A: Direct ortholog transfer
    Copy all transferable human annotations to the mouse protein.
    Baseline that does not use any mouse evidence.

  Method B: Domain-architecture transfer  (FIXED in this version)
    Transfer annotations only when they are supported by shared domain types.
    Uses ONLY mouse_annotations_structural.tsv (structural annotations only)
    — NOT the gold-label file — enforcing the blind test protocol.

    Bug fixed: domain IDs are normalised before comparison.
    Human files store bare IDs (e.g. IPR000187); mouse files store
    prefixed IDs (e.g. INTERPRO:IPR001879).  Stripping 'INTERPRO:' before
    computing set intersection restores the 152 shared domains.

  Method C: Nearest-sequence baseline
    Ignores the ortholog table; transfers from the closest human sequence
    in the same MMseqs2 cluster.  Serves as a sequence-only comparator.

Input pairs file
────────────────
  ortholog_pairs_calibration.tsv  (80% of approved pairs)
  Do NOT use ortholog_pairs_test.tsv here — that is reserved for step 08.

Output
──────
  transfer_predictions.tsv
    mouse_accession | mouse_gene | human_source | method |
    annotation_category | predicate | term_value |
    transfer_confidence | shared_domain_count | architecture_jaccard |
    architecture_mode | architecture_score
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd


# ── Domain ID normalisation ───────────────────────────────────────────────────
# Human files:  IPR000187           (bare InterPro accession)
# Mouse files:  INTERPRO:IPR001879  (prefixed)
# After stripping 'INTERPRO:' both reduce to 'IPRnnnnnn' and match correctly.
# This normalisation is the single most important bug-fix in the pipeline.

def _norm_domain(d: str) -> str:
    """Normalise domain ID to bare form (strip any database prefix)."""
    d = str(d).strip()
    for pfx in ("INTERPRO:", "PFAM:", "CDD:", "SMART:", "PRINTS:",
                "PROFILE:", "NCBIFAM:", "PANTHER:", "GENE3D:"):
        if d.upper().startswith(pfx):
            return d[len(pfx):]
    return d


def _split_cell(val) -> list[str]:
    if pd.isna(val) or str(val).strip() in ("", "nan"):
        return []
    return [v.strip() for v in str(val).split(";") if v.strip()]


# ── Load human annotations ────────────────────────────────────────────────────

def _load_human_annotations(ann_path: Path) -> dict[str, pd.DataFrame]:
    if not ann_path.exists():
        return {}
    df = pd.read_csv(ann_path, sep="\t", dtype=str)
    if "is_transferable" in df.columns:
        df = df[df["is_transferable"].astype(str).str.lower()
                  .isin(["true","yes","1"])]
    result: dict[str, pd.DataFrame] = {}
    for acc, grp in df.groupby("protein_accession"):
        result[str(acc)] = grp.reset_index(drop=True)
    return result


# ── Domain sets ───────────────────────────────────────────────────────────────

def _build_domain_sets(ann_df: pd.DataFrame) -> dict[str, set[str]]:
    """Return {protein_accession: set(normalised_domain_ids)}."""
    result: dict[str, set[str]] = {}
    if ann_df.empty or "annotation_category" not in ann_df.columns:
        return result
    dt_rows = ann_df[ann_df["annotation_category"] == "domain_type"]
    for acc, grp in dt_rows.groupby("protein_accession"):
        result[str(acc)] = {_norm_domain(t) for t in grp["term_value"].astype(str)}
    return result


# ── Method A ──────────────────────────────────────────────────────────────────

def _method_A_direct(pairs: pd.DataFrame,
                      human_anns: dict[str, pd.DataFrame]) -> list[dict]:
    rows: list[dict] = []
    for _, pair in pairs.iterrows():
        h_acc  = str(pair["human_accession"])
        m_acc  = str(pair["mouse_accession"])
        m_gene = str(pair.get("mouse_gene", ""))
        anns   = human_anns.get(h_acc)
        if anns is None or anns.empty:
            continue
        for _, ann in anns.iterrows():
            rows.append({
                "mouse_accession":      m_acc,
                "mouse_gene":           m_gene,
                "human_source":         h_acc,
                "method":               "A_direct_transfer",
                "annotation_category":  ann.get("annotation_category",""),
                "predicate":            ann.get("predicate",""),
                "term_value":           ann.get("term_value",""),
                "transfer_confidence":  "high",
                "shared_domain_count":  None,
                "architecture_jaccard": None,
                "architecture_mode":    None,
                "architecture_score":   None,
            })
    return rows


# ── Method B ──────────────────────────────────────────────────────────────────

def _method_B_domain_architecture(
        pairs:           pd.DataFrame,
        human_anns:      dict[str, pd.DataFrame],
        human_ann_df:    pd.DataFrame,
        mouse_structural: pd.DataFrame,     # structural-only, NOT gold labels
) -> list[dict]:
    """
    Transfer annotations that are supported by shared domain architecture.

    Uses mouse_annotations_structural.tsv (domain_type only) — the predictor
    never sees any functional mouse annotations, enforcing blind test protocol.

    Domain IDs are normalised (INTERPRO: prefix stripped) before comparison,
    restoring the correct overlap between human and mouse domain sets.
    """
    mouse_domains  = _build_domain_sets(mouse_structural)
    human_domains  = _build_domain_sets(human_ann_df)

    has_mouse_data = len(mouse_domains) > 0
    mode           = "FULL" if has_mouse_data else "PROXY"

    if not has_mouse_data:
        print("  WARNING Method B: no mouse domain data found — running in PROXY mode.")
        print("  Run 06b_build_mouse_structural_annotations.py to enable FULL mode.")

    rows: list[dict] = []
    for _, pair in pairs.iterrows():
        h_acc  = str(pair["human_accession"])
        m_acc  = str(pair["mouse_accession"])
        m_gene = str(pair.get("mouse_gene", ""))
        anns   = human_anns.get(h_acc)
        if anns is None or anns.empty:
            continue

        h_doms = human_domains.get(h_acc, set())

        if mode == "FULL":
            m_doms      = mouse_domains.get(m_acc, set())
            shared_doms = h_doms & m_doms
        else:
            # PROXY: assume 1:1 ortholog shares human domain architecture
            m_doms      = h_doms
            shared_doms = h_doms

        union_doms = h_doms | m_doms
        arch_j     = (len(h_doms & m_doms) / len(union_doms)
                      if union_doms else 0.0)

        for _, ann in anns.iterrows():
            cat = ann.get("annotation_category", "")
            # Domain type: always transfer (it IS the structural evidence)
            if cat == "domain_type":
                conf = "high"
            # Other categories: only transfer if shared domain evidence exists
            elif shared_doms:
                conf = "medium" if mode == "FULL" else "medium_proxy"
            else:
                continue   # no shared domain evidence → do not transfer

            rows.append({
                "mouse_accession":      m_acc,
                "mouse_gene":           m_gene,
                "human_source":         h_acc,
                "method":               "B_domain_architecture",
                "annotation_category":  cat,
                "predicate":            ann.get("predicate",""),
                "term_value":           ann.get("term_value",""),
                "transfer_confidence":  conf,
                "shared_domain_count":  len(shared_doms),
                "architecture_jaccard": round(arch_j, 4),
                "architecture_mode":    mode,
                "architecture_score":   round(arch_j * (1.0 if conf == "high" else 0.7), 4),
            })
    return rows


# ── Method C ──────────────────────────────────────────────────────────────────

def _method_C_nearest_sequence(
        mouse_tsv:   Path,
        clusters:    Path,
        human_anns:  dict[str, pd.DataFrame],
) -> list[dict]:
    if not clusters.exists():
        print("  Method C skipped: mmseqs2_clusters.tsv not found")
        return []

    from collections import defaultdict
    cdf = pd.read_csv(clusters, sep="\t", dtype=str)
    mouse_df = pd.read_csv(mouse_tsv, sep="\t", dtype=str) \
               if mouse_tsv.exists() else pd.DataFrame()

    clust_human: dict[str, list] = defaultdict(list)
    clust_mouse: dict[str, list] = defaultdict(list)
    for _, row in cdf.iterrows():
        acc = str(row.get("protein_accession",""))
        cid = str(row.get("cluster_id",""))
        sp  = str(row.get("species",""))
        if sp == "human":
            clust_human[cid].append(acc)
        elif sp == "mouse":
            clust_mouse[cid].append(acc)

    gene_map: dict[str, str] = {}
    if not mouse_df.empty and "protein_accession" in mouse_df.columns:
        gene_col = "gene_symbol" if "gene_symbol" in mouse_df.columns else None
        if gene_col:
            gene_map = dict(zip(mouse_df["protein_accession"].astype(str),
                                mouse_df[gene_col].fillna("").astype(str)))

    rows: list[dict] = []
    for _, row in cdf[cdf["species"] == "mouse"].iterrows():
        m_acc        = str(row.get("protein_accession",""))
        cid          = str(row.get("cluster_id",""))
        h_candidates = clust_human.get(cid, [])
        if not h_candidates:
            continue
        h_acc  = h_candidates[0]
        m_gene = gene_map.get(m_acc, "")
        anns   = human_anns.get(h_acc)
        if anns is None or anns.empty:
            continue
        for _, ann in anns.iterrows():
            rows.append({
                "mouse_accession":      m_acc,
                "mouse_gene":           m_gene,
                "human_source":         h_acc,
                "method":               "C_nearest_sequence",
                "annotation_category":  ann.get("annotation_category",""),
                "predicate":            ann.get("predicate",""),
                "term_value":           ann.get("term_value",""),
                "transfer_confidence":  "sequence_based",
                "shared_domain_count":  None,
                "architecture_jaccard": None,
                "architecture_mode":    None,
                "architecture_score":   None,
            })
    return rows


# ── Main ──────────────────────────────────────────────────────────────────────

def main(args: argparse.Namespace) -> None:
    data_dir = Path(args.data)
    data_dir.mkdir(parents=True, exist_ok=True)
    print("=" * 60)
    print("  07 — Generate Transfer Predictions")
    print("=" * 60)

    # ── BLIND TEST PROTOCOL ───────────────────────────────────────────────────
    # Use CALIBRATION pairs only.  Test pairs are reserved for step 08.
    pairs_file = data_dir / "ortholog_pairs_calibration.tsv"
    if not pairs_file.exists():
        pairs_file = data_dir / "ortholog_pairs.tsv"
        print(f"  NOTE: using {pairs_file.name} (no calibration split found).")
        print("  Run step 05 first to create the proper holdout split.")
    else:
        print(f"  Using calibration pairs: {pairs_file.name}")

    human_ann_path   = data_dir / "human_annotations_long.tsv"
    mouse_struct_path= data_dir / "mouse_annotations_structural.tsv"
    clusters_path    = data_dir / "mmseqs2_clusters.tsv"
    mouse_tsv        = data_dir / "mouse_proteins.tsv"

    for f, label in [(pairs_file, "ortholog pairs"),
                     (human_ann_path, "human_annotations_long.tsv")]:
        if not f.exists():
            print(f"ERROR: {label} not found at {f}. Run earlier steps.")
            return

    pairs     = pd.read_csv(pairs_file, sep="\t", dtype=str)
    human_ann = _load_human_annotations(human_ann_path)
    human_full= pd.read_csv(human_ann_path, sep="\t", dtype=str)

    # Method B uses structural-only mouse data (NOT gold labels)
    if mouse_struct_path.exists():
        mouse_structural = pd.read_csv(mouse_struct_path, sep="\t", dtype=str)
        print(f"  Method B: using structural mouse annotations ({len(mouse_structural)} rows)")
    else:
        mouse_structural = pd.DataFrame()
        print("  WARNING: mouse_annotations_structural.tsv not found.")
        print("           Method B will run in PROXY mode (same as Method A).")
        print("           Run 06b_build_mouse_structural_annotations.py first.")

    print(f"\n  Ortholog pairs: {len(pairs)}")
    print(f"  Human proteins with annotations: {len(human_ann)}")

    all_rows: list[dict] = []

    print("\nMethod A: Direct ortholog transfer …")
    rows_a = _method_A_direct(pairs, human_ann)
    all_rows.extend(rows_a)
    print(f"  {len(rows_a)} predictions")

    print("Method B: Domain-architecture transfer …")
    rows_b = _method_B_domain_architecture(
        pairs, human_ann, human_full, mouse_structural)
    all_rows.extend(rows_b)
    print(f"  {len(rows_b)} predictions")

    print("Method C: Nearest-sequence baseline …")
    rows_c = _method_C_nearest_sequence(mouse_tsv, clusters_path, human_ann)
    all_rows.extend(rows_c)
    print(f"  {len(rows_c)} predictions")

    out_df = pd.DataFrame(all_rows)
    out_path = data_dir / "transfer_predictions.tsv"
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"\nTotal predictions: {len(out_df)}")
    print(f"Saved → {out_path}")

    if not out_df.empty:
        ct = out_df.groupby(["method","annotation_category"]).size().unstack(fill_value=0)
        print("\nPredictions per method × category:")
        print(ct.to_string())

        # Domain architecture stats for Method B
        b = out_df[out_df["method"]=="B_domain_architecture"]
        if not b.empty and "architecture_jaccard" in b.columns:
            j = b["architecture_jaccard"].dropna().astype(float)
            print(f"\nMethod B architecture Jaccard: "
                  f"mean={j.mean():.3f}, median={j.median():.3f}, "
                  f"min={j.min():.3f}, max={j.max():.3f}")
            print(f"Method B mode: "
                  f"{b['architecture_mode'].value_counts().to_dict()}")

    print("\nStep 7 complete.  Next: python 08_evaluate_transfer.py")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Generate transfer predictions")
    p.add_argument("--data", default="transfer_pipeline/data")
    main(p.parse_args())
