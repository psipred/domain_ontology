"""
04_build_ortholog_candidates.py
────────────────────────────────────────────────────────────────────────────
Step 4 — Score candidate human–mouse ortholog pairs.

Four independent signals per pair:
  1. Same MMseqs2 cluster        (weight 3.0 — strongest evidence)
  2. Gene-symbol similarity      (weight 2.5 × similarity score)
  3. Same protein group          (weight 1.0 — family-level)
  4. Similar sequence length     (weight 0.5 × length ratio)

Output
──────
  ortholog_candidates.tsv
    human_accession | human_gene | human_group |
    mouse_accession | mouse_gene | mouse_group |
    same_mmseqs_cluster | gene_symbol_exact | gene_symbol_similar |
    gene_similarity | same_protein_group |
    human_length_est | mouse_length_est | length_ratio |
    candidate_score | confidence_tier | curation_status

MANUAL STEP
───────────
  Open ortholog_candidates.tsv in Excel.
  Fill curation_status column:
    high_confidence   clear 1:1 ortholog
    ambiguous         paralog or unclear — review carefully
    excluded          not orthologs
  Save as ortholog_candidates_curated.tsv (or same filename).
  Then run: python 05_finalize_ortholog_pairs.py
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations
import argparse
from difflib import SequenceMatcher
from pathlib import Path
import pandas as pd


def _gene_sim(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    a, b = a.upper().strip(), b.upper().strip()
    if a == b:
        return 1.0
    a_base, b_base = a.rstrip("0123456789"), b.rstrip("0123456789")
    if a_base == b_base and a_base:
        return 0.9
    return round(SequenceMatcher(None, a, b).ratio(), 4)


def _get_seq_lengths(occ_path: Path) -> dict[str, int]:
    if not occ_path.exists():
        return {}
    try:
        occ = pd.read_csv(occ_path, sep="\t", dtype=str, low_memory=False,
                          usecols=lambda c: c in ["protein_accession","end"])
    except Exception:
        return {}
    if "protein_accession" not in occ.columns or "end" not in occ.columns:
        return {}
    lengths: dict[str, int] = {}
    for acc, grp in occ.groupby("protein_accession"):
        try:
            lengths[str(acc)] = int(grp["end"].dropna().astype(float).max())
        except (ValueError, TypeError):
            pass
    return lengths


def main(args: argparse.Namespace) -> None:
    data_dir = Path(args.data)
    data_dir.mkdir(parents=True, exist_ok=True)
    print("=" * 60)
    print("  04 — Build Ortholog Candidates")
    print("=" * 60)

    human_tsv   = data_dir / "human_proteins.tsv"
    mouse_tsv   = data_dir / "mouse_proteins.tsv"
    cluster_tsv = data_dir / "mmseqs2_clusters.tsv"

    if not human_tsv.exists() or not mouse_tsv.exists():
        print("ERROR: run 01_prepare_protein_tables.py first")
        return

    human = pd.read_csv(human_tsv, sep="\t", dtype=str).fillna("")
    mouse = pd.read_csv(mouse_tsv, sep="\t", dtype=str).fillna("")

    has_clusters = cluster_tsv.exists()
    cluster_map: dict[str, str] = {}
    if has_clusters:
        cdf   = pd.read_csv(cluster_tsv, sep="\t", dtype=str)
        acc_c = next((c for c in ["protein_accession","accession"]
                      if c in cdf.columns), None)
        cid_c = "cluster_id" if "cluster_id" in cdf.columns else None
        if acc_c and cid_c:
            cluster_map = dict(zip(cdf[acc_c].astype(str),
                                   cdf[cid_c].astype(str)))
        print(f"  Loaded {len(cluster_map)} cluster assignments")
    else:
        print("  No mmseqs2_clusters.tsv — cluster signal will be 0")
        print("  (run 03_parse_mmseqs2_results.py first for best results)")

    h_lens = _get_seq_lengths(Path(args.human)) if args.human else {}
    m_lens = _get_seq_lengths(Path(args.mouse)) if args.mouse else {}

    print(f"\nBuilding candidates: {len(human)} human × {len(mouse)} mouse …")
    rows = []

    for _, h in human.iterrows():
        h_acc  = str(h["protein_accession"])
        h_gene = str(h.get("gene_symbol", ""))
        h_grp  = str(h.get("protein_group", "other"))
        h_cid  = cluster_map.get(h_acc, "")
        h_len  = h_lens.get(h_acc, 0)

        for _, m in mouse.iterrows():
            m_acc  = str(m["protein_accession"])
            m_gene = str(m.get("gene_symbol", ""))
            m_grp  = str(m.get("protein_group", "other"))
            m_cid  = cluster_map.get(m_acc, "")
            m_len  = m_lens.get(m_acc, 0)

            same_cluster = (h_cid != "" and h_cid == m_cid)
            s_cluster    = 3.0 if same_cluster else 0.0

            gene_sim    = _gene_sim(h_gene, m_gene)
            same_symbol = gene_sim == 1.0
            sim_symbol  = gene_sim >= 0.8 and not same_symbol
            s_gene      = 2.5 * gene_sim

            same_grp = (h_grp == m_grp and h_grp not in ("other",""))
            s_group  = 1.0 if same_grp else 0.0

            if h_len > 0 and m_len > 0:
                lr    = min(h_len, m_len) / max(h_len, m_len)
                s_len = lr
            else:
                lr, s_len = 0.0, 0.0

            score = s_cluster + s_gene + s_group + s_len * 0.5
            if score < 1.0:
                continue

            if same_cluster and same_symbol:
                tier = "high_confidence"
            elif same_cluster and same_grp:
                tier = "high_confidence"
            elif same_symbol and same_grp:
                tier = "probable"
            elif same_cluster or same_symbol:
                tier = "probable"
            elif same_grp and gene_sim >= 0.7:
                tier = "possible"
            else:
                tier = "low_confidence"

            rows.append({
                "human_accession":     h_acc,
                "human_gene":          h_gene,
                "human_group":         h_grp,
                "mouse_accession":     m_acc,
                "mouse_gene":          m_gene,
                "mouse_group":         m_grp,
                "same_mmseqs_cluster": same_cluster,
                "mmseqs_cluster_id":   h_cid if same_cluster else "",
                "gene_symbol_exact":   same_symbol,
                "gene_symbol_similar": sim_symbol,
                "gene_similarity":     round(gene_sim, 4),   # preserved for step 08 difficulty
                "same_protein_group":  same_grp,
                "human_length_est":    h_len,
                "mouse_length_est":    m_len,
                "length_ratio":        round(lr, 3),
                "candidate_score":     round(score, 3),
                "confidence_tier":     tier,
                "curation_status":     "",   # fill manually
            })

    df = (pd.DataFrame(rows)
            .sort_values("candidate_score", ascending=False)
            .reset_index(drop=True))

    out_path = data_dir / "ortholog_candidates.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    print(f"\nSaved {len(df)} candidate pairs → {out_path}")
    print("\nTier distribution:")
    for tier, grp in df.groupby("confidence_tier"):
        print(f"  {tier:<20} {len(grp)} pairs")

    print(f"\n{'='*60}")
    print("  MANUAL STEP REQUIRED")
    print(f"{'='*60}")
    print(f"""
  1. Open: {out_path}
  2. Fill 'curation_status' for each pair:
       high_confidence  — clear 1:1 ortholog
       ambiguous        — paralog or unclear
       excluded         — not orthologs
  3. Save as (same file or):  {data_dir}/ortholog_candidates_curated.tsv
  4. Run:  python 05_finalize_ortholog_pairs.py
""")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Score ortholog candidate pairs")
    p.add_argument("--data",  default="transfer_pipeline/data")
    p.add_argument("--human", default=None, help="Human occurrence TSV (for length estimates)")
    p.add_argument("--mouse", default=None, help="Mouse occurrence TSV (for length estimates)")
    main(p.parse_args())
