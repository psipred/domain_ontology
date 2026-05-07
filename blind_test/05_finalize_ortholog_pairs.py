"""
05_finalize_ortholog_pairs.py
────────────────────────────────────────────────────────────────────────────
Step 5 — Finalize ortholog pairs and create the blind-test holdout partition.

CRITICAL for a valid blind evaluation
──────────────────────────────────────
  This step creates two separate pair tables that are used differently by
  all downstream scripts:

  ortholog_pairs_calibration.tsv  (80%) — used by:
    • 06_build_annotation_tables.py   to write transferability_rules.tsv
    • 07_generate_transfer_predictions.py  for prediction generation
    • 08_evaluate_transfer.py  for threshold / calibration decisions

  ortholog_pairs_test.tsv          (20%, LOCKED) — used by:
    • 08_evaluate_transfer.py  ONLY, to compute final reported metrics

  Rule: do NOT look at or use the test set file until step 8.
  Tuning transferability_rules.tsv, method thresholds, or any other
  parameter on test-set results invalidates the blind evaluation.

Stratification
──────────────
  The 80/20 split is stratified by confidence_tier so both splits contain
  proportional representations of high_confidence, probable, and possible
  pairs.

Outputs
───────
  ortholog_pairs.tsv              all approved pairs (backward-compatible)
  ortholog_pairs_calibration.tsv  80% subset for development
  ortholog_pairs_test.tsv         20% holdout — BLIND TEST SET
  split_summary.txt               human-readable split report
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import numpy as np


def _stratified_split(df: pd.DataFrame,
                       test_frac: float = 0.20,
                       seed: int = 42) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Stratified split by confidence_tier.
    Ensures both calibration and test sets contain proportional tier counts.
    Falls back to random split when a tier has < 2 pairs.
    """
    rng   = np.random.default_rng(seed)
    cal_idx: list[int]  = []
    test_idx: list[int] = []

    for tier, grp in df.groupby("confidence_tier"):
        idx = grp.index.tolist()
        rng.shuffle(idx)
        n_test = max(1, round(len(idx) * test_frac))
        # Guarantee at least 1 pair in test only if tier has ≥ 2 pairs
        if len(idx) < 2:
            cal_idx.extend(idx)
        else:
            test_idx.extend(idx[:n_test])
            cal_idx.extend(idx[n_test:])

    return df.loc[cal_idx], df.loc[test_idx]


def main(args: argparse.Namespace) -> None:
    data_dir = Path(args.data)
    data_dir.mkdir(parents=True, exist_ok=True)
    print("=" * 60)
    print("  05 — Finalize Ortholog Pairs  (with blind holdout split)")
    print("=" * 60)

    # Accept either curated or raw candidates
    curated_path = data_dir / "ortholog_candidates_curated.tsv"
    raw_path     = data_dir / "ortholog_candidates.tsv"

    if curated_path.exists():
        src = curated_path
    elif raw_path.exists():
        src = raw_path
        print(f"\n  NOTE: using {raw_path.name} (no curated file found).")
        print("  For best results, complete manual curation first.")
    else:
        print(f"\nERROR: {raw_path} not found. Run step 04 first.")
        return

    df = pd.read_csv(src, sep="\t", dtype=str).fillna("")
    print(f"\nLoaded {len(df)} candidate pairs from {src.name}")

    # Filter to approved pairs
    if "curation_status" in df.columns:
        approved_vals = {"high_confidence","high confidence","yes","approved"}
        approved = df[df["curation_status"].str.lower().str.strip()
                        .isin(approved_vals)].copy()
        if approved.empty:
            print(f"  No 'high_confidence' rows. "
                  f"Values found: {df['curation_status'].value_counts().to_dict()}")
            print("  Falling back to auto high_confidence tier pairs.")
            approved = df[df["confidence_tier"].str.lower()
                            == "high_confidence"].copy()
    else:
        print("  No curation_status column — using all pairs above score threshold.")
        approved = df[df["candidate_score"].astype(float) >= 2.0].copy()

    if approved.empty:
        print("ERROR: no approved pairs found.")
        return

    print(f"  Approved pairs: {len(approved)}")

    # Ortholog type
    h_counts = approved["human_accession"].value_counts()
    m_counts = approved["mouse_accession"].value_counts()
    approved["ortholog_type"] = "one2one"
    approved.loc[approved["human_accession"].isin(
        h_counts[h_counts > 1].index), "ortholog_type"] = "one2many"
    approved.loc[approved["mouse_accession"].isin(
        m_counts[m_counts > 1].index), "ortholog_type"] = "many2one"

    many_h = h_counts[h_counts > 1].index.tolist()
    many_m = m_counts[m_counts > 1].index.tolist()
    if many_h or many_m:
        print(f"\n  WARNING: non-1:1 relationships detected:")
        if many_h:
            print(f"    Human with >1 mouse partner: {many_h[:5]}")
        if many_m:
            print(f"    Mouse with >1 human partner: {many_m[:5]}")

    # Save all approved pairs (backward-compatible)
    out_cols = [c for c in
                ["human_accession","human_gene","human_group",
                 "mouse_accession","mouse_gene","mouse_group",
                 "ortholog_type","confidence_tier","candidate_score",
                 "gene_similarity","same_mmseqs_cluster","same_protein_group"]
                if c in approved.columns]
    approved_out = approved[out_cols]
    approved_out.to_csv(data_dir / "ortholog_pairs.tsv", sep="\t", index=False)

    # ── Stratified 80/20 split ────────────────────────────────────────────────
    if len(approved) < 5:
        print("\n  WARNING: fewer than 5 pairs — skipping holdout split.")
        print("  All pairs written to both calibration and test files.")
        cal_df  = approved_out
        test_df = approved_out
    else:
        cal_raw, test_raw = _stratified_split(
            approved, test_frac=float(args.test_frac), seed=int(args.seed))
        cal_df  = cal_raw[out_cols]
        test_df = test_raw[out_cols]

    cal_df.to_csv(data_dir / "ortholog_pairs_calibration.tsv", sep="\t", index=False)
    test_df.to_csv(data_dir / "ortholog_pairs_test.tsv",        sep="\t", index=False)

    # Split summary
    summary_lines = [
        "=" * 60,
        "  BLIND TEST SPLIT SUMMARY",
        "=" * 60,
        f"  Total approved pairs  : {len(approved_out)}",
        f"  Calibration set (80%) : {len(cal_df)}  → ortholog_pairs_calibration.tsv",
        f"  Test set (20%)        : {len(test_df)}  → ortholog_pairs_test.tsv",
        "",
        "  *** DO NOT USE ortholog_pairs_test.tsv until step 8 ***",
        "  *** Tuning any parameter on the test set invalidates the evaluation ***",
        "",
        "  Tier distribution in each split:",
    ]
    for tier in sorted(approved["confidence_tier"].unique()):
        n_cal  = (cal_df.get("confidence_tier","") == tier).sum() if "confidence_tier" in cal_df.columns else "?"
        n_test = (test_df.get("confidence_tier","") == tier).sum() if "confidence_tier" in test_df.columns else "?"
        summary_lines.append(f"    {tier:<20}  cal={n_cal}  test={n_test}")
    summary_lines += [
        "",
        "  Ortholog type distribution (all pairs):",
        f"    one2one  : {(approved_out.get('ortholog_type','') == 'one2one').sum() if 'ortholog_type' in approved_out.columns else '?'}",
        f"    one2many : {(approved_out.get('ortholog_type','') == 'one2many').sum() if 'ortholog_type' in approved_out.columns else '?'}",
        f"    many2one : {(approved_out.get('ortholog_type','') == 'many2one').sum() if 'ortholog_type' in approved_out.columns else '?'}",
        "=" * 60,
    ]
    summary = "\n".join(summary_lines)
    print("\n" + summary)
    (data_dir / "split_summary.txt").write_text(summary, encoding="utf-8")

    print(f"\nFiles saved to {data_dir}:")
    print("  ortholog_pairs.tsv              (all pairs — backward-compatible)")
    print("  ortholog_pairs_calibration.tsv  (80% — for steps 06, 07)")
    print("  ortholog_pairs_test.tsv         (20% — LOCKED until step 08)")
    print("\nStep 5 complete.  Next: python 06_build_annotation_tables.py")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Finalize ortholog pairs with holdout split")
    p.add_argument("--data",      default="transfer_pipeline/data")
    p.add_argument("--test_frac", default=0.20, type=float,
                   help="Fraction reserved for blind test (default 0.20)")
    p.add_argument("--seed",      default=42,   type=int,
                   help="Random seed for split reproducibility")
    main(p.parse_args())
