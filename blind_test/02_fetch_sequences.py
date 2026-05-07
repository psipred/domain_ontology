"""
02_fetch_sequences.py
────────────────────────────────────────────────────────────────────────────
Step 2 — Fetch UniProt sequences for both species.

Outputs
───────
  human_camp_proteins.fasta    human sequences with species/gene/group tags
  mouse_camp_proteins.fasta    mouse sequences with species/gene/group tags
  combined_camp_proteins.fasta both species combined (input for MMseqs2)
  fetch_report.tsv             per-accession status (OK / MISSING / ERROR)

FASTA header format
────────────────────
  >sp|P07550|ADRB2_HUMAN species=human gene=ADRB2 group=GPCR_class_A
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations
import argparse
import time
from pathlib import Path
import pandas as pd
import requests

UNIPROT    = "https://rest.uniprot.org/uniprotkb/accessions"
BATCH_SIZE = 50
SLEEP      = 1.0
MAX_RETRY  = 3


def _fetch_batch(accessions: list[str]) -> dict[str, str]:
    for attempt in range(1, MAX_RETRY + 1):
        try:
            r = requests.get(
                UNIPROT,
                params={"accessions": ",".join(accessions),
                        "format": "fasta", "size": len(accessions)},
                timeout=60,
            )
            r.raise_for_status()
            results: dict[str, str] = {}
            current_acc, current_lines = None, []
            for line in r.text.splitlines():
                if line.startswith(">"):
                    if current_acc and current_lines:
                        results[current_acc] = "\n".join(current_lines)
                    parts       = line.split("|")
                    current_acc = (parts[1] if len(parts) >= 2
                                   else line[1:].split()[0])
                    current_lines = [line]
                elif line.strip():
                    current_lines.append(line)
            if current_acc and current_lines:
                results[current_acc] = "\n".join(current_lines)
            return results
        except Exception as e:
            if attempt == MAX_RETRY:
                raise
            print(f"    Retry {attempt}/{MAX_RETRY} after error: {e}")
            time.sleep(SLEEP * attempt)
    return {}


def _tag_header(fasta_record: str, meta: dict) -> str:
    lines    = fasta_record.split("\n")
    species  = meta.get("species", "unknown")
    gene     = meta.get("gene_symbol", "")
    group    = meta.get("protein_group", "other")
    tag      = f" species={species}"
    if gene:
        tag += f" gene={gene}"
    tag     += f" group={group}"
    lines[0] = lines[0] + tag
    return "\n".join(lines)


def fetch_species(proteins_tsv: Path, species: str,
                  out_fasta: Path) -> tuple[list[dict], dict[str, str]]:
    if not proteins_tsv.exists():
        print(f"  WARNING: {proteins_tsv} not found — skipping {species}")
        return [], {}

    df  = pd.read_csv(proteins_tsv, sep="\t", dtype=str)
    acc_col = next((c for c in ["protein_accession","accession"]
                    if c in df.columns), None)
    if acc_col is None:
        print(f"  ERROR: no accession column in {proteins_tsv}")
        return [], {}

    accessions = df[acc_col].dropna().str.strip().unique().tolist()
    meta_map   = {
        row[acc_col]: {
            "species":       str(row.get("species", species)),
            "gene_symbol":   str(row.get("gene_symbol", "")),
            "protein_group": str(row.get("protein_group", "other")),
        }
        for _, row in df.iterrows()
    }

    print(f"  Fetching {len(accessions)} {species} sequences …")
    all_seqs: dict[str, str] = {}
    report_rows: list[dict] = []

    for i in range(0, len(accessions), BATCH_SIZE):
        batch = accessions[i:i + BATCH_SIZE]
        print(f"    Batch {i//BATCH_SIZE+1}: {batch[0]} … {batch[-1]}", end="  ")
        try:
            fetched = _fetch_batch(batch)
            tagged  = {
                acc: _tag_header(rec, meta_map.get(acc, {}))
                for acc, rec in fetched.items()
            }
            all_seqs.update(tagged)
            print(f"got {len(fetched)}/{len(batch)}")
            for acc in batch:
                report_rows.append({
                    "accession": acc, "species": species,
                    "status": "OK" if acc in fetched else "MISSING",
                    "note":   "" if acc in fetched else "not in UniProt response",
                })
        except Exception as e:
            print(f"ERROR: {e}")
            for acc in batch:
                report_rows.append({"accession": acc, "species": species,
                                     "status": "ERROR", "note": str(e)})
        time.sleep(SLEEP)

    with out_fasta.open("w") as fh:
        for rec in all_seqs.values():
            fh.write(rec + "\n")
    print(f"  Saved {len(all_seqs)} sequences → {out_fasta}")
    return report_rows, all_seqs


def main(args: argparse.Namespace) -> None:
    out_dir  = Path(args.out)
    data_dir = Path(args.data)
    out_dir.mkdir(parents=True, exist_ok=True)
    print("=" * 60)
    print("  02 — Fetch Sequences")
    print("=" * 60)

    all_report:    list[dict]   = []
    combined_seqs: dict[str, str] = {}

    print("\nHuman proteins:")
    report_h, seqs_h = fetch_species(
        data_dir / "human_proteins.tsv", "human",
        out_dir  / "human_camp_proteins.fasta")
    all_report.extend(report_h)
    combined_seqs.update(seqs_h)

    print("\nMouse proteins:")
    report_m, seqs_m = fetch_species(
        data_dir / "mouse_proteins.tsv", "mouse",
        out_dir  / "mouse_camp_proteins.fasta")
    all_report.extend(report_m)
    combined_seqs.update(seqs_m)

    combined_path = out_dir / "combined_camp_proteins.fasta"
    with combined_path.open("w") as fh:
        for rec in combined_seqs.values():
            fh.write(rec + "\n")
    print(f"\nCombined FASTA → {combined_path}  ({len(combined_seqs)} sequences)")

    report_df = pd.DataFrame(all_report)
    report_df.to_csv(out_dir / "fetch_report.tsv", sep="\t", index=False)

    if report_df.empty or "status" not in report_df.columns:
        # No fetch attempts were made — protein TSV files were not found.
        # Likely cause: step 01 not run yet, or --data points to wrong directory.
        print("\nWARNING: No sequences were fetched (fetch report is empty).")
        print(f"  Expected protein tables at: {data_dir}/human_proteins.tsv")
        print("  Fix: run step 01 first:")
        print("    python 01_prepare_protein_tables.py --human <human.tsv> --mouse <mouse.tsv>")
    else:
        missing = report_df[report_df["status"] != "OK"]
        if not missing.empty:
            print(f"\nWarning: {len(missing)} sequences not retrieved:")
            for _, row in missing.iterrows():
                print(f"  {row['accession']} [{row['species']}] {row['note']}")
        else:
            print("All sequences retrieved successfully.")

    print("\nStep 2 complete.")
    print("Next: upload combined_camp_proteins.fasta to Tübingen MMseqs2 toolkit")
    print("      (min seq id = 0.50, coverage = 0.80)")
    print("      then save cluster result and run 03_parse_mmseqs2_results.py")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Fetch UniProt sequences")
    p.add_argument("--data", default="transfer_pipeline/data",
                   help="Directory containing human/mouse_proteins.tsv")
    p.add_argument("--out",  default="transfer_pipeline/data")
    main(p.parse_args())
