#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fetch_sequences_calcium.py
──────────────────────────────────────────────────────────────────────────────
Fetches full UniProt FASTA sequences for every protein in the calcium
signaling pathway domain-occurrence table, then writes a FASTA file ready
for MMseqs2 clustering.

WORKFLOW POSITION
─────────────────
  [4] fetch_sequences_calcium.py  ← YOU ARE HERE
       ↓  calcium_signaling_proteins.fasta
  [5] Upload FASTA to MMseqs2 Tübingen Toolkit  (https://toolkit.bioinformatics.org)
       → Clustal Omega / MMseqs2 → download result TSV
  [6] convert_toolkit_output.py   → camp_family_cluster_ids.tsv
                                     protein_group_labels_sequencing.tsv
  [7] run_pipeline.py --stages 13,16

PATH SETTINGS — edit these three lines only
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations

import time
from pathlib import Path

import pandas as pd
import requests

# ── CONFIGURE THESE ──────────────────────────────────────────────────────────

# Your calcium signaling domain-occurrence table
OCCURRENCE_TABLE = (
    r"C:\Users\32045\Desktop\测试\integrated_pipeline _calcium_human\input\calcium_signaling_pathway_human_MERGED_domain_occurrences.tsv"
)

# Where to write the FASTA output  (feed this to MMseqs2 / Toolkit)
OUTPUT_FASTA   = (
    r"C:\Users\32045\Desktop\测试\integrated_pipeline _calcium_human\input"
    r"\calcium_signaling_proteins.fasta"
)

# Where to write the fetch report  (which accessions succeeded / failed)
OUTPUT_REPORT  = (
    r"C:\Users\32045\Desktop\测试\integrated_pipeline _calcium_human\input"
    r"\calcium_signaling_fetch_report.tsv"
)

# ── TUNEABLE PARAMETERS ───────────────────────────────────────────────────────
BATCH_SIZE         = 50        # accessions per UniProt API call
SLEEP_BETWEEN      = 1.0       # seconds between batches (be polite to UniProt)
MAX_RETRIES        = 3
TIMEOUT            = 60        # seconds per HTTP request

# ── INTERNAL CONSTANTS ────────────────────────────────────────────────────────
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"
HEADERS     = {"User-Agent": "calcium-signaling-seq-fetcher/1.0 (pipeline script)"}


# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def load_accessions(path: Path) -> list[str]:
    """
    Read the domain-occurrence TSV and return a deduplicated, sorted list of
    UniProt accessions from the protein_accession column.
    """
    df = pd.read_csv(path, sep="\t", dtype=str,
                     encoding="utf-8-sig", on_bad_lines="skip").fillna("")

    # Find the protein accession column (handles name variants)
    candidates = ["protein_accession", "accession", "Entry", "uniprot_accession"]
    col = next((c for c in candidates if c in df.columns), None)
    if col is None:
        raise ValueError(
            f"Cannot find a protein accession column in {path.name}.\n"
            f"Tried: {candidates}\nFound: {list(df.columns)}"
        )

    accs = sorted({
        str(v).strip()
        for v in df[col]
        if str(v).strip() and str(v).strip().lower() != "nan"
    })
    return accs


def parse_fasta_records(fasta_text: str) -> dict[str, str]:
    """
    Parse raw FASTA text returned by UniProt and return
    {bare_accession: full_fasta_record_string}.

    Handles both SwissProt (>sp|P07550|ADRB2_HUMAN ...) and
    TrEMBL   (>tr|A0A000|...) headers.
    Isoform suffixes (P12345-2) are stripped to the base accession.
    """
    records: dict[str, str] = {}
    current_acc: str | None = None
    current_lines: list[str] = []

    for line in fasta_text.splitlines():
        if line.startswith(">"):
            # Flush previous record
            if current_acc and current_lines:
                records[current_acc] = "\n".join(current_lines)
            # Parse accession from header
            parts = line.split("|")
            if len(parts) >= 2:
                acc = parts[1].strip()
            else:
                acc = line[1:].split()[0]
            current_acc = acc.split("-")[0]   # strip isoform suffix
            current_lines = [line]
        elif line.strip():
            current_lines.append(line.strip())

    if current_acc and current_lines:
        records[current_acc] = "\n".join(current_lines)

    return records


def fetch_batch(accessions: list[str], session: requests.Session) -> dict[str, str]:
    """
    Fetch FASTA for a batch of accessions from the UniProt REST API.
    Retries up to MAX_RETRIES times on transient errors.
    """
    query  = " OR ".join(f"accession:{a}" for a in accessions)
    params = {"query": query, "format": "fasta", "size": len(accessions)}

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = session.get(UNIPROT_URL, params=params,
                               headers=HEADERS, timeout=TIMEOUT)
            resp.raise_for_status()
            return parse_fasta_records(resp.text)
        except requests.RequestException as exc:
            if attempt == MAX_RETRIES:
                raise RuntimeError(
                    f"UniProt request failed after {MAX_RETRIES} attempts: {exc}"
                ) from exc
            time.sleep(attempt * 2)   # exponential back-off

    return {}   # unreachable


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    occ_path    = Path(OCCURRENCE_TABLE)
    fasta_path  = Path(OUTPUT_FASTA)
    report_path = Path(OUTPUT_REPORT)

    # ── 1. Load accessions ────────────────────────────────────────────────────
    if not occ_path.exists():
        raise FileNotFoundError(f"Occurrence table not found:\n{occ_path}")

    accessions = load_accessions(occ_path)
    print(f"Loaded {len(accessions)} unique protein accessions from {occ_path.name}")
    print(f"Sample: {accessions[:6]} …\n")

    # ── 2. Fetch FASTA sequences in batches ───────────────────────────────────
    all_sequences: dict[str, str] = {}
    report_rows:   list[dict]     = []

    with requests.Session() as session:
        n_batches = (len(accessions) + BATCH_SIZE - 1) // BATCH_SIZE
        for batch_idx, start in enumerate(range(0, len(accessions), BATCH_SIZE)):
            batch   = accessions[start : start + BATCH_SIZE]
            batch_n = batch_idx + 1
            print(f"  Batch {batch_n}/{n_batches}: {batch[0]} … {batch[-1]}", end="  ")

            try:
                fetched = fetch_batch(batch, session)
                all_sequences.update(fetched)
                got = len(fetched)
                print(f"→  {got}/{len(batch)} retrieved")

                for acc in batch:
                    if acc in fetched:
                        report_rows.append({"accession": acc, "status": "OK",   "note": ""})
                    else:
                        # Could be obsolete / merged accession
                        report_rows.append({"accession": acc, "status": "MISSING",
                                            "note": "not returned by UniProt — may be obsolete or merged"})
            except Exception as exc:
                print(f"→  ERROR: {exc}")
                for acc in batch:
                    report_rows.append({"accession": acc, "status": "ERROR", "note": str(exc)})

            if batch_idx < n_batches - 1:
                time.sleep(SLEEP_BETWEEN)

    # ── 3. Write FASTA output ─────────────────────────────────────────────────
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    with fasta_path.open("w", encoding="utf-8") as fh:
        for acc in accessions:
            if acc in all_sequences:
                fh.write(all_sequences[acc])
                fh.write("\n")
                n_written += 1

    print(f"\n{'='*60}")
    print(f"  FASTA written  →  {fasta_path}")
    print(f"  Sequences written : {n_written} / {len(accessions)}")

    # ── 4. Write fetch report ─────────────────────────────────────────────────
    report_df = pd.DataFrame(report_rows)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_df.to_csv(report_path, sep="\t", index=False)
    print(f"  Fetch report      →  {report_path}")

    # ── 5. Print any failures ─────────────────────────────────────────────────
    problems = report_df[report_df["status"] != "OK"]
    if problems.empty:
        print("\n  All accessions retrieved successfully.")
    else:
        print(f"\n  WARNING: {len(problems)} accession(s) not retrieved:")
        for _, row in problems.iterrows():
            print(f"    {row['accession']}  [{row['status']}]  {row['note']}")

    print(f"\n{'='*60}")
    print("  NEXT STEPS")
    print(f"{'='*60}")
    print()
    print("  1. Upload FASTA to the Tübingen Bioinformatics Toolkit:")
    print("       https://toolkit.bioinformatics.org  →  MMseqs2 Cluster")
    print(f"       File to upload: {fasta_path.name}")
    print()
    print("  2. Recommended MMseqs2 settings:")
    print("       Seq. identity threshold  : 0.50  (or 0.40 for broader clusters)")
    print("       Min. alignment coverage  : 0.80")
    print("       Coverage mode            : bidirectional")
    print()
    print("  3. Download the cluster result TSV from the Toolkit job page.")
    print()
    print("  4. Run convert_toolkit_output.py to convert to pipeline format:")
    print("       → camp_family_cluster_ids.tsv")
    print("       → protein_group_labels_sequencing.tsv")
    print()
    print("  5. Run the pipeline:")
    print("       python run_pipeline.py --stages 13,16")
    print()


if __name__ == "__main__":
    main()
