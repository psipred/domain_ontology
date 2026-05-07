#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import io
import re
import requests
import pandas as pd

# Choose ONE query:
# 1) reviewed human proteins
QUERY = "(reviewed:true) AND (organism_id:9606)"

# 2) human proteome
# QUERY = "proteome:UP000005640"

URL = "https://rest.uniprot.org/uniprotkb/stream"
PARAMS = {
    "query": QUERY,
    "format": "tsv",
    "fields": "accession,xref_interpro",
}

OUTFILE = "input/background_protein_domains.tsv"


def normalize_ipr(x: str) -> list[str]:
    if pd.isna(x):
        return []
    s = str(x).strip()
    if not s:
        return []
    # UniProt may return multiple InterPro IDs in one cell
    parts = re.split(r"[;,]\s*", s)
    out = []
    for p in parts:
        p = p.strip()
        if not p:
            continue
        if p.startswith("InterPro:"):
            p = p.replace("InterPro:", "INTERPRO:")
        elif re.fullmatch(r"IPR\d+", p):
            p = f"INTERPRO:{p}"
        out.append(p)
    return sorted(set(out))


def main():
    r = requests.get(URL, params=PARAMS, timeout=300)
    r.raise_for_status()

    df = pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str).fillna("")
    print("Downloaded columns:", list(df.columns))

    # Try to find the accession + InterPro columns robustly
    acc_col = None
    ipr_col = None

    for c in df.columns:
        cl = c.lower()
        if cl in {"entry", "accession"}:
            acc_col = c
        if "interpro" in cl:
            ipr_col = c

    if acc_col is None or ipr_col is None:
        raise ValueError(
            f"Could not find accession/InterPro columns. Found columns: {list(df.columns)}"
        )

    rows = []
    for _, row in df.iterrows():
        acc = str(row[acc_col]).strip()
        if not acc:
            continue
        for ipr in normalize_ipr(row[ipr_col]):
            rows.append((acc, ipr))

    out = pd.DataFrame(rows, columns=["protein_accession", "domain_type_id"])
    out = out.drop_duplicates().sort_values(["protein_accession", "domain_type_id"])
    out.to_csv(OUTFILE, sep="\t", index=False)

    print(f"Saved {len(out)} rows to {OUTFILE}")


if __name__ == "__main__":
    main()