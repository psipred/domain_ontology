#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
01_data_and_merge.py
Domain-Centric Protein Ontology Pipeline  —  Stage 1 + 2

Stage 1 — Data acquisition
    Step1: Queries UniProtKB REST API for reviewed Swiss-Prot entries.
        Fetches per-protein domain occurrence data from the InterPro API.
        Resolves overlapping domain signatures using the longest-representative overlap-group strategy.
        Builds structured ontology-aligned TSV tables per gene.
    Step 2: Pathway-level aggregation  (formerly merged_features_summary.py)
        Merges all per-gene occurrence and domain-type tables into pathway-level TSVs.
        Generates a domain_type_summary.tsv with per-type statistics (occurrence
        counts, length statistics, positional/topology counts, GO/KEGG/Reactome
        unions) — required by downstream Stage 4B (transfer.py).
        Optionally merges link tables (domain order edges, binding, structural,
        and signature links).

Stage 2
    Two columns in the domain occurrence table require post-pipeline manual
    curation before proceeding to Stage 3 (ID injection):
      function_tags   — domain-level functional role  (DOF vocabulary)
      binding_partner — known interaction targets      (DOT vocabulary)

Usage examples
──────────────
  # Human cAMP pathway (batch, gene file)
  python 01_data_and_merge.py \\
      --gene-file genes/cAMP_human.txt \\
      --pathway-name cAMP_human \\
      --organism-id 9606 \\
      --fetch-interpro

  # Mouse orthologs
  python 01_data_and_merge.py \\
      --gene-file genes/cAMP_human.txt \\
      --pathway-name cAMP_mouse \\
      --organism-id 10090 \\
      --fetch-interpro

  # Include link tables in merged output
  python 01_data_and_merge.py \\
      --gene-file genes/cAMP_human.txt \\
      --pathway-name cAMP_human \\
      --organism-id 9606 \\
      --fetch-interpro \\
      --merge-links

  # Single protein (debug / inspection)
  python 01_data_and_merge.py \\
      --name PDE4D \\
      --query "gene:PDE4D* AND organism_id:9606" \\
      --fetch-interpro
"""

from __future__ import annotations

__version__ = "1.2.0"

# ─────────────────────────────────────────────────────────────────────────────
# Standard-library imports
# ─────────────────────────────────────────────────────────────────────────────
import argparse
import ast
import csv
import json
import os
import re
import sys
import time
import traceback
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import urljoin

# ─────────────────────────────────────────────────────────────────────────────
# Third-party imports
# ─────────────────────────────────────────────────────────────────────────────
import requests

try:
    import pandas as pd
    from pandas.errors import EmptyDataError
    _PANDAS_OK = True
except ImportError:
    pd = None             # type: ignore[assignment]
    EmptyDataError = Exception  # type: ignore[misc,assignment]
    _PANDAS_OK = False


# ═══════════════════════════════════════════════════════════════════════════════
# PART 1 — Constants and API endpoints
# ═══════════════════════════════════════════════════════════════════════════════

UNIPROT_SEARCH           = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY_TXT        = "https://rest.uniprot.org/uniprotkb/{acc}.txt"

INTERPRO_BASE            = "https://www.ebi.ac.uk"
INTERPRO_PROT_URL_UNIPROT    = "https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{acc}"
INTERPRO_PROT_URL_REVIEWED   = "https://www.ebi.ac.uk/interpro/api/entry/all/protein/reviewed/{acc}"
INTERPRO_PROT_URL_UNREVIEWED = "https://www.ebi.ac.uk/interpro/api/entry/all/protein/unreviewed/{acc}"

UA               = f"domain-ontology-pipeline/{__version__}"
DEFAULT_DELAY    = 0.34   # seconds between API requests (respects rate limits)
DEFAULT_RETRIES  = 3
DEFAULT_BACKOFF  = 1.7
DEFAULT_ORGANISM = "9606" # NCBI taxonomy ID for Homo sapiens

OCC_RANGE_DELIM  = ".."   # avoids Excel auto-converting "1-46" to "Jan-46"

POSITIONAL_CATEGORIES = ("Nterminal", "Internal", "Cterminal", "unknown")

# Regex helpers (data-acquisition layer)
REGION_ITEM_SPLIT = re.compile(r"\s*\d+\)\s*")
RANGE_RE  = re.compile(r"(\d+)\s*-\s*(\d+)")
LEN_RE    = re.compile(r"len\s*=\s*(\d+)", re.I)
DESC_RE   = re.compile(r"desc\s*=\s*([^;|]+)", re.I)
SEQ_RE    = re.compile(r"seq\s*=\s*([^;|]+)", re.I)
NAME_RE   = re.compile(r"name\s*=\s*([^;|]+)", re.I)

# ── Glob patterns for pathway-level aggregation (Stage 2) ────────────────────
# These match the directory layout produced by Stage 1:
#   <root>/<GENE>_out_final/out/<GENE>_domain_occurrences.tsv  etc.
PAT_OCC    = "*_out_final/out/*_domain_occurrences.tsv"
PAT_TYPES  = "*_out_final/out/*_domain_types.tsv"
PAT_ORDER  = "*_out_final/out/*_domain_order_edges.tsv"
PAT_BIND   = "*_out_final/out/*_occurrence_binding_links.tsv"
PAT_STRUCT = "*_out_final/out/*_occurrence_structural_links.tsv"
PAT_SIG    = "*_out_final/out/*_occurrence_signature_links.tsv"

# ── Deduplication key for occurrence rows ────────────────────────────────────
# The pipeline fetches proteins via gene-name prefix queries (e.g. "ADCY*").
# When multiple queries match the same protein accession (e.g. ADCYAP1R1 is
# returned by both "ADCY*" and "ADCYAP1R1*"), the same domain occurrence row
# appears in multiple per-gene output files and is concatenated into duplicates
# during Stage 2.  DEDUP_KEY uniquely identifies one physical domain occurrence
# — same protein, same coordinates, same InterPro entry — and is used to drop
# these duplicates after the per-gene concat.
DEDUP_KEY: List[str] = [
    "protein_accession",  # UniProt accession
    "start",              # domain start residue (1-based)
    "end",                # domain end residue   (1-based)
    "domain_type_id",     # InterPro entry ID (e.g. INTERPRO:IPR000532)
]

# ── Minimum column sets enforced during merge ─────────────────────────────────
# Missing columns are silently added as empty strings so downstream code
# never crashes on a KeyError.
OCC_MIN_COLS = [
    "protein_accession", "domain_type_id", "label", "label_db", "label_acc",
    "start", "end", "copy_number",
    "positional_category", "topology_context", "proximity_categories",
    "go_terms_all", "go_terms_mf", "go_terms_bp", "go_terms_cc",
    "function_tags", "binding_target_chebi", "signatures", "cath_superfamilies",
    # reactome_ids and kegg_ids removed from occurrence output (protein-level only)
]

TYPE_MIN_COLS = [
    "domain_type_id", "source_db", "accession", "name", "description",
    "go_terms_all", "go_terms_mf", "go_terms_bp", "go_terms_cc",
    "function_tags", "evidence_level",
]

# Multi-value splitter used in aggregation helpers
_SPLIT_MULTI = re.compile(r"[;,|]\s*")
_SPLIT_WS    = re.compile(r"[\s,]+")


# ═══════════════════════════════════════════════════════════════════════════════
# PART 2 — Gene-file loader
# ═══════════════════════════════════════════════════════════════════════════════

def load_gene_file(path: str) -> List[str]:
    """
    Load a gene list from a plain-text or tabular file.

    Formats
    -------
    Plain text (.txt or any extension)
        One HGNC gene symbol per line.  Lines beginning with '#' are comments.

    TSV / CSV
        Must contain a column named 'gene' (case-insensitive).
        All other columns are ignored.

    Returns
    -------
    Deduplicated list of gene name strings, preserving input order.
    """
    p = Path(path)
    if not p.exists():
        raise SystemExit(
            f"[FATAL] Gene file not found: {path}\n"
            "Provide a plain-text file (one gene per line) or a TSV/CSV "
            "with a 'gene' column."
        )

    suffix = p.suffix.lower()
    genes: List[str] = []

    if suffix in (".tsv", ".csv"):
        sep = "\t" if suffix == ".tsv" else ","
        try:
            import csv as _csv
            with open(p, newline="", encoding="utf-8") as fh:
                reader = _csv.DictReader(fh, delimiter=sep)
                fieldnames = reader.fieldnames or []
                gene_col = next(
                    (f for f in fieldnames if f.strip().lower() == "gene"), None
                )
                if gene_col is None:
                    raise SystemExit(
                        f"[FATAL] No 'gene' column found in {path}.\n"
                        f"Available columns: {fieldnames}\n"
                        "Rename the gene column to 'gene' or use a plain-text file."
                    )
                for row in reader:
                    val = (row.get(gene_col) or "").strip()
                    if val and not val.startswith("#"):
                        genes.append(val)
        except SystemExit:
            raise
        except Exception as exc:
            raise SystemExit(f"[FATAL] Could not read {path}: {exc}") from exc
    else:
        with open(p, encoding="utf-8") as fh:
            for line in fh:
                stripped = line.strip()
                if stripped and not stripped.startswith("#"):
                    genes.append(stripped)

    seen: set = set()
    unique: List[str] = []
    for g in genes:
        if g not in seen:
            seen.add(g)
            unique.append(g)

    if not unique:
        raise SystemExit(
            f"[FATAL] No gene names found in {path}.\n"
            "Check the file format — see --help for accepted formats."
        )

    print(f"[INFO] Loaded {len(unique)} gene(s) from {path}")
    return unique


# ═══════════════════════════════════════════════════════════════════════════════
# PART 3 — UniProt + InterPro data-acquisition layer
# ═══════════════════════════════════════════════════════════════════════════════

def _norm_space(s: str) -> str:
    return re.sub(r"\s+", " ", s or "").strip()

def _safe_int(x) -> Optional[int]:
    try:
        return int(x)
    except Exception:
        return None

def _http_get(
        session: requests.Session,
        url: str,
        params: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        timeout: int = 30,
        expect: str = "json",
        max_retries: int = DEFAULT_RETRIES,
        backoff: float = DEFAULT_BACKOFF,
        delay: float = DEFAULT_DELAY,
) -> Any:
    h = headers or {}
    for attempt in range(max_retries):
        try:
            r = session.get(url, params=params, headers=h, timeout=timeout)
            if r.status_code in (429, 500, 502, 503, 504):
                wait = backoff ** attempt
                print(f"[WARN] GET {url} -> {r.status_code}; retry in {wait:.1f}s ({attempt+1}/{max_retries})")
                time.sleep(wait)
                continue
            r.raise_for_status()
            if expect == "json":
                ct = (r.headers.get("Content-Type") or "").lower()
                if "json" not in ct or not r.text.strip():
                    wait = backoff ** attempt
                    print(f"[WARN] Non-JSON/empty for {url} (CT={ct!r}); retry in {wait:.1f}s")
                    time.sleep(wait)
                    continue
                out = r.json()
            else:
                out = r.text
            time.sleep(delay)
            return out
        except (requests.RequestException, ValueError) as e:
            wait = backoff ** attempt
            print(f"[WARN] GET error for {url}: {e}; retry in {wait:.1f}s ({attempt+1}/{max_retries})")
            time.sleep(wait)
    raise RuntimeError(f"Failed GET after {max_retries} retries: {url}")

def _is_reviewed_search_item(item: Dict[str, Any]) -> Optional[bool]:
    t = (item.get("entryType") or "")
    if isinstance(t, str):
        tl = t.lower()
        if "reviewed" in tl and "unreviewed" not in tl:
            return True
        if "unreviewed" in tl:
            return False
    v = item.get("reviewed")
    if isinstance(v, bool):
        return v
    return None

def search_uniprot(
        query: str,
        reviewed: Optional[bool] = None,
        size: int = 300,
        delay: float = DEFAULT_DELAY,
) -> List[Dict[str, Any]]:
    session = requests.Session()
    session.headers.update({"Accept": "application/json", "User-Agent": UA})
    q = query
    if reviewed is True:
        q += " AND reviewed:true"
    elif reviewed is False:
        q += " AND reviewed:false"
    params = {"query": q, "format": "json", "size": size}
    data = _http_get(session, UNIPROT_SEARCH, params=params, expect="json", delay=delay)
    return data.get("results", []) or []

def fetch_uniprot_txt(acc: str, delay: float = DEFAULT_DELAY) -> str:
    session = requests.Session()
    session.headers.update({"Accept": "text/plain", "User-Agent": UA})
    url = UNIPROT_ENTRY_TXT.format(acc=acc)
    return _http_get(session, url, expect="text", delay=delay)

def _parse_sequence(lines: List[str]) -> str:
    seq, in_seq = [], False
    for ln in lines:
        if ln.startswith("SQ "):
            in_seq = True
            continue
        if in_seq:
            if ln.startswith("//"):
                break
            seq.append(re.sub(r"[^A-Za-z]", "", ln))
    return "".join(seq)

def _parse_length_from_sq_header(lines: List[str]) -> Optional[int]:
    for ln in lines:
        if ln.startswith("SQ "):
            m = re.search(r"(\d+)\s+AA;", ln)
            if m:
                return int(m.group(1))
            break
    return None

def _parse_id(lines: List[str]) -> str:
    for ln in lines:
        if ln.startswith("ID "):
            parts = _norm_space(ln[2:]).split()
            return parts[0] if parts else ""
    return ""

def _parse_gn(lines: List[str]) -> str:
    names: List[str] = []
    for ln in lines:
        if ln.startswith("GN   "):
            for m in re.finditer(r"Name=([^;{]+)", ln):
                names.append(m.group(1).strip())
            for m in re.finditer(r"Synonyms=([^;{]+)", ln):
                names.extend([s.strip() for s in m.group(1).split(",")])
    return ";".join([n for n in names if n])

def _parse_protein_name(lines: List[str]) -> str:
    de = " ".join(ln[5:].strip() for ln in lines if ln.startswith("DE   "))
    for pattern in [r"RecName:\s*Full=([^;]+);", r"SubName:\s*Full=([^;]+);", r"Full=([^;]+);"]:
        m = re.search(pattern, de)
        if m:
            return m.group(1)
    return ""

def _parse_dr(lines: List[str]) -> Dict[str, List[str]]:
    xrefs: Dict[str, List[str]] = {
        "Pfam": [], "SMART": [], "PROSITE": [], "CATH": [],
        "GO": [], "GO_MF": [], "GO_BP": [], "GO_CC": [],
        "Reactome": [], "KEGG": [],
    }
    go_all, go_mf, go_bp, go_cc = set(), set(), set(), set()
    for ln in lines:
        if not ln.startswith("DR   "):
            continue
        rest = ln[5:].strip()
        if rest.startswith("GO;"):
            parts = [p.strip() for p in rest.split(";")]
            if len(parts) >= 3:
                go_id = parts[1]
                aspect = parts[2].split(":", 1)[0].strip()
                if go_id.startswith("GO:"):
                    go_all.add(go_id)
                    if aspect == "F":   go_mf.add(go_id)
                    elif aspect == "P": go_bp.add(go_id)
                    elif aspect == "C": go_cc.add(go_id)
            continue
        if rest.startswith("Reactome;"):
            m = re.search(r"Reactome;\s*([^;\s]+)", rest)
            if m: xrefs["Reactome"].append(m.group(1))
            continue
        if rest.startswith("KEGG;"):
            m = re.search(r"KEGG;\s*([^;\s]+)", rest)
            if m: xrefs["KEGG"].append(m.group(1))
            continue
        for db in ("Pfam", "SMART", "PROSITE", "CATH"):
            if rest.startswith(f"{db};"):
                m = re.search(rf"{db};\s*([^;\s]+)", rest)
                if m: xrefs[db].append(m.group(1))
                break
    xrefs["GO"]    = sorted(go_all)
    xrefs["GO_MF"] = sorted(go_mf)
    xrefs["GO_BP"] = sorted(go_bp)
    xrefs["GO_CC"] = sorted(go_cc)
    return xrefs

def _parse_cc_sections(lines: List[str]) -> Dict[str, str]:
    sections: Dict[str, List[str]] = {}
    current, buf = None, []
    for ln in lines:
        if not ln.startswith("CC   "):
            continue
        txt = ln[5:].rstrip()
        if txt.startswith("-!-"):
            if current:
                sections[current] = _norm_space(" ".join(buf))
                buf = []
            m = re.match(r"-!-\s*([^:]+):\s*(.*)", txt)
            current = m.group(1).upper() if m else None
            buf.append(m.group(2) if m else "")
        else:
            if current is not None:
                buf.append(txt)
    if current:
        sections[current] = _norm_space(" ".join(buf))
    return {k: sections.get(k, "") for k in
            ("FUNCTION", "COFACTOR", "ACTIVITY REGULATION",
             "SUBUNIT", "INTERACTION", "SUBCELLULAR LOCATION")}

def _parse_ft_features(lines: List[str]) -> Dict[str, Any]:
    feats: List[Dict[str, Any]] = []
    cur: Optional[Dict[str, Any]] = None
    header_re = re.compile(r"^(\w+)\s+([<>]?\d+)(?:\.\.([<>]?\d+))?\s*(.*)$")

    def _to_int(x: Optional[str]) -> Optional[int]:
        if not x: return None
        m = re.search(r"\d+", x)
        return int(m.group(0)) if m else None

    for ln in lines:
        if not ln.startswith("FT   "):
            continue
        body = ln[5:].rstrip("\n")
        m = header_re.match(body.strip())
        if m:
            ftype = m.group(1).upper()
            s = _to_int(m.group(2))
            e = _to_int(m.group(3)) if m.group(3) else s
            cur = {"type": ftype, "start": s, "end": e,
                   "note": _norm_space(m.group(4) or ""), "ftid": "", "qualifiers": {}}
            feats.append(cur)
            continue
        if body.startswith(" " * 11) and cur is not None:
            kval = body.strip()
            qm = re.match(r'^/([^=]+)=(?:"([^"]*)"|(.+))$', kval)
            if qm:
                key = qm.group(1).strip()
                val = _norm_space(qm.group(2) if qm.group(2) is not None else qm.group(3) or "")
                q = cur["qualifiers"]
                if key in q:
                    q[key] = [q[key], val] if not isinstance(q[key], list) else q[key] + [val]
                else:
                    q[key] = val
                if key.lower() == "note":
                    cur["note"] = _norm_space((cur.get("note") or "") + " " + val)
                if key == "FTId":
                    cur["ftid"] = val
            else:
                if cur:
                    cur["note"] = _norm_space((cur.get("note") or "") + " " + kval)

    return {
        "features": feats,
        "regions":  [f for f in feats if f["type"] == "REGION"],
        "compbias": [f for f in feats if f["type"] == "COMPBIAS"],
        "domains":  [f for f in feats if f["type"] == "DOMAIN"],
        "motifs":   [f for f in feats if f["type"] == "MOTIF"],
        "region_count":  sum(1 for f in feats if f["type"] == "REGION"),
        "compbias_count":sum(1 for f in feats if f["type"] == "COMPBIAS"),
        "domain_count":  sum(1 for f in feats if f["type"] == "DOMAIN"),
        "motif_count":   sum(1 for f in feats if f["type"] == "MOTIF"),
    }

def parse_uniprot_txt(txt: str) -> Dict[str, Any]:
    lines = txt.splitlines()
    seq = _parse_sequence(lines)
    out: Dict[str, Any] = {
        "sequence": seq,
        "protein_length": len(seq) if seq else (_parse_length_from_sq_header(lines) or 0),
        "entry_id_from_txt": _parse_id(lines),
        "protein_name_from_txt": _parse_protein_name(lines),
        "genes_from_gn": _parse_gn(lines),
    }
    xrefs = _parse_dr(lines)
    out["cross_refs"] = {k: v for k, v in xrefs.items() if not k.startswith("GO")}
    out["go_ids"]    = xrefs.get("GO", [])
    out["go_mf_ids"] = xrefs.get("GO_MF", [])
    out["go_bp_ids"] = xrefs.get("GO_BP", [])
    out["go_cc_ids"] = xrefs.get("GO_CC", [])
    out.update(_parse_cc_sections(lines))
    out.update(_parse_ft_features(lines))
    out["cc_function"]             = out.get("FUNCTION", "")
    out["cc_cofactor"]             = out.get("COFACTOR", "")
    out["cc_regulation"]           = out.get("ACTIVITY REGULATION", "")
    out["cc_subunit"]              = out.get("SUBUNIT", "")
    out["cc_interaction"]          = out.get("INTERACTION", "")
    out["cc_subcellular_location"] = out.get("SUBCELLULAR LOCATION", "")
    return out

def extract_protein_name_string(item: Dict[str, Any]) -> str:
    try:
        rec = item.get("proteinDescription", {}).get("recommendedName", {})
        full = rec.get("fullName", {}).get("value")
        if full: return full
    except Exception:
        pass
    try:
        subs = item.get("proteinDescription", {}).get("submissionNames", [])
        if subs:
            full = subs[0].get("fullName", {}).get("value")
            if full: return full
    except Exception:
        pass
    return ""

def _length_from_search_item(item: Dict[str, Any]) -> Optional[int]:
    try:
        l = item.get("sequence", {}).get("length")
        return int(l) if l is not None else None
    except Exception:
        return None

def feature_pretty_list(feats: List[Dict[str, Any]], seq: str,
                        max_items: Optional[int] = None,
                        full_seq: bool = False) -> str:
    if not feats: return "0 items"
    rows = []
    for i, f in enumerate(feats, 1):
        s_i, t_i = _safe_int(f.get("start")), _safe_int(f.get("end"))
        subseq = ""
        if seq and isinstance(s_i, int) and isinstance(t_i, int) and 1 <= s_i <= t_i <= len(seq):
            subseq = seq[s_i-1:t_i]
            if not full_seq and len(subseq) > 40:
                subseq = subseq[:40] + "..."
        length = len(subseq.replace("...", "")) if subseq else (t_i - s_i + 1 if s_i and t_i else 0)
        rows.append(f"{i}) {s_i}-{t_i} | len={length} | seq={subseq} | desc={_norm_space(f.get('note', ''))}")
        if max_items and i >= max_items:
            rows.append(f"... ({len(feats)} total)")
            break
    return " ; ".join(rows)


def write_features_pretty_txt(outdir: str, acc: str,
                               feats: List[Dict[str, Any]], seq: str) -> None:
    path = os.path.join(outdir, f"features_pretty_{acc}.txt")
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"# Pretty feature list for {acc}\n")
        f.write(f"# Protein length: {len(seq) if seq else 0} AA\n\n")
        for ftr in feats:
            s_i, t_i = _safe_int(ftr.get("start")), _safe_int(ftr.get("end"))
            subseq = ""
            if seq and isinstance(s_i, int) and isinstance(t_i, int) and 1 <= s_i <= t_i <= len(seq):
                subseq = seq[s_i - 1:t_i]
            f.write(f"{ftr.get('type','')} | {s_i}-{t_i} | len={len(subseq)} | ftid={ftr.get('ftid','')}\n")
            f.write(f"desc: {ftr.get('note','')}\n")
            q = ftr.get("qualifiers") or {}
            if q:
                f.write(f"qualifiers: {json.dumps(q, ensure_ascii=False)}\n")
            f.write(f"seq : {subseq}\n\n")

# ── InterPro fetch ─────────────────────────────────────────────────────────────

def _interpro_fetch(url: str, page_size: int, delay: float,
                    max_retries: int, backoff: float) -> List[Dict[str, Any]]:
    entries: List[Dict[str, Any]] = []
    session = requests.Session()
    session.headers.update({"Accept": "application/json", "User-Agent": UA})
    params = {"format": "json", "page_size": page_size}
    while True:
        data = None
        for attempt in range(max_retries):
            try:
                r = session.get(url, params=params, timeout=30)
                if r.status_code in (429, 500, 502, 503, 504):
                    wait = backoff ** attempt
                    print(f"[WARN] InterPro {r.status_code}; retry in {wait:.1f}s ({attempt+1}/{max_retries})")
                    time.sleep(wait)
                    continue
                r.raise_for_status()
                ct = (r.headers.get("Content-Type") or "").lower()
                if "json" not in ct or not r.text.strip():
                    time.sleep(backoff ** attempt)
                    continue
                data = r.json()
                break
            except (requests.RequestException, ValueError) as e:
                print(f"[WARN] InterPro error: {e}; retry in {backoff**attempt:.1f}s")
                time.sleep(backoff ** attempt)
        if data is None:
            print(f"[ERROR] InterPro page failed for {url}. Stopping.")
            break
        entries.extend(data.get("results", []) or [])
        next_url = data.get("next") or (data.get("links") or {}).get("next")
        if not next_url:
            break
        if next_url.startswith("/"):
            next_url = urljoin(INTERPRO_BASE, next_url)
        url, params = next_url, {}
        time.sleep(delay)
    return entries

def interpro_entries_for_protein(
        accession: str,
        mode: str = "uniprot",
        page_size: int = 200,
        delay: float = 0.5,
        max_retries: int = DEFAULT_RETRIES,
        backoff: float = DEFAULT_BACKOFF,
) -> List[Dict[str, Any]]:
    kw = dict(page_size=page_size, delay=delay, max_retries=max_retries, backoff=backoff)
    if mode == "uniprot":
        return _interpro_fetch(INTERPRO_PROT_URL_UNIPROT.format(acc=accession), **kw)
    if mode == "reviewed":
        return _interpro_fetch(INTERPRO_PROT_URL_REVIEWED.format(acc=accession), **kw)
    if mode == "unreviewed":
        return _interpro_fetch(INTERPRO_PROT_URL_UNREVIEWED.format(acc=accession), **kw)
    if mode == "both":
        rev   = _interpro_fetch(INTERPRO_PROT_URL_REVIEWED.format(acc=accession), **kw)
        unrev = _interpro_fetch(INTERPRO_PROT_URL_UNREVIEWED.format(acc=accession), **kw)
        seen: set = set()
        out: List[Dict[str, Any]] = []
        for e in rev + unrev:
            meta = e.get("metadata") or {}
            acc0 = meta.get("accession") or ""
            locs = ((e.get("proteins") or [{}])[0].get("entry_protein_locations") or [])
            key = (
                acc0,
                locs[0]["fragments"][0].get("start") if locs and locs[0].get("fragments") else None,
                locs[0]["fragments"][0].get("end")   if locs and locs[0].get("fragments") else None,
            )
            if key not in seen:
                seen.add(key)
                out.append(e)
        return out
    return _interpro_fetch(INTERPRO_PROT_URL_UNIPROT.format(acc=accession), **kw)

# ── Interval helpers ───────────────────────────────────────────────────────────

def _interval_overlap_len(a: Tuple[int,int], b: Tuple[int,int]) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)

def _interval_distance(a: Tuple[int,int], b: Tuple[int,int]) -> int:
    if a[1] < b[0]: return b[0] - a[1]
    if b[1] < a[0]: return a[0] - b[1]
    return 0

def _segments_distance(segments: List[Tuple[int,int]], interval: Tuple[int,int]) -> Optional[int]:
    if not segments: return None
    return min(_interval_distance(seg, interval) for seg in segments)

def _parse_ranges_to_intervals(ranges_str: str) -> List[Tuple[int,int]]:
    out: List[Tuple[int,int]] = []
    for part in (ranges_str or "").split(";"):
        part = part.strip()
        m = re.match(r"^\s*(\d+)\s*(?:\.\.|-)\s*(\d+)\s*$", part)
        if m:
            s, t = int(m.group(1)), int(m.group(2))
            if s <= t: out.append((s, t))
    return out

def _collect_intervals(features: List[Dict[str,Any]], types: set) -> List[Tuple[int,int]]:
    out = []
    for f in features or []:
        if (f.get("type") or "").upper() in types:
            s, e = _safe_int(f.get("start")), _safe_int(f.get("end"))
            if isinstance(s, int) and isinstance(e, int) and s <= e:
                out.append((s, e))
    return out

def _signal_cleavage_pos(features: List[Dict[str,Any]]) -> Optional[int]:
    sigs = _collect_intervals(features, {"SIGNAL"})
    if not sigs: return None
    sigs.sort(key=lambda x: x[0])
    return sigs[0][1]

# ── InterPro parsing helpers ───────────────────────────────────────────────────

def _safe_meta(e: Dict[str, Any]) -> Dict[str, Any]:
    return (e.get("metadata") or e.get("entry_metadata") or e.get("integrated") or {}) or {}

def _entry_locations_for_protein(e: Dict[str, Any], uniprot_acc: str) -> List[Dict[str, Any]]:
    if isinstance(e.get("entry_protein_locations"), list):
        return e["entry_protein_locations"]
    want = (uniprot_acc or "").upper()
    for p in (e.get("proteins") or []):
        pacc = (p.get("accession") or p.get("protein_accession") or "").upper()
        if pacc == want:
            return p.get("entry_protein_locations") or []
    return []

def _parse_interpro_go(go_terms: Any) -> Tuple[List[str], List[str], List[str], List[str], List[str]]:
    go_all, go_mf, go_bp, go_cc = [], [], [], []
    func_tags: set = set()
    if not isinstance(go_terms, list):
        return [], [], [], [], []
    for g in go_terms:
        if not isinstance(g, dict): continue
        go_id = g.get("identifier") or g.get("accession") or g.get("id") or ""
        cat   = g.get("category") or g.get("aspect") or ""
        code  = ""
        if isinstance(cat, dict):
            code = (cat.get("code") or cat.get("name") or "").upper()
        elif isinstance(cat, str):
            code = cat.upper()
        if go_id.startswith("GO:"):
            go_all.append(go_id)
            if code in ("F", "MOLECULAR_FUNCTION", "MF"):
                go_mf.append(go_id); func_tags.add(go_id)
            elif code in ("P", "BIOLOGICAL_PROCESS", "BP"):
                go_bp.append(go_id)
            elif code in ("C", "CELLULAR_COMPONENT", "CC"):
                go_cc.append(go_id)
    return go_all, go_mf, go_bp, go_cc, list(func_tags)


# ── Delegate to new_pipeline.py for the heavy InterPro annotation logic ────────
# This keeps the file at a manageable size. To produce a truly single-file
# distribution, replace this block with the verbatim function bodies from
# new_pipeline.py (interpro_domain_occurrences, build_domain_instances_from_interpro,
# annotate_* functions, run_uniprot_summary, etc.).

# ── Fully inlined from new_pipeline.py ────────────────────────────────────────

def _member_db_rows_from_member_databases(interpro_ipr, interpro_name, member_databases):
    rows = []
    if not isinstance(member_databases, dict):
        return rows
    for db, mapping in member_databases.items():
        if isinstance(mapping, dict):
            for acc, nm in mapping.items():
                rows.append({"ipr": interpro_ipr, "interpro_name": interpro_name,
                             "member_db": db, "member_accession": acc,
                             "member_name": (nm if isinstance(nm, str) else ""),
                             "member_source_database": db, "member_type": ""})
    return rows


def interpro_domain_occurrences(entries, uniprot_acc, protein_seq):
    dom_rows: List[Dict[str, Any]] = []
    member_rows: List[Dict[str, Any]] = []
    for e in (entries or []):
        meta = _safe_meta(e)
        etype = (meta.get("type") or e.get("type") or e.get("entry_type") or "")
        if str(etype).strip().lower() != "domain":
            continue
        entry_acc = meta.get("accession") or meta.get("uid") or e.get("ac") or e.get("entry_acc") or ""
        name = meta.get("name") or e.get("name") or ""
        srcdb = (meta.get("source_database") or e.get("source_database") or "").strip()
        integrated = meta.get("integrated")
        desc = meta.get("description") or meta.get("abstract") or ""
        go_all, go_mf, go_bp, go_cc, func_tags = _parse_interpro_go(meta.get("go_terms") or e.get("go_terms"))
        if str(srcdb).strip().lower() == "interpro":
            member_rows.extend(_member_db_rows_from_member_databases(
                str(entry_acc), str(name), meta.get("member_databases")))
        locs = _entry_locations_for_protein(e, uniprot_acc)
        if not locs:
            continue
        for li, loc in enumerate(locs, 1):
            frags = loc.get("fragments") or []
            segs: List[Tuple[int, int]] = []
            sequences: List[str] = []
            for fr in frags:
                s = fr.get("start") or fr.get("begin") or fr.get("from")
                t = fr.get("end") or fr.get("to")
                if isinstance(s, int) and isinstance(t, int) and s <= t:
                    segs.append((s, t))
                    if protein_seq and 1 <= s <= t <= len(protein_seq):
                        sequences.append(protein_seq[s - 1:t])
            if not segs:
                continue
            start = min(s for s, _ in segs)
            end   = max(t for _, t in segs)
            ranges = ";".join([f"{s}..{t}" for s, t in segs])
            dom_rows.append({
                "ipr": str(entry_acc), "name": str(name), "type": "Domain",
                "entry_type": "domain", "source_database": str(srcdb),
                "integrated": (str(integrated) if integrated else ""),
                "occurrence_index": li, "n_fragments": len(segs),
                "ranges": ranges, "segments": segs, "start": start, "end": end,
                "sequences": ";".join(sequences),
                "description": _norm_space(str(desc)),
                "go_terms": ";".join(go_all), "go_terms_mf": ";".join(go_mf),
                "go_terms_bp": ";".join(go_bp), "go_terms_cc": ";".join(go_cc),
                "function_tags": ";".join(func_tags), "evidence_level": "InterProEntry",
            })
    return dom_rows, member_rows


def interpro_structural_hits(entries, uniprot_acc):
    def _frags(e):
        out = []
        for loc in (_entry_locations_for_protein(e, uniprot_acc) or []):
            for fr in (loc.get("fragments") or []):
                s = fr.get("start") or fr.get("begin") or fr.get("from")
                t = fr.get("end") or fr.get("to")
                if isinstance(s, int) and isinstance(t, int):
                    out.append((s, t))
        return out
    def _parse_cath(code):
        parts = code.split(".")
        if len(parts) >= 4:
            return parts[0], parts[1], parts[2], ".".join(parts[3:])
        return "", "", "", ""
    rows = []
    for e in entries or []:
        meta = _safe_meta(e)
        srcdb = (meta.get("source_database") or e.get("source_database") or "").strip().lower()
        acc = meta.get("accession") or meta.get("uid") or e.get("ac") or ""
        if srcdb in {"cathgene3d", "cath-gene3d"} or (isinstance(acc, str) and acc.startswith("G3DSA:")):
            rngs = _frags(e)
            cath_code = acc.split(":", 1)[1] if ":" in acc else acc
            c, a, t, h = _parse_cath(cath_code)
            rows.append({"struct_db": "CATH-Gene3D", "struct_id": acc, "cath_code": cath_code,
                         "cath_class": c, "cath_arch": a, "cath_topology": t,
                         "cath_superfamily": cath_code,
                         "ranges": ";".join([f"{s}..{t}" for s, t in rngs]),
                         "n_fragments": len(rngs),
                         "entry_type": (meta.get("type") or e.get("type") or "").strip(),
                         "evidence": "InterProAPI"})
    return rows


def positional_category(start, end, protein_length,
                        n_start_cut=0.2, c_end_cut=0.8):
    """
    Classify a domain occurrence positionally within its host protein.

    Thresholds (applied to 0-to-1 relative coordinates):
      Nterminal : relative_start < n_start_cut  (default 0.2)
                  → domain begins in the first 20 % of the sequence
      Cterminal : relative_end   > c_end_cut    (default 0.8)
                  → domain ends in the last 20 % of the sequence
      Internal  : all other cases
                  → domain lies fully within the middle 60 % of the sequence

    When both N- and C-terminal criteria are satisfied simultaneously
    (very short proteins or very long domains), Nterminal takes precedence.
    """
    if not protein_length or protein_length <= 0:
        return {"positional_category": "unknown", "relative_start": None,
                "relative_end": None, "relative_midpoint": None}
    rel_start = start / protein_length
    rel_end   = end   / protein_length
    rel_mid   = ((start + end) / 2.0) / protein_length
    if rel_start < n_start_cut:   cat = "Nterminal"
    elif rel_end  > c_end_cut:    cat = "Cterminal"
    else:                         cat = "Internal"
    return {"positional_category": cat, "relative_start": round(rel_start, 4),
            "relative_end": round(rel_end, 4), "relative_midpoint": round(rel_mid, 4)}


def build_domain_instances_from_interpro(dom_occ_rows, protein_length=None,
                                          assign_min_overlap_aa=10, assign_min_overlap_frac=0.2):
    instances: List[Dict[str, Any]] = []
    canon_by_ipr: Dict[str, List[Dict[str, Any]]] = {}
    for r in dom_occ_rows:
        if str(r.get("source_database", "")).strip().lower() != "interpro":
            continue
        ipr = r.get("ipr", "")
        if not ipr:
            continue
        segs = r.get("segments") or _parse_ranges_to_intervals(r.get("ranges", ""))
        segs = [(int(a), int(b)) for a, b in segs if isinstance(a, int) and isinstance(b, int) and a <= b]
        if not segs:
            continue
        s = min(a for a, _ in segs)
        e = max(b for _, b in segs)
        cont = "DiscontinuousDomain" if len(segs) > 1 else "ContinuousDomain"
        inst: Dict[str, Any] = {
            "start": s, "end": e, "segments": sorted(segs),
            "label": r.get("name") or ipr, "label_db": "INTERPRO", "label_acc": ipr,
            "signatures": set([f"INTERPRO:{ipr}"]),
            "go_terms": r.get("go_terms", ""), "go_terms_mf": r.get("go_terms_mf", ""),
            "go_terms_bp": r.get("go_terms_bp", ""), "go_terms_cc": r.get("go_terms_cc", ""),
            "function_tags": r.get("function_tags", ""),
            "evidence_level_domain": r.get("evidence_level", "InterProEntry"),
            "continuity_category": cont, "occurrence_index": r.get("occurrence_index", ""),
            "position_evidence": "FromInterProLocations",
        }
        if protein_length:
            inst.update(positional_category(s, e, protein_length))
        else:
            inst["positional_category"] = "unknown"
            inst["relative_start"] = inst["relative_end"] = inst["relative_midpoint"] = None
        instances.append(inst)
        canon_by_ipr.setdefault(ipr, []).append(inst)
    for r in dom_occ_rows:
        srcdb = str(r.get("source_database", "")).strip()
        if not srcdb or srcdb.lower() == "interpro":
            continue
        integrated = r.get("integrated", "")
        if not integrated:
            continue
        cands = canon_by_ipr.get(integrated, [])
        if not cands:
            continue
        segs = r.get("segments") or _parse_ranges_to_intervals(r.get("ranges", ""))
        segs = [(int(a), int(b)) for a, b in segs if isinstance(a, int) and isinstance(b, int) and a <= b]
        if not segs:
            continue
        sig_s = min(a for a, _ in segs)
        sig_e = max(b for _, b in segs)
        sig_len = sig_e - sig_s + 1
        best = None
        best_score = -1.0
        for inst in cands:
            ov = _interval_overlap_len((sig_s, sig_e), (inst["start"], inst["end"]))
            if ov <= 0:
                continue
            inst_len = inst["end"] - inst["start"] + 1
            frac = ov / float(min(sig_len, inst_len))
            score = frac * 1000.0 + ov
            if ov >= assign_min_overlap_aa and frac >= assign_min_overlap_frac and score > best_score:
                best_score = score
                best = inst
        if best is None:
            for inst in cands:
                if _interval_overlap_len((sig_s, sig_e), (inst["start"], inst["end"])) > 0:
                    best = inst; break
        if best is not None:
            token = f"{srcdb.upper()}:{r.get('ipr', '')}"
            if token.endswith(":") or token.endswith(":None"):
                continue
            best["signatures"].add(token)
    counts = Counter([inst["label_acc"] for inst in instances])
    for inst in instances:
        inst["signatures"] = sorted(inst["signatures"])
        inst["copy_number"] = counts[inst["label_acc"]]
    all_ids = sorted({f"INTERPRO:{inst['label_acc']}" for inst in instances})
    for inst in instances:
        this_id = f"INTERPRO:{inst['label_acc']}"
        inst["cooccurring_domain_types"] = ";".join([i for i in all_ids if i != this_id])
        inst["cooccurrence_evidence"] = "ComputedFromProteinDomainSet"
    return instances


def annotate_disorder_status_for_instances(domain_instances, features):
    disordered: List[Tuple[int, int]] = []
    for f in features or []:
        ftype = (f.get("type") or "").upper()
        if ftype in {"DISORDER", "COMPBIAS"}:
            s = _safe_int(f.get("start")); e = _safe_int(f.get("end"))
            if isinstance(s, int) and isinstance(e, int) and s <= e:
                disordered.append((s, e))
    for cl in domain_instances:
        s = _safe_int(cl.get("start")); e = _safe_int(cl.get("end"))
        status = "unknown"
        if isinstance(s, int) and isinstance(e, int) and e >= s:
            length = e - s + 1
            if length > 0 and disordered:
                overlap = sum(_interval_overlap_len((s, e), (ds, de)) for ds, de in disordered)
                frac = overlap / float(length)
                status = "disordered" if frac >= 0.7 else ("structured" if frac <= 0.3 else "mixed")
        cl["disorder_status"] = status
        cl["disorder_evidence"] = "FromUniProtFT" if disordered else "NoDisorderFeatures"


def annotate_topology_and_binding_for_instances(domain_instances, features,
        protein_subcellular_location="", protein_go_cc_ids=None,
        protein_reactome_ids=None, protein_kegg_ids=None):
    protein_go_cc_ids    = protein_go_cc_ids    or []
    protein_reactome_ids = protein_reactome_ids or []
    protein_kegg_ids     = protein_kegg_ids     or []
    topo: Dict[str, List[Tuple[int, int]]] = {
        "cytoplasmic": [], "extracellular": [], "luminal": [],
        "transmembrane": [], "signal_peptide": [],
    }
    binding_sites: List[Dict[str, Any]] = []
    for f in features or []:
        ftype = (f.get("type") or "").upper()
        s = _safe_int(f.get("start")); e = _safe_int(f.get("end"))
        if not isinstance(s, int) or not isinstance(e, int) or s > e:
            continue
        note = (f.get("note") or "").lower()
        if ftype == "TOPO_DOM":
            if "cytoplasm" in note or "cytoplasmic" in note:
                topo["cytoplasmic"].append((s, e))
            elif any(w in note for w in ["extracellular", "cell surface", "periplasmic"]):
                topo["extracellular"].append((s, e))
            elif any(w in note for w in ["lumenal", "lumen", "luminal"]):
                topo["luminal"].append((s, e))
        elif ftype == "TRANSMEM":
            topo["transmembrane"].append((s, e))
        elif ftype == "SIGNAL":
            topo["signal_peptide"].append((s, e))
        elif ftype in {"BINDING", "NP_BIND", "METAL"}:
            binding_sites.append({"type": ftype, "start": s, "end": e,
                                   "note": f.get("note", ""), "qualifiers": f.get("qualifiers", {}) or {}})
    def _overlaps_any(dom_iv, intervals):
        return any(_interval_overlap_len(dom_iv, itv) > 0 for itv in intervals)
    for cl in domain_instances:
        s = _safe_int(cl.get("start")); e = _safe_int(cl.get("end"))
        if not isinstance(s, int) or not isinstance(e, int) or s > e:
            cl["topology_context"] = "unknown"
            cl["topology_evidence"] = "NoTopologyFeatures"
            cl["binding_site_types"] = ""
            cl["binding_site_evidence"] = "NoBindingSites"
        else:
            dom_iv = (s, e)
            ctxs = {ctx for ctx, ivs in topo.items() if _overlaps_any(dom_iv, ivs)}
            if ctxs:
                cl["topology_context"] = ";".join(sorted(ctxs))
                cl["topology_evidence"] = "FromUniProtFT"
            else:
                cl["topology_context"] = "unknown"
                cl["topology_evidence"] = "NoTopologyFeatures"
            bs_tokens = [
                f"{site['type']}:{site['start']}..{site['end']}:{site.get('note','')}"
                for site in binding_sites
                if _interval_overlap_len(dom_iv, (site["start"], site["end"])) > 0
            ]
            cl["binding_site_types"] = ";".join(bs_tokens)
            cl["binding_site_evidence"] = "FromUniProtFT" if bs_tokens else "NoBindingSites"
        cl["protein_subcellular_location"] = protein_subcellular_location
        cl["protein_go_cc_ids"]    = ";".join(protein_go_cc_ids)
        cl["protein_reactome_ids"] = ";".join(protein_reactome_ids)
        cl["protein_kegg_ids"]     = ";".join(protein_kegg_ids)


def annotate_proximity_for_instances(domain_instances, features,
        near_tm_thresh=30, near_cleavage_thresh=20, near_active_thresh=10):
    tm_intervals  = _collect_intervals(features, {"TRANSMEM"})
    cleavage      = _signal_cleavage_pos(features)
    act_intervals = _collect_intervals(features, {"ACT_SITE", "BINDING", "METAL", "NP_BIND"})
    for cl in domain_instances:
        segs = cl.get("segments") or [(cl.get("start"), cl.get("end"))]
        segs = [(int(a), int(b)) for a, b in segs if isinstance(a, int) and isinstance(b, int) and a <= b]
        cats = []
        d_tm = None
        if tm_intervals and segs:
            d_tm = min(_segments_distance(segs, iv) for iv in tm_intervals)
            if d_tm is not None and d_tm <= near_tm_thresh:
                cats.append("NearTransmembraneHelix")
        cl["nearest_tm_distance"] = d_tm
        if cleavage is not None and segs:
            s = min(a for a, _ in segs); e = max(b for _, b in segs)
            d_clv = 0 if (s <= cleavage <= e) else min(abs(cleavage - s), abs(cleavage - e))
            cl["signal_cleavage_pos"] = cleavage
            cl["nearest_signal_cleavage_distance"] = d_clv
            if d_clv <= near_cleavage_thresh:
                cats.append("NearSignalPeptideCleavage")
        else:
            cl["signal_cleavage_pos"] = None
            cl["nearest_signal_cleavage_distance"] = None
        d_act = None
        if act_intervals and segs:
            d_act = min(_segments_distance(segs, iv) for iv in act_intervals)
            if d_act is not None and d_act <= near_active_thresh:
                cats.append("NearActiveSiteRegion")
        cl["nearest_active_site_distance"] = d_act
        cl["proximity_categories"] = ";".join(cats)
        cl["proximity_evidence"] = "ComputedFromUniProtFT"


def annotate_chebi_targets_for_instances(domain_instances, features):
    chebi_sites: List[Tuple[int, int, str]] = []
    for f in features or []:
        ftype = (f.get("type") or "").upper()
        if ftype not in {"BINDING", "NP_BIND", "METAL"}:
            continue
        s = _safe_int(f.get("start")); e = _safe_int(f.get("end"))
        if not isinstance(s, int) or not isinstance(e, int) or s > e:
            continue
        q = f.get("qualifiers") or {}
        ligand_id = q.get("ligand_id") or q.get("Ligand_id") or q.get("ligandId")
        if isinstance(ligand_id, list):
            for lid in ligand_id:
                if isinstance(lid, str) and "CHEBI:" in lid.upper():
                    chebi_sites.append((s, e, lid.replace("ChEBI:", "CHEBI:")))
        elif isinstance(ligand_id, str) and "CHEBI:" in ligand_id.upper():
            chebi_sites.append((s, e, ligand_id.replace("ChEBI:", "CHEBI:")))
    for cl in domain_instances:
        segs = cl.get("segments") or [(cl.get("start"), cl.get("end"))]
        segs = [(int(a), int(b)) for a, b in segs if isinstance(a, int) and isinstance(b, int) and a <= b]
        hits = {chebi for bs, be, chebi in chebi_sites
                if any(_interval_overlap_len(seg, (bs, be)) > 0 for seg in segs)}
        cl["binding_target_chebi"] = ";".join(sorted(hits))
        cl["binding_target_class_evidence"] = "UniProtFT:ligand_id" if hits else "None"


def _inst_id(d):
    return f"{(d.get('label_db') or '').upper()}:{d.get('label_acc') or ''}@{d.get('start')}-{d.get('end')}"

def _type_id(d):
    db = (d.get("label_db") or "").upper()
    acc = d.get("label_acc") or ""
    return f"{db}:{acc}" if (db and acc) else (d.get("label") or "")


def annotate_domain_order_plan_b(domain_instances):
    domain_instances.sort(key=lambda x: (_safe_int(x.get("start")) or 10**9,
                                          _safe_int(x.get("end"))   or 10**9))
    for i, cl in enumerate(domain_instances, 1):
        cl["ordinal_position_all"] = i
    groups: List[List[Dict[str, Any]]] = []
    cur_group: List[Dict[str, Any]] = []
    cur_end: Optional[int] = None
    for cl in domain_instances:
        s = _safe_int(cl.get("start")); e = _safe_int(cl.get("end"))
        if not isinstance(s, int) or not isinstance(e, int):
            continue
        if not cur_group:
            cur_group = [cl]; cur_end = e; continue
        if s <= (cur_end if cur_end is not None else -1):
            cur_group.append(cl); cur_end = max(cur_end, e) if cur_end is not None else e
        else:
            groups.append(cur_group); cur_group = [cl]; cur_end = e
    if cur_group:
        groups.append(cur_group)
    reps: List[Dict[str, Any]] = []
    rep_id_by_group: Dict[int, str] = {}
    for gi, g in enumerate(groups, 1):
        for cl in g:
            cl["overlap_group_id"] = gi; cl["overlap_group_size"] = len(g)
        def rep_key(x):
            s = _safe_int(x.get("start")) or 10**9; e = _safe_int(x.get("end")) or -1
            return (-(e - s + 1) if e >= s else 0, s, x.get("label_acc") or "")
        rep = sorted(g, key=rep_key)[0]
        reps.append(rep)
        rep_id = _inst_id(rep); rep_id_by_group[gi] = rep_id
        rep_s = _safe_int(rep.get("start")) or 10**9; rep_e = _safe_int(rep.get("end")) or -1
        for cl in g:
            cl["representative_id"] = rep_id
            cl["is_representative"] = "yes" if cl is rep else "no"
            s = _safe_int(cl.get("start")) or 10**9; e = _safe_int(cl.get("end")) or -1
            if cl is rep:
                cl["overlap_relation_to_rep"] = "self"
            elif rep_s <= s and e <= rep_e:
                cl["overlap_relation_to_rep"] = "contained_in_rep"
            elif s <= rep_s and rep_e <= e:
                cl["overlap_relation_to_rep"] = "contains_rep"
            else:
                cl["overlap_relation_to_rep"] = "overlaps_rep"
            cl["domain_order_evidence"] = "ComputedFromSortedCoordinates;PlanBLongestRepresentative"
    reps.sort(key=lambda x: (_safe_int(x.get("start")) or 10**9, _safe_int(x.get("end")) or 10**9))
    rep_ids = [_inst_id(r) for r in reps]
    rep_to_ord  = {rid: i + 1 for i, rid in enumerate(rep_ids)}
    rep_to_type = {rid: _type_id(r) for rid, r in zip(rep_ids, reps)}
    for cl in domain_instances:
        rid = cl.get("representative_id", "")
        cl["ordinal_position"]  = rep_to_ord.get(rid, "")
        cl["representative_type"] = rep_to_type.get(rid, "")
    prev_type: Dict[str, str] = {}; next_type: Dict[str, str] = {}
    prev_id:   Dict[str, str] = {}; next_id:   Dict[str, str] = {}
    for i, rid in enumerate(rep_ids):
        prev_id[rid]   = rep_ids[i - 1] if i > 0 else ""
        next_id[rid]   = rep_ids[i + 1] if i < len(rep_ids) - 1 else ""
        prev_type[rid] = rep_to_type.get(prev_id[rid], "") if prev_id[rid] else ""
        next_type[rid] = rep_to_type.get(next_id[rid], "") if next_id[rid] else ""
    for cl in domain_instances:
        rid = cl.get("representative_id", "")
        cl["prev_domain_type"] = prev_type.get(rid, ""); cl["next_domain_type"] = next_type.get(rid, "")
        cl["prev_representative_id"] = prev_id.get(rid, ""); cl["next_representative_id"] = next_id.get(rid, "")
        if cl.get("is_representative") == "yes" and rid in rep_ids:
            idx = rep_ids.index(rid); this_rep = reps[idx]
            if idx > 0:
                prev_rep = reps[idx - 1]
                gap = (_safe_int(this_rep["start"]) or 0) - (_safe_int(prev_rep["end"]) or 0)
                cl["gap_from_prev"] = gap
                cl["adjacent_to_prev"] = "yes" if isinstance(gap, int) and gap <= 1 else "no"
            else:
                cl["gap_from_prev"] = ""; cl["adjacent_to_prev"] = ""
            if idx < len(reps) - 1:
                nxt_rep = reps[idx + 1]
                gap2 = (_safe_int(nxt_rep["start"]) or 0) - (_safe_int(this_rep["end"]) or 0)
                cl["gap_to_next"] = gap2
                cl["adjacent_to_next"] = "yes" if isinstance(gap2, int) and gap2 <= 1 else "no"
            else:
                cl["gap_to_next"] = ""; cl["adjacent_to_next"] = ""
        else:
            cl["gap_from_prev"] = ""; cl["adjacent_to_prev"] = ""
            cl["gap_to_next"] = "";   cl["adjacent_to_next"] = ""


def annotate_insertion_discontinuity(domain_instances):
    for cl in domain_instances:
        segs = cl.get("segments") or []
        if len(segs) <= 1:
            cl["insertion_category"] = ""; continue
        segs_sorted = sorted(segs)
        inserted = False
        for (s1, e1), (s2, e2) in zip(segs_sorted, segs_sorted[1:]):
            if s2 <= e1 + 1: continue
            gap_iv = (e1 + 1, s2 - 1)
            for other in domain_instances:
                if other is cl: continue
                if _interval_overlap_len((other["start"], other["end"]), gap_iv) >= 10:
                    inserted = True; break
            if inserted: break
        cl["insertion_category"] = "InsertedDomain" if inserted else ""



def annotate_cooccurring_occurrences(domain_instances):
    """
    Populate the 'cooccurring_domain_occurrences' field on each domain instance.

    For every domain occurrence in the protein, this records the other domain
    occurrences that co-exist in the same protein — not just the type IDs
    (already stored in 'cooccurring_domain_types') but the specific occurrence
    with its coordinates and ordinal position.

    Format of each entry in the semicolon-delimited list:
        <type_id>|<start>..<end>|ord<ordinal>

    Example:
        INTERPRO:IPR003439|93..377|ord1;INTERPRO:IPR011527|386..633|ord2

    This column is useful for querying which specific occurrence of a
    multi-copy domain co-occurs with a given target — something that the
    type-only 'cooccurring_domain_types' column cannot answer.
    """
    for inst in domain_instances:
        this_id = "INTERPRO:" + str(inst.get("label_acc", ""))
        entries = []
        for other in domain_instances:
            if other is inst:
                continue
            other_id  = "INTERPRO:" + str(other.get("label_acc", ""))
            o_start   = other.get("start", "")
            o_end     = other.get("end", "")
            o_ord     = other.get("ordinal_position", "")
            coord_str = f"{o_start}..{o_end}" if (o_start != "" and o_end != "") else ""
            ord_str   = f"ord{o_ord}" if o_ord != "" else ""
            parts     = [p for p in [other_id, coord_str, ord_str] if p]
            entries.append("|".join(parts))
        inst["cooccurring_domain_occurrences"] = ";".join(entries)

def annotate_structural_context(domain_instances, structural_hits,
                                 min_overlap_aa=30, min_overlap_frac=0.5):
    hit_intervals = []
    for h in structural_hits:
        for (hs, he) in _parse_ranges_to_intervals(h.get("ranges", "")):
            hit_intervals.append((hs, he, h))
    for cl in domain_instances:
        s = _safe_int(cl.get("start")); e = _safe_int(cl.get("end"))
        if not isinstance(s, int) or not isinstance(e, int) or s > e: continue
        dom_iv = (s, e); dom_len = e - s + 1
        cath_sfs = set(); cath_debug = []
        for hs, he, h in hit_intervals:
            if h.get("struct_db") != "CATH-Gene3D": continue
            ov = _interval_overlap_len(dom_iv, (hs, he))
            if ov <= 0: continue
            frac = ov / float(min(dom_len, he - hs + 1))
            if ov >= min_overlap_aa and frac >= min_overlap_frac:
                sf = h.get("cath_superfamily", "")
                if sf: cath_sfs.add(sf); cath_debug.append(f"{sf}:{ov}:{frac:.2f}")
        cl["cath_superfamilies"]         = ";".join(sorted(cath_sfs))
        cl["cath_overlap_debug"]         = ";".join(cath_debug)
        cl["funfam_ids"]                 = ""
        cl["structural_coverage_quality"]= "predicted_gene3d" if cath_sfs else ""
        cl["structural_context_evidence"]= "InterProStructuralHitsOverlap" if cath_sfs else "None"


def summarize_domain_instances(instances):
    if not instances: return ""
    rows = []
    for i, cl in enumerate(sorted(instances, key=lambda x: (x.get("start", 10**9), x.get("end", 10**9))), 1):
        sigs = ", ".join(cl.get("signatures", []))
        cat  = cl.get("positional_category", "unknown")
        rows.append(f"{i}) {cl['start']}-{cl['end']} | {cl['label']} "
                    f"[{cl['label_db']}:{cl['label_acc']}] | pos={cat} | sigs={sigs}")
    return " ; ".join(rows)


def _build_interpro_summaries(dom_occ_rows):
    if not dom_occ_rows: return "", "", ""
    interpro_domains = [d for d in dom_occ_rows
                        if str(d.get("type","")).lower()=="domain"
                        and str(d.get("source_database","")).lower()=="interpro"]
    seen_names: set = set()
    ipr_names_list: List[str] = []
    for d in interpro_domains:
        nm = d.get("name", "")
        if nm and nm not in seen_names:
            seen_names.add(nm); ipr_names_list.append(nm)
    ipr_domain_names = "; ".join(ipr_names_list)
    ipr_rows = [f"{i}) {d.get('name','')} ({d.get('ipr','')}) | ranges={d.get('ranges','')}"
                for i, d in enumerate(interpro_domains, 1)]
    ipr_domains_pretty = " ; ".join(ipr_rows)
    sig_rows = []
    for d in dom_occ_rows:
        db = (d.get("source_database") or "").upper(); acc = d.get("ipr",""); nm = d.get("name","")
        label = f"{db}:{acc}" if db else acc
        _nm_part  = (' (' + nm + ')') if nm else ''
        _rng_part = (' | ranges=' + str(d.get('ranges', ''))) if d.get('ranges') else ''
        sig_rows.append(label + _nm_part + _rng_part)
    domain_signatures_pretty = " ; ".join(sig_rows)
    return ipr_domain_names, ipr_domains_pretty, domain_signatures_pretty


# Global state lists for run_uniprot_summary
_drops:    List[Tuple[str, str, str, str]] = []
_kept:     List[str] = []
_excluded: List[Tuple[str, str, str, str]] = []

def _drop(acc: str, entry_id: str, reason: str, detail: str = "") -> bool:
    _drops.append((acc or "", entry_id or "", reason, (detail or "")[:300]))
    return True


def run_uniprot_summary(
        query: str,
        reviewed: str = "yes",
        outdir: str = "uniprot_out",
        size: int = 300,
        sleep: float = DEFAULT_DELAY,
        exclude_regex: str = r"(hypothetical|putative|\blike\b)",
        fetch_interpro: bool = True,
        interpro_mode: str = "auto",
        save_interpro_json: bool = False,
) -> str:
    _drops.clear(); _kept.clear(); _excluded.clear()
    reviewed_flag: Optional[bool] = None
    rlow = (reviewed or "").lower().strip()
    if rlow in {"yes", "true", "reviewed"}:     reviewed_flag = True
    elif rlow in {"no", "false", "unreviewed"}: reviewed_flag = False
    Path(outdir).mkdir(parents=True, exist_ok=True)
    if save_interpro_json:
        Path(os.path.join(outdir, "interpro_json")).mkdir(parents=True, exist_ok=True)
    print(f"[INFO] Searching UniProtKB: {query} (reviewed={reviewed_flag})")
    items = search_uniprot(query, reviewed=reviewed_flag, size=size, delay=sleep)
    print(f"[INFO] Hits: {len(items)}")
    ex_re = re.compile(exclude_regex, flags=re.I)
    seen: set = set()
    cooccurrence_counter: Counter = Counter()
    repeat_counter: Counter = Counter()
    domain_type_labels: Dict[str, str] = {}
    summary_path = os.path.join(outdir, "proteins_summary.tsv")
    with open(summary_path, "w", newline="", encoding="utf-8") as fsum:
        w = csv.writer(fsum, delimiter="\t")
        w.writerow(["accession","entry_id","genes","protein_name","protein_length",
                    "go_ids_all","go_ids_mf","go_ids_bp","go_ids_cc",
                    "subcellular_location_text","reactome_ids","kegg_ids",
                    "region_count","compbias_count","domain_count","motif_count",
                    "regions_pretty","compbias_pretty","domains_pretty","motifs_pretty",
                    "function_text","interaction_text","cofactor_text","regulation_text",
                    "pfam_ids_protein","smart_ids_protein","prosite_ids_protein","cath_ids_protein",
                    "interpro_domain_instance_count","domain_instances_pretty",
                    "ipr_domain_names","ipr_domains_pretty","domain_signatures_pretty"])
        for it in items:
            acc      = it.get("primaryAccession") or it.get("accession") or ""
            entry_id = it.get("uniProtkbId") or it.get("id") or ""
            if not acc and _drop("", entry_id, "missing_accession", "no primaryAccession"): continue
            if acc in seen and _drop(acc, entry_id, "duplicate_accession", "already processed"): continue
            seen.add(acc)
            pname_search = extract_protein_name_string(it)
            m = ex_re.search(pname_search or "")
            if m:
                _excluded.append((acc, entry_id, pname_search or "", m.group(0)))
                _drop(acc, entry_id, "excluded_by_name", pname_search or "")
                print(f"[SKIP] {acc} excluded: {pname_search}")
                continue
            try:
                txt = fetch_uniprot_txt(acc, delay=sleep)
            except Exception as e:
                _drop(acc, entry_id, "fetch_failed", str(e))
                print(f"[WARN] fetch failed for {acc}: {e}")
                continue
            parsed = parse_uniprot_txt(txt)
            genes       = parsed.get("genes_from_gn") or ""
            pname       = parsed.get("protein_name_from_txt") or pname_search
            entry_id_txt = parsed.get("entry_id_from_txt")
            if entry_id_txt and not entry_id: entry_id = entry_id_txt
            pfam    = ";".join(parsed["cross_refs"].get("Pfam", []))
            smart   = ";".join(parsed["cross_refs"].get("SMART", []))
            prosite = ";".join(parsed["cross_refs"].get("PROSITE", []))
            cath    = ";".join(parsed["cross_refs"].get("CATH", []))
            go_ids    = ";".join(parsed.get("go_ids", []))
            go_mf_ids = ";".join(parsed.get("go_mf_ids", []))
            go_bp_ids = ";".join(parsed.get("go_bp_ids", []))
            go_cc_ids = ";".join(parsed.get("go_cc_ids", []))
            subcell_loc  = parsed.get("cc_subcellular_location", "")
            reactome_ids = ";".join(parsed.get("cross_refs", {}).get("Reactome", []))
            kegg_ids     = ";".join(parsed.get("cross_refs", {}).get("KEGG", []))
            seq      = parsed.get("sequence", "")
            prot_len = parsed.get("protein_length") or (len(seq) if seq else None)
            if prot_len is None: prot_len = _length_from_search_item(it) or 0
            regions  = parsed.get("regions",  [])
            compbias = parsed.get("compbias", [])
            domains  = parsed.get("domains",  [])
            motifs   = parsed.get("motifs",   [])
            regions_pretty  = feature_pretty_list(regions,  seq, max_items=12, full_seq=False)
            compbias_pretty = feature_pretty_list(compbias, seq, max_items=12, full_seq=False)
            domains_pretty  = feature_pretty_list(domains,  seq, max_items=12, full_seq=False)
            motifs_pretty   = feature_pretty_list(motifs,   seq, max_items=12, full_seq=False)
            feat_path = os.path.join(outdir, f"features_{acc}.tsv")
            with open(feat_path, "w", newline="", encoding="utf-8") as ff:
                wf = csv.writer(ff, delimiter="\t")
                wf.writerow(["type","ftid","start","end","length","description","qualifiers_json","sequence"])
                for f in parsed.get("features", []):
                    if f.get("type") not in {"REGION","COMPBIAS","DOMAIN","MOTIF","TRANSMEM","SIGNAL",
                                             "TOPO_DOM","BINDING","NP_BIND","METAL","ACT_SITE"}:
                        continue
                    s_i, t_i = _safe_int(f.get("start")), _safe_int(f.get("end"))
                    subseq = ""
                    if seq and isinstance(s_i,int) and isinstance(t_i,int) and 1<=s_i<=t_i<=len(seq):
                        subseq = seq[s_i-1:t_i]
                    length = len(subseq) if subseq else (t_i-s_i+1 if s_i and t_i else "")
                    wf.writerow([f.get("type",""),f.get("ftid",""),s_i or "",t_i or "",length,
                                 f.get("note",""),json.dumps(f.get("qualifiers") or {}, ensure_ascii=False),subseq])
            write_features_pretty_txt(outdir, acc, parsed.get("features", []), seq)
            interpro_instance_count = 0; domain_instances_pretty = ""
            ipr_domain_names = ""; ipr_domains_pretty = ""; domain_signatures_pretty = ""
            if fetch_interpro:
                mode_eff = interpro_mode
                if mode_eff == "auto":
                    rev = _is_reviewed_search_item(it)
                    if rev is True:     mode_eff = "reviewed"
                    elif rev is False:  mode_eff = "unreviewed"
                    else:               mode_eff = "uniprot"
                ipr_entries = interpro_entries_for_protein(acc, mode=mode_eff, page_size=200, delay=0.5)
                if not ipr_entries and mode_eff not in ("uniprot", "both"):
                    ipr_entries = interpro_entries_for_protein(acc, mode="uniprot", page_size=200, delay=0.5)
                dom_occ_rows, member_rows = interpro_domain_occurrences(ipr_entries, uniprot_acc=acc, protein_seq=seq)
                struct_hits = interpro_structural_hits(ipr_entries, uniprot_acc=acc)
                domain_instances = build_domain_instances_from_interpro(dom_occ_rows, protein_length=prot_len)
                domain_type_counts: Counter = Counter()
                for cl in domain_instances:
                    dt_id = f"{(cl.get('label_db') or '').upper()}:{cl.get('label_acc')}"
                    domain_type_counts[dt_id] += 1
                    if dt_id not in domain_type_labels:
                        domain_type_labels[dt_id] = cl.get("label", "")
                dt_list = sorted(domain_type_counts.keys())
                for dt_id, cnt in domain_type_counts.items():
                    if cnt > 1: repeat_counter[dt_id] += 1
                for i1 in range(len(dt_list)):
                    for i2 in range(i1+1, len(dt_list)):
                        cooccurrence_counter[(dt_list[i1], dt_list[i2])] += 1
                annotate_disorder_status_for_instances(domain_instances, parsed.get("features", []))
                annotate_topology_and_binding_for_instances(
                    domain_instances, parsed.get("features", []),
                    protein_subcellular_location=parsed.get("cc_subcellular_location",""),
                    protein_go_cc_ids=parsed.get("go_cc_ids",[]),
                    protein_reactome_ids=parsed.get("cross_refs",{}).get("Reactome",[]),
                    protein_kegg_ids=parsed.get("cross_refs",{}).get("KEGG",[]),
                )
                annotate_proximity_for_instances(domain_instances, parsed.get("features", []))
                annotate_chebi_targets_for_instances(domain_instances, parsed.get("features", []))
                annotate_domain_order_plan_b(domain_instances)
                annotate_insertion_discontinuity(domain_instances)
                annotate_cooccurring_occurrences(domain_instances)  # must run after Plan B assigns ordinals
                annotate_structural_context(domain_instances, struct_hits, min_overlap_aa=30, min_overlap_frac=0.5)
                for cl in domain_instances:
                    cl["domain_evidence_level"] = cl.get("evidence_level_domain") or "InterProEntry"
                interpro_instance_count = len(domain_instances)
                domain_instances_pretty = summarize_domain_instances(domain_instances)
                ipr_domain_names, ipr_domains_pretty, domain_signatures_pretty = _build_interpro_summaries(dom_occ_rows)
                ipr_path = os.path.join(outdir, f"interpro_domains_{acc}.tsv")
                with open(ipr_path, "w", newline="", encoding="utf-8") as fi:
                    wi = csv.writer(fi, delimiter="\t")
                    wi.writerow(["entry_accession","name","source_database","integrated_ipr",
                                 "occurrence_index","start","end","n_fragments","ranges",
                                 "go_terms_all","go_terms_mf","go_terms_bp","go_terms_cc",
                                 "function_tags","evidence_level"])
                    for d in dom_occ_rows:
                        wi.writerow([d.get("ipr",""),d.get("name",""),d.get("source_database",""),
                                     d.get("integrated",""),d.get("occurrence_index",""),
                                     d.get("start",""),d.get("end",""),
                                     d.get("n_fragments",0),d.get("ranges",""),
                                     d.get("go_terms",""),d.get("go_terms_mf",""),d.get("go_terms_bp",""),
                                     d.get("go_terms_cc",""),d.get("function_tags",""),d.get("evidence_level","")])
                if member_rows:
                    mem_path = os.path.join(outdir, f"interpro_domain_members_{acc}.tsv")
                    with open(mem_path, "w", newline="", encoding="utf-8") as fm:
                        wm = csv.writer(fm, delimiter="\t")
                        wm.writerow(["ipr","interpro_name","member_db","member_accession",
                                     "member_name","member_source_database","member_type"])
                        for r in member_rows:
                            wm.writerow([r["ipr"],r["interpro_name"],r["member_db"],r["member_accession"],
                                         r["member_name"],r["member_source_database"],r["member_type"]])
                inst_path = os.path.join(outdir, f"domain_instances_{acc}.tsv")
                with open(inst_path, "w", newline="", encoding="utf-8") as fd:
                    wd = csv.writer(fd, delimiter="\t")
                    wd.writerow(["protein_accession","start","end",
                                 "start_position","end_position",
                                 "start_residue","end_residue","domain_sequence",
                                 "segments","label","label_db","label_acc",
                                 "ordinal_position","ordinal_position_all","overlap_group_id","overlap_group_size",
                                 "is_representative","representative_id","representative_type","overlap_relation_to_rep",
                                 "prev_domain_type","gap_from_prev","adjacent_to_prev",
                                 "next_domain_type","gap_to_next","adjacent_to_next",
                                 "prev_representative_id","next_representative_id",
                                 "positional_category","relative_start","relative_end","relative_midpoint",
                                 "continuity_category","insertion_category","copy_number","signatures",
                                 "go_terms_all","go_terms_mf","go_terms_bp","go_terms_cc","function_tags",
                                 "disorder_status","topology_context","proximity_categories",
                                 "nearest_tm_distance","signal_cleavage_pos","nearest_signal_cleavage_distance",
                                 "nearest_active_site_distance","binding_site_types","binding_target_chebi",
                                 "cooccurring_domain_types","cooccurring_domain_occurrences",
                                 "cath_superfamilies","cath_overlap_debug",
                                 "funfam_ids","structural_coverage_quality","position_evidence",
                                 "cooccurrence_evidence","domain_order_evidence","disorder_evidence",
                                 "topology_evidence","proximity_evidence","binding_site_evidence",
                                 "binding_target_class_evidence","structural_context_evidence",
                                 "domain_evidence_level","protein_subcellular_location",
                                 "protein_go_cc_ids",
                                 "function_tags","binding_partner","binding_partner_specific",
                                 "note","reference"])
                    for cl in domain_instances:
                        # ── five sequence columns computed here in Stage 1 ──────────
                        # seq is already in memory from the UniProt flat-file parse
                        # above (line: seq = parsed.get("sequence", "")).
                        # No extra API call needed — this is why Stage 2.5 is removed.
                        _s = cl.get("start")
                        _e = cl.get("end")
                        try:
                            _si = int(float(str(_s))) if _s not in (None, "") else 0
                        except (ValueError, TypeError):
                            _si = 0
                        try:
                            _ei = int(float(str(_e))) if _e not in (None, "") else 0
                        except (ValueError, TypeError):
                            _ei = 0
                        _seq_ok = (seq and _si >= 1 and _ei >= _si
                                   and _ei <= len(seq))
                        _start_res  = seq[_si - 1]      if _seq_ok else ""
                        _end_res    = seq[_ei - 1]      if _seq_ok else ""
                        _dom_seq    = seq[_si - 1: _ei] if _seq_ok else ""
                        # ────────────────────────────────────────────────────────────
                        wd.writerow([acc,
                            cl.get("start"), cl.get("end"),
                            _si or "", _ei or "",
                            _start_res, _end_res, _dom_seq,
                            ";".join([str(a)+".."+str(b) for a,b in (cl.get("segments") or [])]),
                            cl.get("label"),cl.get("label_db"),cl.get("label_acc"),
                            cl.get("ordinal_position",""),cl.get("ordinal_position_all",""),
                            cl.get("overlap_group_id",""),cl.get("overlap_group_size",""),
                            cl.get("is_representative",""),cl.get("representative_id",""),
                            cl.get("representative_type",""),cl.get("overlap_relation_to_rep",""),
                            cl.get("prev_domain_type",""),cl.get("gap_from_prev",""),cl.get("adjacent_to_prev",""),
                            cl.get("next_domain_type",""),cl.get("gap_to_next",""),cl.get("adjacent_to_next",""),
                            cl.get("prev_representative_id",""),cl.get("next_representative_id",""),
                            cl.get("positional_category",""),cl.get("relative_start",""),
                            cl.get("relative_end",""),cl.get("relative_midpoint",""),
                            cl.get("continuity_category",""),cl.get("insertion_category",""),
                            cl.get("copy_number"),",".join(cl.get("signatures",[])),
                            cl.get("go_terms",""),cl.get("go_terms_mf",""),cl.get("go_terms_bp",""),
                            cl.get("go_terms_cc",""),cl.get("function_tags",""),
                            cl.get("disorder_status",""),cl.get("topology_context",""),
                            cl.get("proximity_categories",""),cl.get("nearest_tm_distance",""),
                            cl.get("signal_cleavage_pos",""),cl.get("nearest_signal_cleavage_distance",""),
                            cl.get("nearest_active_site_distance",""),cl.get("binding_site_types",""),
                            cl.get("binding_target_chebi",""),
                            cl.get("cooccurring_domain_types",""),
                            cl.get("cooccurring_domain_occurrences",""),
                            cl.get("cath_superfamilies",""),cl.get("cath_overlap_debug",""),
                            cl.get("funfam_ids",""),cl.get("structural_coverage_quality",""),
                            cl.get("position_evidence",""),cl.get("cooccurrence_evidence",""),
                            cl.get("domain_order_evidence",""),cl.get("disorder_evidence",""),
                            cl.get("topology_evidence",""),cl.get("proximity_evidence",""),
                            cl.get("binding_site_evidence",""),cl.get("binding_target_class_evidence",""),
                            cl.get("structural_context_evidence",""),cl.get("domain_evidence_level",""),
                            cl.get("protein_subcellular_location",""),cl.get("protein_go_cc_ids",""),
                            "",  # function_tags  — fill manually in Stage 3
                            "",  # binding_partner
                            "",  # binding_partner_specific
                            "",  # note
                            "",  # reference
                            ])
            interaction_text = _norm_space(
                ((parsed.get("cc_subunit","") + " " + parsed.get("cc_interaction","")).strip()))
            w.writerow([acc,entry_id,genes,pname,prot_len,
                        go_ids,go_mf_ids,go_bp_ids,go_cc_ids,subcell_loc,reactome_ids,kegg_ids,
                        len(regions),len(compbias),len(domains),len(motifs),
                        regions_pretty,compbias_pretty,domains_pretty,motifs_pretty,
                        parsed.get("cc_function",""),interaction_text,
                        parsed.get("cc_cofactor",""),parsed.get("cc_regulation",""),
                        pfam,smart,prosite,cath,interpro_instance_count,domain_instances_pretty,
                        ipr_domain_names,ipr_domains_pretty,domain_signatures_pretty])
            _kept.append(acc)
    # write auxiliary log files
    if _drops:
        drops_path = os.path.join(outdir, "dropped_accessions.tsv")
        with open(drops_path, "w", newline="", encoding="utf-8") as fd:
            wr = csv.writer(fd, delimiter="\t")
            wr.writerow(["accession","entry_id","reason","detail"])
            wr.writerows(_drops)
    if _kept:
        kept_path = os.path.join(outdir, "kept_accessions.tsv")
        with open(kept_path, "w", newline="", encoding="utf-8") as fk:
            wr = csv.writer(fk, delimiter="\t")
            wr.writerow(["accession"])
            wr.writerows([[a] for a in _kept])
    if cooccurrence_counter:
        co_path = os.path.join(outdir, "domain_type_cooccurrence.tsv")
        with open(co_path, "w", newline="", encoding="utf-8") as fc:
            wc = csv.writer(fc, delimiter="\t")
            wc.writerow(["domain_type_A","domain_label_A","domain_type_B","domain_label_B","n_proteins_cooccurring"])
            for (dtA, dtB), cnt in sorted(cooccurrence_counter.items(), key=lambda x: -x[1]):
                wc.writerow([dtA,domain_type_labels.get(dtA,""),dtB,domain_type_labels.get(dtB,""),cnt])
    print(f"[INFO] Done: {len(_kept)} proteins processed, {len(_drops)} dropped.")
    return outdir


# ═══════════════════════════════════════════════════════════════════════════════
# PART 4 — Ontology export builder
# ═══════════════════════════════════════════════════════════════════════════════

# ── Fully inlined from new_domain_integrated.py ───────────────────────────────

def _ndi_split_ids(cell: Any) -> List[str]:
    if not isinstance(cell, str): return []
    return [p2 for p2 in (p.strip() for p in re.split(r"[;,]", cell)) if p2]

def ensure_columns(df: "pd.DataFrame", columns: List[str]) -> "pd.DataFrame":
    if df is None or (df.empty and len(df.columns) == 0):
        return pd.DataFrame(columns=columns)
    out = df.copy()
    for c in columns:
        if c not in out.columns: out[c] = ""
    return out

def empty_domain_occurrence_columns() -> List[str]:
    return [
        # ── Identity ────────────────────────────────────────────────────
        "occurrence_id", "protein_accession", "domain_type_id",
        "label", "label_db", "label_acc",
        # ── Coordinates ─────────────────────────────────────────────────
        "start", "end", "segments",
        "start_position", "end_position",      # 1-based ints, mirrors start/end
        "start_residue", "end_residue",         # single amino acid at boundary
        "domain_sequence",                      # full aa sequence of domain
        # ── Positional context ──────────────────────────────────────────
        "positional_category",                  # Nterminal / Internal / Cterminal
        "relative_start", "relative_end", "relative_midpoint",
        # ── Ordinal / copy-number ────────────────────────────────────────
        "ordinal_position_group", "ordinal_position_all", "ordinal_label_nc",
        "copy_number",
        # ── Continuity / overlap ────────────────────────────────────────
        "continuity_category", "insertion_category",
        "is_representative", "representative_id", "representative_type",
        "overlap_group_id", "overlap_group_size", "overlap_relation_to_rep",
        # ── Domain order ────────────────────────────────────────────────
        "prev_domain_type", "next_domain_type",
        "gap_from_prev", "gap_to_next",
        "adjacent_to_prev", "adjacent_to_next",
        # ── Co-occurrence ────────────────────────────────────────────────
        "cooccurring_domain_types",
        "cooccurring_domain_occurrences",       # type|start..end|ordinal for each co-occurring occ
        # ── Signatures / annotations ────────────────────────────────────
        "signatures",
        "disorder_status", "topology_context", "proximity_categories",
        "nearest_tm_distance", "signal_cleavage_pos",
        "nearest_signal_cleavage_distance", "nearest_active_site_distance",
        "go_terms_all", "go_terms_mf", "go_terms_bp", "go_terms_cc",
        "binding_site_types", "binding_target_chebi",
        "cath_superfamilies", "structural_coverage_quality",
        # ── Evidence provenance ──────────────────────────────────────────
        "position_evidence", "domain_order_evidence", "cooccurrence_evidence",
        "disorder_evidence", "topology_evidence", "proximity_evidence",
        "binding_site_evidence", "binding_target_class_evidence",
        "structural_context_evidence", "domain_evidence_level",
        "protein_subcellular_location", "protein_go_cc_ids",
        # ── Manual annotation columns (Stage 3 — fill in manually) ──────
        "function_tags",            # DOF vocabulary terms
        "binding_partner",          # DOT vocabulary terms (general)
        "binding_partner_specific", # DOT vocabulary terms (specific molecule)
        "note",                     # free-text curation note
        "reference",                # PMID or DOI supporting the annotation
    ]

def empty_domain_order_edge_columns() -> List[str]:
    return ["protein_accession","occurrence_id","prev_occurrence_id","next_occurrence_id",
            "gap_from_prev","gap_to_next","adjacent_to_prev","adjacent_to_next",
            "edge_scope","domain_order_evidence"]

def empty_domain_type_cooccurrence_columns() -> List[str]:
    return ["domain_type_id_A","domain_type_name_A","domain_type_id_B","domain_type_name_B",
            "supportCount_proteins","proteins_with_pair","association_type","association_evidence"]

def empty_occurrence_structural_link_columns() -> List[str]:
    return ["occurrence_id","structural_classification_id","structural_db","link_type","evidence"]

def empty_occurrence_binding_link_columns() -> List[str]:
    return ["occurrence_id","binding_target_id","binding_db","link_type","evidence"]

def empty_occurrence_signature_link_columns() -> List[str]:
    return ["occurrence_id","signature_id","signature_db","link_type","evidence"]

def load_proteins_summary(proteins_path: Path) -> "pd.DataFrame":
    df = pd.read_csv(proteins_path, sep="\t", dtype=str, encoding="utf-8").fillna("")
    if "accession" not in df.columns:
        for alt in ("acc","Entry","entry","uniprot","primary_accession"):
            if alt in df.columns: df = df.rename(columns={alt:"accession"}); break
    return df

def build_protein_table(df: "pd.DataFrame") -> "pd.DataFrame":
    cols = {"protein_accession":"accession","entry_id":"entry_id","genes":"genes",
            "protein_name":"protein_name","protein_length":"protein_length",
            "protein_go_ids_all":"go_ids_all","protein_go_ids_mf":"go_ids_mf",
            "protein_go_ids_bp":"go_ids_bp","protein_go_ids_cc":"go_ids_cc",
            "subcellular_location_text":"subcellular_location_text",
            "reactome_ids":"reactome_ids","kegg_ids":"kegg_ids"}
    out = {oc: (df[ic] if ic in df.columns else "") for oc, ic in cols.items()}
    proteins = pd.DataFrame(out).fillna("")
    proteins["protein_accession"] = proteins["protein_accession"].astype(str).str.strip()
    return proteins[proteins["protein_accession"] != ""].drop_duplicates("protein_accession")

def parse_regions_pretty_block(block: str) -> List[Dict[str, Any]]:
    if not isinstance(block, str) or not block.strip(): return []
    items = REGION_ITEM_SPLIT.split(block)
    out = []
    for it in items:
        it = it.strip()
        if not it: continue
        m_range = RANGE_RE.search(it)
        if not m_range: continue
        start, end = int(m_range.group(1)), int(m_range.group(2))
        m_len  = LEN_RE.search(it);  m_desc = DESC_RE.search(it)
        m_seq  = SEQ_RE.search(it);  m_name = NAME_RE.search(it)
        out.append({"start": start, "end": end,
                    "length": int(m_len.group(1)) if m_len else (end - start + 1),
                    "desc": (m_desc.group(1).strip() if m_desc else ""),
                    "seq":  (m_seq.group(1).strip()  if m_seq  else ""),
                    "name": (m_name.group(1).strip() if m_name else "")})
    return out

def render_regions_pretty(regs):
    parts = []
    for i, r in enumerate(regs, 1):
        fields = [f"{r['start']}-{r['end']}", f"len={r.get('length', r['end']-r['start']+1)}"]
        if r.get("seq"):  fields.append(f"seq={r['seq']}")
        if r.get("name"): fields.append(f"name={r['name']}")
        fields.append(f"desc={r.get('desc','')}")
        parts.append(f"{i}) " + " | ".join(fields))
    return " ; ".join(parts)

def safe_parse_list(value: Any) -> List[dict]:
    if not isinstance(value, str): return []
    txt = value.strip()
    if not txt: return []
    if txt.startswith("["):
        for parser in (json.loads, ast.literal_eval):
            try:
                parsed = parser(txt)
                if isinstance(parsed, list): return parsed
            except Exception: pass
    return []

def filter_disordered_regions(df: "pd.DataFrame") -> "pd.DataFrame":
    df = df.copy()
    if "regions_pretty" not in df.columns:
        df["regions_filtered_pretty"] = ""; df["region_remained_count"] = ""; return df
    filtered_pretty = []; region_counts = []
    for _, row in df.iterrows():
        raw = row.get("regions_pretty", "")
        regs = parse_regions_pretty_block(raw)
        if not regs:
            lst = safe_parse_list(raw); regs = []
            for r in lst:
                try:
                    start = int(r.get("start") or r.get("begin") or r.get("from"))
                    end   = int(r.get("end") or r.get("to"))
                    regs.append({"start": start, "end": end, "length": int(r.get("length") or end-start+1),
                                 "desc": (r.get("desc") or r.get("description") or r.get("name") or "").strip(),
                                 "seq": r.get("seq",""), "name": r.get("name","")})
                except Exception: continue
        kept = [r for r in regs if "disordered" not in (r.get("desc","").lower())]
        filtered_pretty.append(render_regions_pretty(kept)); region_counts.append(len(kept))
    df["regions_filtered_pretty"] = filtered_pretty; df["region_remained_count"] = region_counts
    return df

def normalize_db_label(raw: str) -> str:
    """Normalise a member-database label to a canonical uppercase short form.
    e.g. 'pfam' -> 'PFAM', 'smart' -> 'SMART', 'gene3d' -> 'CATH-Gene3D'
    Unknown values are returned as uppercase stripped strings.
    """
    _MAP = {
        "pfam":          "PFAM",
        "smart":         "SMART",
        "cdd":           "CDD",
        "prints":        "PRINTS",
        "panther":       "PANTHER",
        "tigrfam":       "TIGRFAM",
        "hamap":         "HAMAP",
        "pirsf":         "PIRSF",
        "superfamily":   "SUPERFAMILY",
        "gene3d":        "CATH-Gene3D",
        "cath-gene3d":   "CATH-Gene3D",
        "cathgene3d":    "CATH-Gene3D",
        "ssf":           "SUPERFAMILY",
        "prosite":       "PROSITE",
        "prosite profiles": "PROSITE",
        "prosite patterns": "PROSITE",
        "ncbifam":       "NCBIfam",
        "sfld":          "SFLD",
        "funfam":        "FunFam",
        "profile":       "PROFILE",
        "mobidblt":      "MobiDB-lite",
        "mobidb":        "MobiDB-lite",
        "interpro":      "INTERPRO",
        "":              "",
    }
    key = str(raw).strip().lower()
    return _MAP.get(key, str(raw).strip().upper())


def load_domain_instances(domain_dir: Path) -> "pd.DataFrame":
    files = sorted(domain_dir.glob("domain_instances_*.tsv"))
    if not files: return pd.DataFrame()
    frames = []
    for f in files:
        dfi = pd.read_csv(f, sep="\t", dtype=str, encoding="utf-8").fillna("")
        if "protein_accession" not in dfi.columns:
            dfi["protein_accession"] = f.stem.replace("domain_instances_","")
        frames.append(dfi)
    df_occ = pd.concat(frames, ignore_index=True).fillna("")
    df_occ["protein_accession"] = df_occ["protein_accession"].astype(str).str.strip()
    return df_occ[df_occ["protein_accession"] != ""]

def _ensure_numeric_coords(df: "pd.DataFrame") -> "pd.DataFrame":
    df = df.copy()
    df["start_i"] = pd.to_numeric(df.get("start",""), errors="coerce")
    df["end_i"]   = pd.to_numeric(df.get("end",""),   errors="coerce")
    df = df[df["start_i"].notna() & df["end_i"].notna()].copy()
    df["start_i"] = df["start_i"].astype(int); df["end_i"] = df["end_i"].astype(int)
    return df

def _ensure_domain_type_id(df: "pd.DataFrame") -> "pd.DataFrame":
    df = df.copy()
    if "domain_type_id" not in df.columns or (df["domain_type_id"].astype(str).str.strip()=="").all():
        db  = df.get("label_db","").astype(str).map(normalize_db_label)
        acc = df.get("label_acc","").astype(str).str.strip()
        df["domain_type_id"] = (db+":"+acc).str.replace(r"^:$","",regex=True)
    df["domain_type_id"] = df["domain_type_id"].astype(str).str.strip()
    return df[df["domain_type_id"] != ""].copy()

def assign_occurrence_ids_preserve_plan_b(df_occ_raw: "pd.DataFrame") -> Tuple["pd.DataFrame","pd.DataFrame"]:
    if df_occ_raw.empty:
        return df_occ_raw.assign(occurrence_id="",domain_type_id="",
                                  ordinal_position_group="",ordinal_label_nc=""), \
               pd.DataFrame(columns=empty_domain_order_edge_columns())
    df = df_occ_raw.copy().fillna("")
    df = _ensure_numeric_coords(df); df = _ensure_domain_type_id(df)
    if "positional_category" in df.columns:
        df["positional_category"] = df["positional_category"].astype(str).map(normalize_positional_category)
    else:
        df["positional_category"] = "unknown"
    df["ordinal_position_group"] = df["ordinal_position"].astype(str) if "ordinal_position" in df.columns else ""
    if "ordinal_position_all" not in df.columns: df["ordinal_position_all"] = ""
    occ_ids = []
    for prot, sub in df.groupby("protein_accession", sort=False):
        sub2 = sub.sort_values(["start_i","end_i","domain_type_id"]).copy()
        if (sub2["ordinal_position_all"].astype(str).str.strip()!="").any():
            conv = pd.to_numeric(sub2["ordinal_position_all"], errors="coerce")
            local_idx = conv.astype(int).tolist() if conv.notna().all() else list(range(1,len(sub2)+1))
        else:
            local_idx = list(range(1,len(sub2)+1))
        for (idx, row), i_local in zip(sub2.iterrows(), local_idx):
            s = int(row["start_i"]); e = int(row["end_i"])
            occ_ids.append((idx, f"{prot}|{row['domain_type_id']}|{s}{OCC_RANGE_DELIM}{e}|{i_local}"))
    df["occurrence_id"] = df.index.map(dict(occ_ids).get)
    has_plan_b = all(c in df.columns for c in ["is_representative","representative_id",
                                                 "prev_representative_id","next_representative_id"])
    df["ordinal_label_nc"] = ""
    edges: List[Dict[str, Any]] = []
    for prot, sub in df.groupby("protein_accession", sort=False):
        sub2 = sub.sort_values(["start_i","end_i","domain_type_id"]).copy()
        if has_plan_b and (sub2["is_representative"].astype(str).str.lower()=="yes").any():
            reps = sub2[sub2["is_representative"].astype(str).str.lower()=="yes"].copy()
            rep_map = {str(r.get("representative_id","")).strip(): str(r.get("occurrence_id","")).strip()
                       for _, r in reps.iterrows() if str(r.get("representative_id","")).strip()}
            reps = reps.sort_values(["start_i","end_i","domain_type_id"]).copy()
            rep_occ_ids = reps["occurrence_id"].astype(str).tolist()
            nrep = len(rep_occ_ids)
            rep_labels = (["N"] if nrep==1 else
                          ["N"] + [f"N+{i}" for i in range(1,nrep-1)] + ["C"])
            rep_occ_to_label = dict(zip(rep_occ_ids, rep_labels))
            repid_to_repocc  = {str(r.get("representative_id","")).strip(): str(r.get("occurrence_id","")).strip()
                                 for _,r in reps.iterrows()}
            for idx, r in sub2.iterrows():
                rid = str(r.get("representative_id","")).strip()
                rep_occ = repid_to_repocc.get(rid,"")
                df.at[idx,"ordinal_label_nc"] = rep_occ_to_label.get(rep_occ,"")
            for _, r in reps.iterrows():
                cur_occ  = str(r.get("occurrence_id","")).strip()
                prev_rid = str(r.get("prev_representative_id","")).strip()
                next_rid = str(r.get("next_representative_id","")).strip()
                edges.append({"protein_accession": prot, "occurrence_id": cur_occ,
                               "prev_occurrence_id": rep_map.get(prev_rid,"") if prev_rid else "",
                               "next_occurrence_id": rep_map.get(next_rid,"") if next_rid else "",
                               "gap_from_prev": str(r.get("gap_from_prev","")).strip(),
                               "gap_to_next":   str(r.get("gap_to_next","")).strip(),
                               "adjacent_to_prev": str(r.get("adjacent_to_prev","")).strip(),
                               "adjacent_to_next": str(r.get("adjacent_to_next","")).strip(),
                               "edge_scope": "representatives",
                               "domain_order_evidence": str(r.get("domain_order_evidence","")).strip()
                                   or "ComputedFromSortedCoordinates;PlanBLongestRepresentative"})
        else:
            n = len(sub2)
            labels = (["N"] if n==1 else ["N"]+[f"N+{i}" for i in range(1,n-1)]+["C"])
            for (idx, r), lbl in zip(sub2.iterrows(), labels):
                df.at[idx,"ordinal_label_nc"] = lbl
            occs = sub2["occurrence_id"].astype(str).tolist()
            for i, (_, r) in enumerate(sub2.iterrows()):
                prev_gap = next_gap = adj_prev = adj_next = ""
                if i > 0:
                    prev_row = sub2.iloc[i-1]
                    pv = int(r["start_i"]) - int(prev_row["end_i"])
                    prev_gap = str(pv); adj_prev = "yes" if pv<=1 else "no"
                if i < n-1:
                    nxt_row = sub2.iloc[i+1]
                    nv = int(nxt_row["start_i"]) - int(r["end_i"])
                    next_gap = str(nv); adj_next = "yes" if nv<=1 else "no"
                edges.append({"protein_accession": prot, "occurrence_id": occs[i],
                               "prev_occurrence_id": occs[i-1] if i>0 else "",
                               "next_occurrence_id": occs[i+1] if i<n-1 else "",
                               "gap_from_prev": prev_gap, "gap_to_next": next_gap,
                               "adjacent_to_prev": adj_prev, "adjacent_to_next": adj_next,
                               "edge_scope": "all_occurrences",
                               "domain_order_evidence": str(r.get("domain_order_evidence","")).strip()
                                   or "ComputedFromSortedCoordinates"})
    edge_df = pd.DataFrame(edges).fillna("")
    edge_df = ensure_columns(edge_df, empty_domain_order_edge_columns())
    return df, edge_df

def load_domain_type_info(domain_dir: Path) -> "pd.DataFrame":
    files = sorted(domain_dir.glob("interpro_domains_*.tsv"))
    if not files: return pd.DataFrame()
    return pd.concat([pd.read_csv(f,sep="\t",dtype=str,encoding="utf-8").fillna("") for f in files],
                     ignore_index=True).fillna("")

def _ndi_pick_col(r: "pd.Series", *names: str) -> str:
    for n in names:
        if n in r and str(r.get(n,"")).strip(): return str(r.get(n,"")).strip()
    return ""

def build_domain_types_table(df_occ: "pd.DataFrame", df_ipr: "pd.DataFrame") -> "pd.DataFrame":
    rows: Dict[str, Dict[str, Any]] = {}
    if not df_ipr.empty:
        for _, r in df_ipr.iterrows():
            acc = _ndi_pick_col(r,"entry_accession","ipr","accession")
            db  = normalize_db_label(_ndi_pick_col(r,"source_database","source_db","db"))
            if not acc or not db: continue
            dt_id = f"{db}:{acc}"
            if dt_id in rows: continue
            rows[dt_id] = {"domain_type_id": dt_id, "source_db": db, "accession": acc,
                           "name": _ndi_pick_col(r,"name"),
                           "integrated_interpro_id": _ndi_pick_col(r,"integrated_ipr","integrated"),
                           "description": _ndi_pick_col(r,"description"),
                           "go_terms_all": _ndi_pick_col(r,"go_terms_all","go_terms"),
                           "go_terms_mf": _ndi_pick_col(r,"go_terms_mf"),
                           "go_terms_bp": _ndi_pick_col(r,"go_terms_bp"),
                           "go_terms_cc": _ndi_pick_col(r,"go_terms_cc"),
                           "function_tags": _ndi_pick_col(r,"function_tags"),
                           "evidence_level": _ndi_pick_col(r,"evidence_level")}
    if not df_occ.empty:
        for _, r in df_occ.iterrows():
            dt_id = str(r.get("domain_type_id","")).strip()
            if not dt_id or ":" not in dt_id: continue
            if dt_id not in rows:
                db, acc = dt_id.split(":",1)
                rows[dt_id] = {"domain_type_id": dt_id, "source_db": normalize_db_label(db),
                                "accession": acc, "name": str(r.get("label","")).strip(),
                                "integrated_interpro_id": "", "description": "",
                                "go_terms_all": str(r.get("go_terms_all","")).strip(),
                                "go_terms_mf": str(r.get("go_terms_mf","")).strip(),
                                "go_terms_bp": str(r.get("go_terms_bp","")).strip(),
                                "go_terms_cc": str(r.get("go_terms_cc","")).strip(),
                                "function_tags": str(r.get("function_tags","")).strip(),
                                "evidence_level": str(r.get("domain_evidence_level","")).strip()}
            else:
                if not rows[dt_id].get("name"):
                    rows[dt_id]["name"] = str(r.get("label","")).strip()
                if not rows[dt_id].get("function_tags"):
                    rows[dt_id]["function_tags"] = str(r.get("function_tags","")).strip()
    df_types = pd.DataFrame(list(rows.values())).fillna("")
    if not df_types.empty:
        df_types = df_types.sort_values(["source_db","accession","domain_type_id"])
    return df_types

def build_domain_type_cooccurrence(df_occ: "pd.DataFrame", df_types: "pd.DataFrame") -> "pd.DataFrame":
    if df_occ.empty: return pd.DataFrame(columns=empty_domain_type_cooccurrence_columns())
    name_map = ({str(r["domain_type_id"]): str(r.get("name","")).strip() or str(r["domain_type_id"])
                 for _, r in df_types.iterrows()} if not df_types.empty else {})
    grouped = df_occ.groupby("protein_accession")["domain_type_id"].apply(
        lambda s: sorted({d for d in s.tolist() if str(d).strip()}))
    co_counts: Dict = defaultdict(int); co_proteins: Dict = defaultdict(set)
    for prot, dlist in grouped.items():
        if len(dlist) < 2: continue
        for a, b in combinations(dlist, 2):
            key = tuple(sorted((a,b))); co_counts[key] += 1; co_proteins[key].add(prot)
    rows = []
    for (a,b), _ in sorted(co_counts.items(), key=lambda x: (-x[1],x[0][0],x[0][1])):
        prots = sorted(co_proteins[(a,b)])
        rows.append({"domain_type_id_A": a, "domain_type_name_A": name_map.get(a,a),
                     "domain_type_id_B": b, "domain_type_name_B": name_map.get(b,b),
                     "supportCount_proteins": len(prots), "proteins_with_pair": ", ".join(prots),
                     "association_type": "CoOccurrencePatternContext",
                     "association_evidence": "ComputedFromProteinDomainSet"})
    return ensure_columns(pd.DataFrame(rows).fillna(""), empty_domain_type_cooccurrence_columns())

def explode_occurrence_links(df_occ: "pd.DataFrame") -> Tuple["pd.DataFrame","pd.DataFrame","pd.DataFrame"]:
    if df_occ.empty:
        return (pd.DataFrame(columns=empty_occurrence_structural_link_columns()),
                pd.DataFrame(columns=empty_occurrence_binding_link_columns()),
                pd.DataFrame(columns=empty_occurrence_signature_link_columns()))
    struct_rows: List[Dict] = []; bind_rows: List[Dict] = []; sig_rows: List[Dict] = []
    if "cath_superfamilies" in df_occ.columns:
        for _, r in df_occ.iterrows():
            occ = str(r.get("occurrence_id","")).strip()
            if not occ: continue
            for c in _ndi_split_ids(str(r.get("cath_superfamilies","")).replace("|",";")):
                if c: struct_rows.append({"occurrence_id": occ, "structural_classification_id": c,
                                           "structural_db": "CATH-Gene3D", "link_type": "hasStructuralClassification",
                                           "evidence": str(r.get("structural_context_evidence","")).strip() or "InterProStructuralHitsOverlap"})
    if "binding_target_chebi" in df_occ.columns:
        for _, r in df_occ.iterrows():
            occ = str(r.get("occurrence_id","")).strip()
            if not occ: continue
            for cid in _ndi_split_ids(str(r.get("binding_target_chebi","")).replace("|",";")):
                if cid: bind_rows.append({"occurrence_id": occ, "binding_target_id": cid,
                                           "binding_db": "ChEBI", "link_type": "hasBindingTarget",
                                           "evidence": str(r.get("binding_target_class_evidence","")).strip() or "UniProtFT:ligand_id"})
    if "signatures" in df_occ.columns:
        for _, r in df_occ.iterrows():
            occ = str(r.get("occurrence_id","")).strip()
            if not occ: continue
            for s in _ndi_split_ids(str(r.get("signatures","")).replace("|",";")):
                if s:
                    db = s.split(":",1)[0] if ":" in s else ""
                    sig_rows.append({"occurrence_id": occ, "signature_id": s,
                                     "signature_db": normalize_db_label(db),
                                     "link_type": "hasSignature", "evidence": "InterProMemberOverlapAssignment"})
    df_struct = ensure_columns(pd.DataFrame(struct_rows).fillna(""), empty_occurrence_structural_link_columns())
    df_bind   = ensure_columns(pd.DataFrame(bind_rows).fillna(""),   empty_occurrence_binding_link_columns())
    df_sig    = ensure_columns(pd.DataFrame(sig_rows).fillna(""),    empty_occurrence_signature_link_columns())
    return df_struct, df_bind, df_sig


def _export_main(argv: Optional[List[str]] = None) -> None:
    """Ontology export builder — builds TSV tables from per-protein pipeline outputs."""
    import argparse as _ap
    ap = _ap.ArgumentParser(description="Build ontology-aligned TSV exports (no API calls).")
    ap.add_argument("--proteins",   required=True)
    ap.add_argument("--domain-dir", required=False)
    ap.add_argument("--outdir",     required=True)
    ap.add_argument("--pname",      required=False)
    args = ap.parse_args(argv)
    proteins_path = Path(args.proteins)
    domain_dir    = Path(args.domain_dir) if args.domain_dir else proteins_path.parent
    out_dir       = Path(args.outdir);  out_dir.mkdir(parents=True, exist_ok=True)
    pname = (args.pname or out_dir.name).strip()
    df_sum  = load_proteins_summary(proteins_path)
    df_prot = build_protein_table(df_sum)
    df_prot.to_csv(out_dir/f"{pname}_proteins.tsv", sep="\t", index=False, encoding="utf-8")
    print(f"[INFO] Wrote proteins: {len(df_prot)} rows")
    df_filt = filter_disordered_regions(df_sum)
    keep_cols = [c for c in ["accession","entry_id","genes","protein_name","protein_length",
                              "go_ids_all","go_ids_mf","go_ids_bp","go_ids_cc",
                              "subcellular_location_text","reactome_ids","kegg_ids",
                              "regions_pretty","regions_filtered_pretty","region_remained_count",
                              "domains_pretty","motifs_pretty","interpro_domain_instance_count",
                              "domain_instances_pretty"] if c in df_filt.columns]
    df_filt[keep_cols].to_csv(out_dir/f"{pname}_proteins_filtered.tsv", sep="\t", index=False, encoding="utf-8")
    df_occ_raw = load_domain_instances(domain_dir)
    if df_occ_raw.empty:
        print(f"[WARN] No domain_instances_*.tsv in {domain_dir}")
        df_occ   = pd.DataFrame(columns=empty_domain_occurrence_columns())
        df_order = pd.DataFrame(columns=empty_domain_order_edge_columns())
    else:
        # Deduplicate at the per-gene level too (defensive — normally each gene
        # has a single protein set, but guard against incremental re-runs that
        # accidentally left old domain_instances files from an earlier wider query).
        dup_key_local = [c for c in ["protein_accession","start","end","label_acc"]
                         if c in df_occ_raw.columns]
        if dup_key_local:
            n_before = len(df_occ_raw)
            df_occ_raw = (df_occ_raw
                          .drop_duplicates(subset=dup_key_local, keep="first")
                          .reset_index(drop=True))
            if len(df_occ_raw) < n_before:
                print(f"[DEDUP] Per-gene: removed {n_before - len(df_occ_raw)} "
                      f"duplicate instance(s) in {domain_dir.name}")
        df_occ, df_order = assign_occurrence_ids_preserve_plan_b(df_occ_raw)
    df_ipr   = load_domain_type_info(domain_dir)
    df_types = build_domain_types_table(df_occ, df_ipr)
    df_types.to_csv(out_dir/f"{pname}_domain_types.tsv", sep="\t", index=False, encoding="utf-8")
    print(f"[INFO] Wrote domain types: {len(df_types)} rows")
    if "binding_parnter" not in df_occ.columns: df_occ["binding_parnter"] = ""
    df_occ_export = ensure_columns(df_occ, empty_domain_occurrence_columns())[empty_domain_occurrence_columns()]
    df_occ_export.to_csv(out_dir/f"{pname}_domain_occurrences.tsv", sep="\t", index=False, encoding="utf-8")
    print(f"[INFO] Wrote domain occurrences: {len(df_occ_export)} rows")
    df_order_export = ensure_columns(df_order, empty_domain_order_edge_columns())[empty_domain_order_edge_columns()]
    df_order_export.to_csv(out_dir/f"{pname}_domain_order_edges.tsv", sep="\t", index=False, encoding="utf-8")
    co_df = build_domain_type_cooccurrence(df_occ, df_types)
    ensure_columns(co_df, empty_domain_type_cooccurrence_columns())[empty_domain_type_cooccurrence_columns()]\
        .to_csv(out_dir/f"{pname}_domain_type_cooccurrence.tsv", sep="\t", index=False, encoding="utf-8")
    df_struct, df_bind, df_sig = explode_occurrence_links(df_occ)
    ensure_columns(df_struct, empty_occurrence_structural_link_columns())[empty_occurrence_structural_link_columns()]\
        .to_csv(out_dir/f"{pname}_occurrence_structural_links.tsv", sep="\t", index=False, encoding="utf-8")
    ensure_columns(df_bind, empty_occurrence_binding_link_columns())[empty_occurrence_binding_link_columns()]\
        .to_csv(out_dir/f"{pname}_occurrence_binding_links.tsv", sep="\t", index=False, encoding="utf-8")
    ensure_columns(df_sig, empty_occurrence_signature_link_columns())[empty_occurrence_signature_link_columns()]\
        .to_csv(out_dir/f"{pname}_occurrence_signature_links.tsv", sep="\t", index=False, encoding="utf-8")
    print(f"[INFO] Export complete -> {out_dir}")

def normalize_positional_category(x: str) -> str:
    """
    Normalise a positional category string to one of the three canonical
    vocabulary terms: Nterminal / Internal / Cterminal (or unknown).

    Accepts common variants and the legacy 'Central' label, which is
    mapped to 'Internal' (Central implied geometric midpoint; Internal
    correctly means between the two termini regardless of exact position).
    """
    s = (x or "").strip().lower()
    if not s:
        return "unknown"
    # Exact matches first
    if s in ("nterminal", "n-terminal", "n_terminal", "nterm", "n"):
        return "Nterminal"
    if s in ("cterminal", "c-terminal", "c_terminal", "cterm", "c"):
        return "Cterminal"
    if s in ("internal", "central", "middle", "intr", "i"):
        return "Internal"
    # Prefix fallback — "cent..." must map to Internal, not Cterminal
    if s.startswith("n"):
        return "Nterminal"
    if s.startswith("c") and not s.startswith("cent"):
        return "Cterminal"
    if s.startswith("i") or s.startswith("cent"):
        return "Internal"
    return "unknown"


# ═══════════════════════════════════════════════════════════════════════════════
# PART 5 — Per-gene orchestration
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class RunResult:
    gene: str
    query: str
    base_dir: Path
    proteins_summary: Path
    export_dir: Path
    exports: Dict[str, Path]


def _assert_exists(p: Path, what: str) -> None:
    if not p.exists():
        raise SystemExit(f"[FATAL] Missing {what}: {p}")

def _first_hit(root: Path, patterns: List[str]) -> Optional[Path]:
    for pat in patterns:
        hits = sorted(root.glob(pat))
        if hits: return hits[0]
    return None

def _find_export_files(export_dir: Path, pname: str) -> Dict[str, Path]:
    expected = {
        "proteins":                    f"{pname}_proteins.tsv",
        "proteins_filtered":           f"{pname}_proteins_filtered.tsv",
        "domain_types":                f"{pname}_domain_types.tsv",
        "domain_occurrences":          f"{pname}_domain_occurrences.tsv",
        "domain_order_edges":          f"{pname}_domain_order_edges.tsv",
        "domain_type_cooccurrence":    f"{pname}_domain_type_cooccurrence.tsv",
        "occurrence_structural_links": f"{pname}_occurrence_structural_links.tsv",
        "occurrence_binding_links":    f"{pname}_occurrence_binding_links.tsv",
        "occurrence_signature_links":  f"{pname}_occurrence_signature_links.tsv",
    }
    out: Dict[str, Path] = {}
    for k, fname in expected.items():
        p = export_dir / fname
        if p.exists(): out[k] = p
    if "domain_occurrences" not in out:
        c = _first_hit(export_dir, [f"*{pname}*domain_occurrences*.tsv", "*domain_occurrences*.tsv"])
        if c: out["domain_occurrences"] = c
    if "domain_types" not in out:
        c = _first_hit(export_dir, [f"*{pname}*domain_types*.tsv", "*domain_types*.tsv"])
        if c: out["domain_types"] = c
    return out

def _safe_read_tsv(path: Path) -> Optional["pd.DataFrame"]:
    """Read a TSV file safely; return None on empty/unreadable files."""
    if not _PANDAS_OK: return None
    try:
        if not path.exists() or path.stat().st_size == 0:
            return None
        df = pd.read_csv(path, sep="\t", dtype=str, encoding="utf-8")
        if df is None or len(df.columns) == 0:
            return None
        return df.fillna("")
    except EmptyDataError:
        print(f"[WARN] Skipping blank TSV: {path}")
        return None
    except Exception as e:
        print(f"[WARN] Skipping unreadable TSV: {path} | {e}")
        return None

def run_pipeline_for_gene(
        gene: str,
        query: str,
        reviewed: str,
        base_dir: Path,
        fetch_interpro: bool,
        interpro_mode: str,
        save_interpro_json: bool,
        force: bool,
        pname: Optional[str] = None,
) -> RunResult:
    """
    Run Step A (UniProt + InterPro acquisition) then Step B (ontology export)
    for a single gene query. Skips either step if outputs already exist and
    force is False.
    """
    pname = pname or gene
    base_dir.mkdir(parents=True, exist_ok=True)

    proteins_summary = base_dir / "proteins_summary.tsv"
    export_dir = base_dir / "out"
    export_dir.mkdir(parents=True, exist_ok=True)

    need_uniprot = force or (not proteins_summary.exists())
    existing_occ = _first_hit(export_dir, [f"{pname}_domain_occurrences.tsv",
                                            "*domain_occurrences*.tsv"])
    # Re-run Step B if the file is missing OR if it exists but is empty/header-only.
    # Empty files are left behind by previously-failed Step B runs (e.g. NameError).
    # A header-only TSV has <=200 bytes; real data is always larger.
    _occ_is_empty = (existing_occ is not None and
                     existing_occ.exists() and
                     existing_occ.stat().st_size <= 200)
    need_export  = force or (existing_occ is None) or _occ_is_empty
    if _occ_is_empty:
        print(f"[B] {gene}: existing domain_occurrences.tsv is empty — will regenerate")

    # Step A — data acquisition
    if need_uniprot:
        print(f"[A] {gene}: fetching UniProt + InterPro -> {base_dir}")
        ret = run_uniprot_summary(
            query=query,
            reviewed=reviewed,
            outdir=str(base_dir),
            fetch_interpro=fetch_interpro,
            interpro_mode=interpro_mode,
            save_interpro_json=save_interpro_json,
        )
        if ret:
            base_dir = Path(ret)
            proteins_summary = base_dir / "proteins_summary.tsv"
            export_dir = base_dir / "out"
            export_dir.mkdir(parents=True, exist_ok=True)
    else:
        print(f"[A] {gene}: SKIP (found {proteins_summary.name})")

    _assert_exists(proteins_summary, "proteins_summary.tsv")

    # Step B — ontology export
    if need_export:
        print(f"[B] {gene}: building ontology tables -> {export_dir}")
        argv = [
            "--proteins",   str(proteins_summary),
            "--domain-dir", str(base_dir),
            "--outdir",     str(export_dir),
            "--pname",      pname,
        ]
        try:
            _export_main(argv)
        except TypeError:
            saved = sys.argv
            try:
                sys.argv = ["01_data_and_merge.py"] + argv
                _export_main()
            finally:
                sys.argv = saved
    else:
        print(f"[B] {gene}: SKIP export (outputs already in {export_dir})")

    exports = _find_export_files(export_dir, pname=pname)
    if "domain_occurrences" not in exports:
        raise SystemExit(
            f"[FATAL] {gene}: domain_occurrences table not found in {export_dir}"
        )

    print(f"[OK] {gene}: {exports['domain_occurrences'].name}")
    return RunResult(
        gene=gene, query=query, base_dir=base_dir,
        proteins_summary=proteins_summary,
        export_dir=export_dir, exports=exports,
    )


# ═══════════════════════════════════════════════════════════════════════════════
# PART 6 — Pathway-level aggregation
#           Integrated from merged_features_summary.py
#           Replaces the simple merge_batch_tables() with a comprehensive
#           version that also produces the domain_type_summary.tsv required
#           by downstream Stage 4B (transfer.py).
# ═══════════════════════════════════════════════════════════════════════════════

# ── Aggregation helpers ────────────────────────────────────────────────────────

def _ensure_cols(df: "pd.DataFrame", cols: List[str]) -> "pd.DataFrame":
    """Add any missing columns as empty strings without modifying existing ones."""
    df = df.copy()
    for c in cols:
        if c not in df.columns:
            df[c] = ""
    return df

def _infer_gene_from_filename(fn: str) -> str:
    """
    Extract the gene name from an output filename.
    e.g. 'ADCY1_domain_occurrences.tsv' -> 'ADCY1'
    """
    m = re.match(
        r"^(.+?)_(domain_occurrences|domain_types|domain_order_edges|occurrence_.*_links)\.tsv$",
        fn, flags=re.I,
    )
    return m.group(1) if m else Path(fn).stem

def _split_multi_value(cell: str) -> List[str]:
    """Split a multi-value cell (separated by ; , | or whitespace) into tokens."""
    if cell is None: return []
    s = str(cell).strip()
    if not s: return []
    parts = []
    for chunk in _SPLIT_MULTI.split(s):
        chunk = chunk.strip()
        if not chunk: continue
        for tok in _SPLIT_WS.split(chunk):
            tok = tok.strip()
            if tok: parts.append(tok)
    return parts

def _union_cells(series: "pd.Series", max_items: int = 200) -> str:
    """Collect all unique tokens across a column and return a sorted '; '-joined string."""
    vals: set = set()
    for cell in series.fillna(""):
        for tok in _split_multi_value(cell):
            vals.add(tok)
    out = sorted(vals)
    if max_items and len(out) > max_items:
        out = out[:max_items]
    return "; ".join(out)

def _count_multicats(series: "pd.Series") -> str:
    """
    Tally token frequency across a column.
    Returns e.g. 'Cterminal:22; Nterminal:4; unknown:1'
    """
    counts: Dict[str, int] = {}
    for cell in series.fillna(""):
        toks = _split_multi_value(cell)
        if not toks:
            toks = ["unknown"]
        for t in toks:
            counts[t] = counts.get(t, 0) + 1
    items = sorted(counts.items(), key=lambda x: (-x[1], x[0]))
    return "; ".join([f"{k}:{v}" for k, v in items])

def _series_to_int(series: "pd.Series") -> "pd.Series":
    """Coerce a Series to int, filling non-numeric values with 0."""
    return pd.to_numeric(series, errors="coerce").fillna(0).astype(int)

def _domain_length_stats(sub: "pd.DataFrame") -> Tuple[str, str, str, str]:
    """Return (mean, median, min, max) domain length as strings."""
    s = pd.to_numeric(sub["start"], errors="coerce")
    e = pd.to_numeric(sub["end"],   errors="coerce")
    lens = (e - s + 1).dropna()
    if lens.empty:
        return "", "", "", ""
    return (f"{lens.mean():.2f}", f"{lens.median():.2f}",
            str(int(lens.min())), str(int(lens.max())))

def _copy_number_profile(sub: "pd.DataFrame") -> str:
    """
    Per-protein max copy-number distribution for this domain type.
    Returns e.g. '1:18; 2:3' meaning 18 proteins have 1 copy, 3 have 2 copies.
    """
    tmp = sub.copy()
    tmp["_cn"] = _series_to_int(tmp.get("copy_number", pd.Series(dtype=str)))
    per_prot = tmp.groupby("protein_accession")["_cn"].max()
    counts = per_prot.value_counts().to_dict()
    items = sorted(counts.items(), key=lambda x: x[0])
    return "; ".join([f"{k}:{v}" for k, v in items])

def _sample_accessions(accs: List[str], n: int) -> str:
    accs = sorted({a for a in accs if a})
    if len(accs) <= n:
        return ", ".join(accs)
    return ", ".join(accs[:n]) + f", ... (+{len(accs) - n} more)"

def _pick_col(df: "pd.DataFrame", candidates: List[str]) -> Optional[str]:
    """Return the first candidate column name that exists in df."""
    return next((c for c in candidates if c in df.columns), None)


# ── Core merge functions ───────────────────────────────────────────────────────

def _glob_files(root: Path, pattern: str) -> List[Path]:
    return sorted(root.glob(pattern))

def _dedup_occurrences(df: "pd.DataFrame", context: str = "") -> "pd.DataFrame":
    """
    Remove duplicate domain occurrence rows that arise when the same protein
    accession is returned by more than one gene-query in Stage 1.

    Duplication mechanism
    ─────────────────────
    Stage 1 queries UniProt with patterns like 'gene:ADCY* AND organism_id:9606'.
    A protein such as ADCYAP1R1 (P41586) is returned by three separate queries:
      • gene:ADCY*       (prefix match on the gene name)
      • gene:ADCYAP1*
      • gene:ADCYAP1R1*
    Each query writes its own *_out_final/ directory.  When Stage 2 concatenates
    all per-gene domain_occurrences.tsv files, P41586's domain rows appear three
    times — identical in every column except 'gene' and 'source_file'.

    Strategy
    ────────
    Duplicate rows are identified by DEDUP_KEY (protein_accession + start + end +
    domain_type_id).  Among duplicates the "richest" row is kept — the one with
    the highest count of non-empty cells — so the maximum available metadata is
    preserved.  In practice all duplicates are fully identical, so 'first' would
    work too; the richness sort is a safety net for edge cases.

    Parameters
    ──────────
    df      : concatenated occurrence DataFrame (may contain duplicates)
    context : label printed in log messages (e.g. pathway name)

    Returns
    ───────
    Deduplicated DataFrame with the same column set as the input.
    """
    if not _PANDAS_OK:
        return df

    # Build dedup key from whichever DEDUP_KEY columns exist in df
    key = [c for c in DEDUP_KEY if c in df.columns]
    if not key:
        print(f"[DEDUP] {context}No DEDUP_KEY columns found — skipping deduplication.")
        return df

    before = len(df)
    dup_mask = df.duplicated(subset=key, keep=False)
    n_dups = int(dup_mask.sum())

    if n_dups == 0:
        print(f"[DEDUP] {context}{before} rows — no duplicates found.")
        return df

    # ── Report which proteins and genes caused duplicates ─────────────────────
    dup_df = df[dup_mask]
    affected_proteins = sorted(dup_df["protein_accession"].unique()) if "protein_accession" in dup_df.columns else []
    affected_genes    = sorted(dup_df["gene"].unique())              if "gene"              in dup_df.columns else []
    print(f"[DEDUP] {context}{before} rows — found {n_dups} duplicate rows "
          f"across {len(affected_proteins)} protein(s).")
    if affected_genes:
        print(f"[DEDUP]   Caused by overlapping gene queries: {affected_genes}")
    if affected_proteins:
        # Per-protein duplicate counts
        per_prot = (dup_df.groupby("protein_accession").size()
                    if "protein_accession" in dup_df.columns else {})
        for prot, cnt in sorted(per_prot.items(), key=lambda x: -x[1]):
            genes_for_prot = (sorted(dup_df[dup_df["protein_accession"] == prot]["gene"].unique())
                              if "gene" in dup_df.columns else [])
            print(f"[DEDUP]     {prot}: {cnt} duplicate rows  "
                  f"(fetched by genes: {genes_for_prot})")

    # ── Keep the richest row per key (max non-empty cells) ────────────────────
    df = df.copy()
    df["_richness"] = (df != "").sum(axis=1)
    df = (df.sort_values("_richness", ascending=False)
            .drop_duplicates(subset=key, keep="first")
            .drop(columns=["_richness"])
            .reset_index(drop=True))

    after = len(df)
    print(f"[DEDUP]   Removed {before - after} duplicate row(s) → {after} rows remaining.")
    return df


def _merge_occurrences_from_dir(root: Path) -> "pd.DataFrame":
    """
    Glob all per-gene domain_occurrences.tsv files under root, concatenate
    them, and enforce OCC_MIN_COLS.

    Fallback: if no domain_occurrences.tsv files exist or all are empty,
    attempt to rebuild from domain_instances_*.tsv files directly using
    _export_main.  This handles the case where Step B previously failed
    (e.g. due to a NameError) and left behind empty domain_occurrences.tsv
    files that blocked Step B from re-running.
    """
    files = _glob_files(root, PAT_OCC)

    dfs = []
    if files:
        for p in files:
            df = _safe_read_tsv(p)
            if df is None or df.empty:
                continue
            df = _ensure_cols(df, OCC_MIN_COLS)
            gene = _infer_gene_from_filename(p.name)
            df.insert(0, "gene",        gene)
            df.insert(1, "source_file", str(p))
            dfs.append(df)

    if not dfs:
        # ── Fallback: rebuild from domain_instances_*.tsv files ──────────────
        # This happens when Step B previously failed and left empty files, or
        # was never run.  We locate each *_out_final/ directory, find its
        # domain_instances_*.tsv files, and run _export_main to produce the
        # domain_occurrences.tsv before merging.
        print("[WARN] No non-empty domain_occurrences.tsv files found.")
        print("[WARN] Attempting fallback rebuild from domain_instances_*.tsv files...")
        gene_dirs = sorted([d for d in root.iterdir()
                            if d.is_dir() and d.name.endswith("_out_final")])
        if not gene_dirs:
            raise SystemExit(
                f"[FATAL] No domain_occurrences or domain_instances files found under {root}.\n"
                "  Run Stage 1 first with --fetch-interpro."
            )
        rebuilt = 0
        for gdir in gene_dirs:
            psum = gdir / "proteins_summary.tsv"
            if not psum.exists():
                print(f"  [SKIP] {gdir.name}: no proteins_summary.tsv")
                continue
            gene = gdir.name.replace("_out_final", "")
            export_d = gdir / "out"
            export_d.mkdir(parents=True, exist_ok=True)
            print(f"  [REBUILD] {gene}: running Step B -> {export_d}")
            argv = ["--proteins", str(psum),
                    "--domain-dir", str(gdir),
                    "--outdir",     str(export_d),
                    "--pname",      gene]
            try:
                _export_main(argv)
                rebuilt += 1
            except Exception as _e:
                print(f"  [ERROR] {gene}: Step B failed: {_e}")
                continue
        if rebuilt == 0:
            raise SystemExit(
                "[FATAL] Fallback rebuild failed for all genes.\n"
                "  Check that domain_instances_*.tsv files exist in the *_out_final/ directories."
            )
        # Re-glob after rebuild
        files = _glob_files(root, PAT_OCC)
        for p in files:
            df = _safe_read_tsv(p)
            if df is None or df.empty:
                continue
            df = _ensure_cols(df, OCC_MIN_COLS)
            gene = _infer_gene_from_filename(p.name)
            df.insert(0, "gene",        gene)
            df.insert(1, "source_file", str(p))
            dfs.append(df)

    if not dfs:
        raise SystemExit("[FATAL] All domain_occurrences files were empty or unreadable "
                         "even after fallback rebuild.")

    merged = pd.concat(dfs, ignore_index=True).fillna("")
    print(f"[MERGE] Occurrences: {len(files)} gene file(s) -> {len(merged)} rows")
    merged = _dedup_occurrences(merged, context="")
    return merged

def _merge_types_from_dir(root: Path) -> "pd.DataFrame":
    """
    Glob all per-gene domain_types.tsv files under root, concatenate them,
    and deduplicate by domain_type_id (keeping the row with the most metadata).
    """
    files = _glob_files(root, PAT_TYPES)
    if not files:
        return pd.DataFrame(columns=["gene", "source_file"] + TYPE_MIN_COLS)

    dfs = []
    for p in files:
        df = _safe_read_tsv(p)
        if df is None or df.empty:
            continue
        df = _ensure_cols(df, TYPE_MIN_COLS)
        gene = _infer_gene_from_filename(p.name)
        df.insert(0, "gene",        gene)
        df.insert(1, "source_file", str(p))
        dfs.append(df)

    if not dfs:
        return pd.DataFrame(columns=["gene", "source_file"] + TYPE_MIN_COLS)

    merged = pd.concat(dfs, ignore_index=True).fillna("")
    # Sort so the richest row (most metadata) survives dedup
    merged = merged.sort_values(
        by=["domain_type_id", "name", "description"],
        ascending=[True, False, False],
        kind="mergesort",
    )
    merged = merged.drop_duplicates(subset=["domain_type_id"], keep="first")
    return merged

def _merge_link_table(root: Path, pattern: str, out_path: Path) -> int:
    """Merge one category of link tables (order, binding, structural, signature)."""
    files = _glob_files(root, pattern)
    dfs = []
    for p in files:
        df = _safe_read_tsv(p)
        if df is None or df.empty: continue
        gene = _infer_gene_from_filename(p.name)
        df.insert(0, "gene",        gene)
        df.insert(1, "source_file", str(p))
        dfs.append(df)
    if not dfs:
        return 0
    pd.concat(dfs, ignore_index=True).fillna("").to_csv(
        out_path, sep="\t", index=False, encoding="utf-8"
    )
    return sum(len(df) for df in dfs)


# ── Domain type summary ────────────────────────────────────────────────────────

def _summarize_domain_types(
        occ: "pd.DataFrame",
        types: Optional["pd.DataFrame"],
        full_accessions: bool = False,
        accession_sample_n: int = 50,
        union_max_items: int = 200,
) -> "pd.DataFrame":
    """
    Build the domain_type_summary.tsv: one row per unique domain type with
    aggregated statistics across all proteins in the pathway.

    This is the primary input for Stage 4B (transfer.py / type annotation).
    """
    if occ.empty:
        return pd.DataFrame()

    occ = _ensure_cols(occ, OCC_MIN_COLS + ["gene"])

    react_col = _pick_col(occ, ["protein_reactome_ids", "reactome_ids"])
    kegg_col  = _pick_col(occ, ["protein_kegg_ids", "kegg_ids"])

    # Build name/description/evidence lookup from the types table
    type_lookup: Dict[str, Dict[str, str]] = {}
    if types is not None and not types.empty and "domain_type_id" in types.columns:
        for _, r in types.iterrows():
            dt = str(r.get("domain_type_id", "")).strip()
            if dt:
                type_lookup[dt] = {
                    "domain_name":       str(r.get("name", "")).strip(),
                    "type_description":  str(r.get("description", "")).strip(),
                    "type_evidence_level": str(r.get("evidence_level", "")).strip(),
                }

    rows = []
    for dt_id, sub in occ.groupby("domain_type_id", dropna=False):
        dt_id = str(dt_id).strip()
        if not dt_id:
            continue

        proteins = sorted({str(x).strip() for x in sub["protein_accession"] if str(x).strip()})
        genes    = sorted({str(x).strip() for x in sub["gene"] if str(x).strip()})

        # Best available name
        name_guess = type_lookup.get(dt_id, {}).get("domain_name", "")
        if not name_guess:
            vc = sub["label"].astype(str).str.strip()
            vc = vc[vc != ""]
            if not vc.empty:
                name_guess = vc.value_counts().index[0]

        mean_len, med_len, min_len, max_len = _domain_length_stats(sub)

        reactome_union = _union_cells(sub[react_col], max_items=union_max_items) if react_col else ""
        kegg_union     = _union_cells(sub[kegg_col],  max_items=union_max_items) if kegg_col  else ""

        rows.append({
            "domain_type_id":   dt_id,
            "domain_name":      name_guess,

            "n_occurrences":    int(len(sub)),
            "n_proteins":       int(len(proteins)),
            "n_genes":          int(len(genes)),
            "genes":            ", ".join(genes),

            "copy_number_profile_per_protein": _copy_number_profile(sub),

            "mean_domain_length":   mean_len,
            "median_domain_length": med_len,
            "min_domain_length":    min_len,
            "max_domain_length":    max_len,

            "positional_category_counts": _count_multicats(sub["positional_category"]),
            "topology_context_counts":    _count_multicats(sub["topology_context"]),
            "proximity_category_counts":  _count_multicats(sub["proximity_categories"]),

            "reactome_ids_union": reactome_union,
            "kegg_ids_union":     kegg_union,

            "go_terms_all_union": _union_cells(sub["go_terms_all"], max_items=union_max_items),
            "go_terms_mf_union":  _union_cells(sub["go_terms_mf"],  max_items=union_max_items),
            "go_terms_bp_union":  _union_cells(sub["go_terms_bp"],  max_items=union_max_items),
            "go_terms_cc_union":  _union_cells(sub["go_terms_cc"],  max_items=union_max_items),

            "function_tags_union":        _union_cells(sub["function_tags"],        max_items=union_max_items),
            "binding_target_chebi_union": _union_cells(sub["binding_target_chebi"], max_items=union_max_items),
            "signatures_union":           _union_cells(sub["signatures"],           max_items=union_max_items),

            # Manual column (filled during Stage 3 curation):
            # protein / nucleic_acid / lipid / metal / small_molecule / etc.
            "bind_target": "",

            "cath_superfamilies_union": _union_cells(sub["cath_superfamilies"], max_items=union_max_items),

            "protein_accessions": (
                ", ".join(proteins)
                if full_accessions
                else _sample_accessions(proteins, accession_sample_n)
            ),

            "type_description":    type_lookup.get(dt_id, {}).get("type_description", ""),
            "type_evidence_level": type_lookup.get(dt_id, {}).get("type_evidence_level", ""),
        })

    out = pd.DataFrame(rows).fillna("")
    out = out.sort_values(
        by=["n_proteins", "n_occurrences", "domain_type_id"],
        ascending=[False, False, True],
    )
    return out


# ═══════════════════════════════════════════════════════════════════════════════
# PART 6.5 — Domain sequence enrichment
#
# Adds five sequence-derived columns to the merged domain occurrence table.
# Sequences are fetched from the UniProt FASTA endpoint once per unique
# protein accession and cached in memory — no extra files are written.
#
# New columns:
#   start_position  int  Absolute start residue (1-based). Mirrors 'start'.
#                        Named explicitly for unambiguous OWL property mapping.
#   end_position    int  Absolute end residue (1-based).   Mirrors 'end'.
#   start_residue   str  Single-letter amino acid at start_position.
#   end_residue     str  Single-letter amino acid at end_position.
#   domain_sequence str  Full aa sequence of the domain (seq[start-1 : end]).
# ═══════════════════════════════════════════════════════════════════════════════

UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"


def _fetch_fasta_sequence(acc: str, delay: float = DEFAULT_DELAY) -> str:
    """
    Fetch the protein FASTA sequence for a UniProt accession.
    Returns the raw amino-acid sequence string (no header, no newlines).
    Returns an empty string on any error.
    """
    try:
        session = requests.Session()
        session.headers.update({"Accept": "text/plain", "User-Agent": UA})
        url = UNIPROT_FASTA_URL.format(acc=acc)
        txt = _http_get(session, url, expect="text", delay=delay)
        lines = txt.splitlines()
        return "".join(ln for ln in lines if ln and not ln.startswith(">")).strip()
    except Exception as e:
        print(f"[WARN] Could not fetch FASTA for {acc}: {e}")
        return ""


def _build_sequence_cache(
        accessions: List[str],
        delay: float = DEFAULT_DELAY,
) -> Dict[str, str]:
    """
    Fetch FASTA sequences for a list of unique UniProt accessions.
    Returns dict: accession -> sequence string (empty string on failure).
    """
    cache: Dict[str, str] = {}
    total = len(accessions)
    for i, acc in enumerate(accessions, 1):
        seq = _fetch_fasta_sequence(acc, delay=delay)
        cache[acc] = seq
        status = f"{len(seq)} aa" if seq else "FAILED"
        print(f"[SEQ] ({i}/{total}) {acc}: {status}")
    return cache


def enrich_domain_sequences(
        occ_df: "pd.DataFrame",
        seq_cache: Optional[Dict[str, str]] = None,
        delay: float = DEFAULT_DELAY,
) -> "pd.DataFrame":
    """
    Add five sequence-derived columns to a domain occurrence DataFrame.

    The five columns are inserted immediately after the 'end' column:
      start_position  — absolute start residue (1-based int, mirrors 'start')
      end_position    — absolute end residue   (1-based int, mirrors 'end')
      start_residue   — single amino acid at start_position
      end_residue     — single amino acid at end_position
      domain_sequence — full domain amino acid sequence (seq[start-1 : end])

    Parameters
    ----------
    occ_df : pd.DataFrame
        Merged domain occurrence table. Must have 'protein_accession',
        'start', and 'end' columns.
    seq_cache : dict, optional
        Pre-built {accession: sequence} mapping. Fetched automatically if None.
    delay : float
        Inter-request delay (seconds) used when fetching sequences.
    """
    if not _PANDAS_OK:
        print("[WARN] pandas unavailable — skipping sequence enrichment.")
        return occ_df

    df = occ_df.copy()

    if seq_cache is None:
        unique_accs = [
            acc for acc in df["protein_accession"].dropna().unique()
            if str(acc).strip()
        ]
        if not unique_accs:
            print("[WARN] No protein accessions found — skipping sequence enrichment.")
            return df
        print(f"[SEQ] Fetching FASTA for {len(unique_accs)} unique protein(s)...")
        seq_cache = _build_sequence_cache(unique_accs, delay=delay)

    start_positions: List[int] = []
    end_positions:   List[int] = []
    start_residues:  List[str] = []
    end_residues:    List[str] = []
    domain_seqs:     List[str] = []

    for _, row in df.iterrows():
        acc = str(row.get("protein_accession", "")).strip()
        seq = seq_cache.get(acc, "")

        try:
            s = int(float(str(row.get("start") or 0)))
        except (ValueError, TypeError):
            s = 0
        try:
            e = int(float(str(row.get("end") or 0)))
        except (ValueError, TypeError):
            e = 0

        start_positions.append(s)
        end_positions.append(e)

        if seq and s >= 1 and e >= s and e <= len(seq):
            start_residues.append(seq[s - 1])
            end_residues.append(seq[e - 1])
            domain_seqs.append(seq[s - 1 : e])
        else:
            start_residues.append("")
            end_residues.append("")
            domain_seqs.append("")

    new_cols = {
        "start_position":  start_positions,
        "end_position":    end_positions,
        "start_residue":   start_residues,
        "end_residue":     end_residues,
        "domain_sequence": domain_seqs,
    }

    # Insert immediately after 'end' column for readability
    if "end" in df.columns:
        idx = df.columns.tolist().index("end") + 1
        for offset, (col_name, col_data) in enumerate(new_cols.items()):
            df.insert(idx + offset, col_name, col_data)
    else:
        for col_name, col_data in new_cols.items():
            df[col_name] = col_data

    filled = sum(1 for s in domain_seqs if s)
    total  = len(domain_seqs)
    print(
        f"[SEQ] Enrichment complete: {filled}/{total} occurrences annotated "
        f"({total - filled} blank — check FASTA fetch warnings above)."
    )
    return df


# ── Master aggregation entry point ─────────────────────────────────────────────

def comprehensive_pathway_merge(
        root_outdir: Path,
        pathway_name: str,
        merge_links: bool = False,
        full_accessions: bool = False,
        accession_sample_n: int = 50,
        union_max_items: int = 200,
        skip_summary: bool = False,
) -> None:
    """
    Stage 2 — Pathway-level aggregation.

    Scans root_outdir for all per-gene output directories produced by Stage 1,
    then writes pathway-level tables.

    Always written
    ──────────────
    {pathway_name}_MERGED_domain_occurrences.tsv
        One row per domain occurrence across all genes, with five
        sequence-derived columns appended (unless --no-sequences):
          start_position, end_position, start_residue, end_residue,
          domain_sequence.
        Primary input for Stage 4A (ID injection).

    {pathway_name}_MERGED_domain_types.tsv
        One row per unique InterPro domain type (deduplicated).

    {pathway_name}_domain_type_summary.tsv
        Per-type statistics. Primary input for Stage 4B (transfer.py).

    Written only with --merge-links
    ────────────────────────────────
    {pathway_name}_MERGED_domain_order_edges.tsv
    {pathway_name}_MERGED_occurrence_binding_links.tsv
    {pathway_name}_MERGED_occurrence_structural_links.tsv
    {pathway_name}_MERGED_occurrence_signature_links.tsv
    """
    if not _PANDAS_OK:
        print("[WARN] pandas not available; skipping pathway-level aggregation.")
        return

    # 1 — Merge occurrence rows
    # The five sequence columns (start_position, end_position, start_residue,
    # end_residue, domain_sequence) are now computed in Stage 1 for each protein
    # while the sequence and domain coordinates are both in memory.  No separate
    # Stage 2.5 enrichment step or extra FASTA API calls are needed.
    occ = _merge_occurrences_from_dir(root_outdir)
    # _dedup_occurrences is already called inside _merge_occurrences_from_dir;
    # this second call is a safety net for any code path that bypasses it.
    occ = _dedup_occurrences(occ, context=f"{pathway_name}: ")

    seq_cols = ["start_position","end_position","start_residue","end_residue","domain_sequence"]
    present  = [c for c in seq_cols if c in occ.columns
                and occ[c].astype(str).str.strip().ne("").any()]
    missing  = [c for c in seq_cols if c not in present]
    if present:
        print(f"[SEQ] Sequence columns already populated from Stage 1: {present}")
    if missing:
        print(f"[SEQ] Note — the following sequence columns are absent/empty "
              f"(re-run Stage 1 with --fetch-interpro to populate): {missing}")

    occ_path = root_outdir / f"{pathway_name}_MERGED_domain_occurrences.tsv"
    occ.to_csv(occ_path, sep="\t", index=False, encoding="utf-8")
    print(f"[MERGE] -> {occ_path.name}  ({len(occ)} rows)")

    # 2 — Merge domain type rows
    types = _merge_types_from_dir(root_outdir)
    types_path = root_outdir / f"{pathway_name}_MERGED_domain_types.tsv"
    types.to_csv(types_path, sep="\t", index=False, encoding="utf-8")
    print(f"[MERGE] -> {types_path.name}  ({len(types)} rows)")

    # 3 — Domain type summary (needed by Stage 4B)
    if not skip_summary:
        summary = _summarize_domain_types(
            occ=occ, types=types,
            full_accessions=full_accessions,
            accession_sample_n=accession_sample_n,
            union_max_items=union_max_items,
        )
        sum_path = root_outdir / f"{pathway_name}_domain_type_summary.tsv"
        summary.to_csv(sum_path, sep="\t", index=False, encoding="utf-8")
        print(f"[MERGE] -> {sum_path.name}  ({len(summary)} domain types)")

    # 4 — Optional link tables
    if merge_links:
        for pat, label in [
            (PAT_ORDER,  "domain_order_edges"),
            (PAT_BIND,   "occurrence_binding_links"),
            (PAT_STRUCT, "occurrence_structural_links"),
            (PAT_SIG,    "occurrence_signature_links"),
        ]:
            out_path = root_outdir / f"{pathway_name}_MERGED_{label}.tsv"
            n = _merge_link_table(root_outdir, pat, out_path)
            if n:
                print(f"[MERGE] -> {out_path.name}  ({n} rows)")
            else:
                print(f"[MERGE] No {label} files found — skipped.")


# ═══════════════════════════════════════════════════════════════════════════════
# PART 7 — Batch runner
# ═══════════════════════════════════════════════════════════════════════════════

def run_pathway_batch(
        pathway_name: str,
        genes: List[str],
        reviewed: str,
        organism_id: str,
        root_outdir: Path,
        fetch_interpro: bool,
        interpro_mode: str,
        save_interpro_json: bool,
        force: bool,
        do_merge: bool,
        merge_links: bool = False,
        full_accessions: bool = False,
        accession_sample_n: int = 50,
        union_max_items: int = 200,
        skip_summary: bool = False,
) -> None:
    root_outdir.mkdir(parents=True, exist_ok=True)

    ok_results: List[RunResult] = []
    failures:   List[Tuple[str, str]] = []

    for gene in genes:
        query    = f"gene:{gene}* AND organism_id:{organism_id}"
        base_dir = root_outdir / f"{gene}_out_final"
        try:
            res = run_pipeline_for_gene(
                gene=gene, query=query, reviewed=reviewed,
                base_dir=base_dir,
                fetch_interpro=fetch_interpro,
                interpro_mode=interpro_mode,
                save_interpro_json=save_interpro_json,
                force=force, pname=gene,
            )
            ok_results.append(res)
        except Exception as e:
            failures.append((gene, str(e)))
            print(f"[ERROR] {gene} failed: {e}")
            traceback.print_exc()

    # Write batch status log
    if _PANDAS_OK:
        rows = [
            {"gene": r.gene, "query": r.query, "base_dir": str(r.base_dir),
             "domain_occurrences": str(r.exports.get("domain_occurrences", ""))}
            for r in ok_results
        ]
        rows += [
            {"gene": g, "query": f"gene:{g}* AND organism_id:{organism_id}",
             "base_dir": "", "domain_occurrences": "", "error": err}
            for g, err in failures
        ]
        pd.DataFrame(rows).fillna("").to_csv(
            root_outdir / f"{pathway_name}_BATCH_STATUS.tsv", sep="\t", index=False
        )

    if failures:
        fail_path = root_outdir / f"{pathway_name}_FAILED_GENES.tsv"
        fail_path.write_text(
            "\n".join([f"{g}\t{err}" for g, err in failures]), encoding="utf-8"
        )
        print(f"[BATCH] {len(failures)} gene(s) failed — see {fail_path}")

    # Stage 2 — comprehensive pathway-level aggregation
    if do_merge and ok_results:
        print(f"\n[STAGE 2] Pathway-level aggregation -> {root_outdir}")
        comprehensive_pathway_merge(
            root_outdir=root_outdir,
            pathway_name=pathway_name,
            merge_links=merge_links,
            full_accessions=full_accessions,
            accession_sample_n=accession_sample_n,
            union_max_items=union_max_items,
            skip_summary=skip_summary,
        )

    print(
        f"\n[DONE] {len(ok_results)} gene(s) completed, "
        f"{len(failures)} failed.\n"
        f"Output root: {root_outdir}\n"
        "\n*** NEXT STEP (Stage 3) ***\n"
        f"Open {pathway_name}_MERGED_domain_occurrences.tsv and fill:\n"
        "  function_tags   (DOF vocabulary terms)\n"
        "  binding_partner (DOT vocabulary terms)\n"
        "Then run 02_id_injection.py."
    )


# ═══════════════════════════════════════════════════════════════════════════════
# PART 8 — Unified CLI
# ═══════════════════════════════════════════════════════════════════════════════

def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="01_data_and_merge.py",
        description=(
            f"Domain-Centric Protein Ontology Pipeline  v{__version__}  —  Stage 1 + 2\n\n"
            "Stage 1: Fetches UniProtKB + InterPro data for each gene.\n"
            "Stage 2: Merges per-gene tables into pathway-level TSVs and generates\n"
            "         the domain_type_summary.tsv required by Stage 4B.\n\n"
            "MODES\n"
            "  Batch  (recommended): --gene-file FILE  or  --genes GENE ...\n"
            "  Single (debug):       --name NAME --query QUERY"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "EXAMPLES\n"
            "  # Human cAMP pathway\n"
            "  python 01_data_and_merge.py \\\n"
            "      --gene-file genes/cAMP_human.txt \\\n"
            "      --pathway-name cAMP_human --organism-id 9606 --fetch-interpro\n\n"
            "  # Mouse orthologs (same gene file, different organism)\n"
            "  python 01_data_and_merge.py \\\n"
            "      --gene-file genes/cAMP_human.txt \\\n"
            "      --pathway-name cAMP_mouse --organism-id 10090 --fetch-interpro\n\n"
            "  # Include all link tables in merged output\n"
            "  python 01_data_and_merge.py \\\n"
            "      --gene-file genes/cAMP_human.txt \\\n"
            "      --pathway-name cAMP_human --organism-id 9606 \\\n"
            "      --fetch-interpro --merge-links\n\n"
            "  # Re-run everything from scratch\n"
            "  python 01_data_and_merge.py \\\n"
            "      --gene-file genes/cAMP_human.txt \\\n"
            "      --pathway-name cAMP_human --organism-id 9606 \\\n"
            "      --fetch-interpro --force\n\n"
            "  # Single protein (debug / inspection)\n"
            "  python 01_data_and_merge.py \\\n"
            "      --name PDE4D --query \"gene:PDE4D* AND organism_id:9606\"\\\n"
            "      --fetch-interpro\n\n"
            "OUTPUT FILES (in --root-outdir / --pathway-name/)\n"
            "  *_MERGED_domain_occurrences.tsv  -> input for Stage 4A\n"
            "  *_MERGED_domain_types.tsv\n"
            "  *_domain_type_summary.tsv        -> input for Stage 4B\n"
            "  *_BATCH_STATUS.tsv\n"
            "  *_MERGED_domain_order_edges.tsv  (--merge-links only)\n"
            "  *_MERGED_occurrence_*.tsv         (--merge-links only)\n"
        ),
    )

    # ── Input: gene list (mutually exclusive) ─────────────────────────────────
    gene_group = ap.add_mutually_exclusive_group()
    gene_group.add_argument(
        "--gene-file", metavar="FILE",
        help=(
            "Path to a gene-list file (recommended for reproducibility).\n"
            "  Plain text (.txt): one gene per line; lines starting with # are comments.\n"
            "  TSV/CSV: must contain a column named 'gene'."
        ),
    )
    gene_group.add_argument(
        "--genes", nargs="+", metavar="GENE",
        help="Batch mode: space-separated gene names.",
    )

    # ── Batch options ──────────────────────────────────────────────────────────
    ap.add_argument("--pathway-name", default="",
                    help="Label used for output filenames and batch logs.")
    ap.add_argument("--root-outdir", default="",
                    help="Root output directory (default: --pathway-name).")
    ap.add_argument("--no-merge", action="store_true",
                    help="Skip Stage 2 aggregation entirely.")

    # ── Single mode options ────────────────────────────────────────────────────
    ap.add_argument("--name",   default="", help="Single mode: output file prefix.")
    ap.add_argument("--query",  default="",
                    help='Single mode: UniProtKB query, e.g. "gene:PDE4D* AND organism_id:9606".')
    ap.add_argument("--outdir", default="",
                    help="Single mode: output directory (default: <name>_out_final).")

    # ── Organism / filter ──────────────────────────────────────────────────────
    ap.add_argument(
        "--organism-id", default=DEFAULT_ORGANISM, metavar="TAXID",
        help=(
            f"NCBI taxonomy ID (default: {DEFAULT_ORGANISM} = Homo sapiens).\n"
            "Common: 9606 (human), 10090 (mouse), 10116 (rat), 7955 (zebrafish)."
        ),
    )
    ap.add_argument("--reviewed", default="yes", choices=["yes", "no", ""],
                    help="Restrict to Swiss-Prot reviewed entries (default: yes).")

    # ── InterPro options ───────────────────────────────────────────────────────
    ap.add_argument("--fetch-interpro", action="store_true",
                    help="Fetch InterPro domain occurrence data (recommended).")
    ap.add_argument("--interpro-mode", default="auto",
                    choices=["auto", "reviewed", "unreviewed", "uniprot", "both"],
                    help="InterPro API query mode (default: auto).")
    ap.add_argument("--save-interpro-json", action="store_true",
                    help="Cache raw InterPro JSON to disk (increases disk usage).")

    # ── Stage 2 aggregation options ────────────────────────────────────────────
    ap.add_argument("--merge-links", action="store_true",
                    help=(
                        "Also merge domain order edges and occurrence link tables\n"
                        "(binding, structural, signature) into pathway-level TSVs."
                    ))
    ap.add_argument("--full-accessions", action="store_true",
                    help="Write the full accession list in the domain type summary\n"
                         "(default: sample of first 50 to avoid Excel row-length limits).")
    ap.add_argument("--accession-sample-n", type=int, default=50, metavar="N",
                    help="Number of accessions to show per domain type when not --full-accessions (default: 50).")
    ap.add_argument("--union-max-items", type=int, default=200, metavar="N",
                    help="Maximum number of items in any union column in the summary (default: 200).")
    ap.add_argument("--no-summary", action="store_true",
                    help="Skip generating the domain_type_summary.tsv (Stage 4B will not work without it).")
    # --no-sequences and --seq-delay removed:
    # sequence columns are now computed in Stage 1 with no extra API calls.

    # ── Execution ──────────────────────────────────────────────────────────────
    ap.add_argument("--force", action="store_true",
                    help="Re-run all steps even if output files already exist.")
    ap.add_argument("--version", action="version",
                    version=f"01_data_and_merge.py {__version__}")

    return ap


def main(argv: Optional[List[str]] = None) -> None:
    ap   = build_arg_parser()
    args = ap.parse_args(argv)

    # ── Resolve gene list ──────────────────────────────────────────────────────
    genes: List[str] = []
    if args.gene_file:
        genes = load_gene_file(args.gene_file)
    elif args.genes:
        genes = [g.strip() for g in args.genes if g.strip()]

    # ── Batch mode ────────────────────────────────────────────────────────────
    if genes:
        pathway_name = (
            args.pathway_name.strip()
            or (Path(args.gene_file).stem if args.gene_file else "pathway")
        )
        root = Path(args.root_outdir.strip() or pathway_name)
        run_pathway_batch(
            pathway_name=pathway_name,
            genes=genes,
            reviewed=args.reviewed,
            organism_id=args.organism_id,
            root_outdir=root,
            fetch_interpro=args.fetch_interpro,
            interpro_mode=args.interpro_mode,
            save_interpro_json=args.save_interpro_json,
            force=args.force,
            do_merge=not args.no_merge,
            merge_links=args.merge_links,
            full_accessions=args.full_accessions,
            accession_sample_n=args.accession_sample_n,
            union_max_items=args.union_max_items,
            skip_summary=args.no_summary,
        )
        return

    # ── Single mode ───────────────────────────────────────────────────────────
    name  = args.name.strip()
    query = args.query.strip()
    if not name or not query:
        ap.print_help()
        raise SystemExit(
            "\n[FATAL] Provide --gene-file or --genes for batch mode,\n"
            "        or --name and --query for single-protein mode."
        )

    base_dir = Path(args.outdir.strip() or f"{name}_out_final")
    res = run_pipeline_for_gene(
        gene=name, query=query,
        reviewed=args.reviewed,
        base_dir=base_dir,
        fetch_interpro=args.fetch_interpro,
        interpro_mode=args.interpro_mode,
        save_interpro_json=args.save_interpro_json,
        force=args.force,
        pname=name,
    )
    print("\n[DONE] Outputs:")
    print(f"  Base:   {res.base_dir}")
    print(f"  Export: {res.export_dir}")
    for k, p in res.exports.items():
        print(f"   - {k}: {p.name}")


# ═══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    main()
