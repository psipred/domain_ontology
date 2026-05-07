"""
scripts/utils.py
────────────────────────────────────────────────────────────────────────────
Shared utilities for the ontology testing pipeline.
Provides: config loading, logging setup, file I/O helpers, ID normalisation,
multi-value parsing, column alias resolution, and regex validation.

All pipeline modules import from here — adapt this file if your column
names, namespaces, or file layouts change.
────────────────────────────────────────────────────────────────────────────
"""

from __future__ import annotations

import logging
import re
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

# ══════════════════════════════════════════════════════════════════════════
# Logging
# ══════════════════════════════════════════════════════════════════════════

def setup_logging(log_dir: str | Path, stage_name: str,
                  level: int = logging.INFO) -> logging.Logger:
    """
    Create a logger that writes to both stdout and a file in log_dir.
    Each stage gets its own log file.
    """
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{stage_name}.log"

    logger = logging.getLogger(stage_name)
    logger.setLevel(level)

    # Avoid duplicate handlers on re-import
    if logger.handlers:
        return logger

    fmt = logging.Formatter(
        "%(asctime)s  %(levelname)-8s  [%(name)s]  %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    # File handler
    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setFormatter(fmt)
    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(fmt)

    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


# ══════════════════════════════════════════════════════════════════════════
# Config loading
# ══════════════════════════════════════════════════════════════════════════

def load_config(config_path: str | Path) -> dict:
    """Load and return a YAML config file as a plain dict."""
    with open(config_path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    return cfg


def get_input_path(cfg: dict, key: str) -> Path | None:
    """
    Return the full path for an input file key.
    Returns None if the key is absent from config or the file doesn't exist.
    """
    fname = cfg.get(key)
    if not fname:
        return None
    p = Path(cfg.get("input_dir", "input")) / fname
    return p if p.exists() else None


def require_input(cfg: dict, key: str, logger: logging.Logger) -> Path:
    """Like get_input_path but raises FileNotFoundError if absent."""
    p = get_input_path(cfg, key)
    if p is None:
        msg = f"Required input '{key}' not found in config or on disk."
        logger.error(msg)
        raise FileNotFoundError(msg)
    return p


def make_output_dir(cfg: dict, sub: str = "") -> Path:
    """Create and return the output directory (optionally a subdirectory)."""
    base = Path(cfg.get("output_dir", "output"))
    d = base / sub if sub else base
    d.mkdir(parents=True, exist_ok=True)
    return d


# ══════════════════════════════════════════════════════════════════════════
# File loading
# ══════════════════════════════════════════════════════════════════════════

def load_table(path: str | Path, **kwargs) -> pd.DataFrame:
    """
    Load a TSV or XLSX file into a DataFrame.
    Detects format from the file extension.
    All columns are loaded as strings (dtype=str) by default for safe
    ID handling; pass dtype=None to override.
    """
    path = Path(path)
    kw = {"dtype": str, "low_memory": False}
    kw.update(kwargs)

    if path.suffix.lower() in (".xlsx", ".xls"):
        df = pd.read_excel(path, **{k: v for k, v in kw.items()
                                    if k not in ("low_memory", "sep")})
    else:
        kw.setdefault("sep", "\t")
        df = pd.read_csv(path, **kw)

    # Strip whitespace from all string columns
    for col in df.select_dtypes(include="object").columns:
        df[col] = df[col].str.strip()
    return df


def save_tsv(df: pd.DataFrame, path: str | Path,
             logger: logging.Logger | None = None) -> None:
    """Save a DataFrame to a TSV file, creating parent directories."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    if logger:
        logger.info(f"Saved {len(df)} rows → {path}")


# ══════════════════════════════════════════════════════════════════════════
# Column alias resolution
# ══════════════════════════════════════════════════════════════════════════

def resolve_column(df: pd.DataFrame, canonical: str,
                   aliases: dict[str, list[str]]) -> str | None:
    """
    Return the actual column name in df that corresponds to canonical.
    Checks the canonical name first, then all aliases.
    Returns None if nothing found.
    """
    if canonical in df.columns:
        return canonical
    for alias in aliases.get(canonical, []):
        if alias in df.columns:
            return alias
    return None


def rename_to_canonical(df: pd.DataFrame,
                        aliases: dict[str, list[str]]) -> pd.DataFrame:
    """
    Rename columns in df from any alias to the canonical name.
    Canonical names already present are left unchanged.
    """
    rename_map: dict[str, str] = {}
    for canonical, alias_list in aliases.items():
        if canonical in df.columns:
            continue
        for alias in alias_list:
            if alias in df.columns:
                rename_map[alias] = canonical
                break
    if rename_map:
        df = df.rename(columns=rename_map)
    return df


# ══════════════════════════════════════════════════════════════════════════
# Multi-value cell parsing
# ══════════════════════════════════════════════════════════════════════════

def split_cell(val: Any, delimiter: str = ";") -> list[str]:
    """
    Split a multi-value cell string by delimiter.
    Returns [] for null / empty values.
    Strips whitespace from each token.
    """
    if pd.isna(val) or str(val).strip() == "":
        return []
    return [t.strip() for t in str(val).split(delimiter) if t.strip()]


def get_delimiter(cfg: dict, col: str) -> str:
    """Look up the configured delimiter for a column, defaulting to ';'."""
    mv = cfg.get("multi_value_columns", {})
    overrides = mv.get("column_overrides", {})
    return overrides.get(col, mv.get("default_delimiter", ";"))


def explode_column(df: pd.DataFrame, col: str,
                   delimiter: str = ";") -> pd.DataFrame:
    """
    Return a DataFrame with one token per row for the given column.
    Original row index is preserved in '_original_index'.
    """
    rows = []
    for idx, val in df[col].items():
        for tok in split_cell(val, delimiter):
            rows.append({"_original_index": idx, col: tok})
    return pd.DataFrame(rows)


# ══════════════════════════════════════════════════════════════════════════
# ID normalisation
# ══════════════════════════════════════════════════════════════════════════

# Matches "DOP:000003 Cterminal" or "DOF:000012 kinase_activity" etc.
# _ID_LABEL_RE: matches 'ID description' tokens in annotation columns.
# Full database prefixes (INTERPRO:IPR..., PFAM:PF..., etc.) are included
# so normalise_id_label correctly handles them when a description is present.
_ID_LABEL_RE = re.compile(
    r'^((INTERPRO|PFAM|CDD|SMART|NCBIFAM|PRINTS|PROFILE'
    r'|DOP|DOF|DOT|GO|IPR|CHEBI|PR|MOD|R-[A-Z]{3}):[^\s]+)\s+(.+)$'
)


def normalise_id_label(value: str) -> tuple[str, str | None]:
    """
    Split a combined "ID label" token into (id_part, label_part).
    Returns (original_value, None) if no label detected.
    """
    m = _ID_LABEL_RE.match(value.strip())
    if m:
        return m.group(1), m.group(3)
    return value.strip(), None


def extract_id(value: str) -> str:
    """Extract only the ID part from a combined 'ID label' token."""
    id_part, _ = normalise_id_label(value)
    return id_part


# ══════════════════════════════════════════════════════════════════════════
# ID format validation
# ══════════════════════════════════════════════════════════════════════════

def compile_namespace_patterns(cfg: dict) -> dict[str, re.Pattern]:
    """Compile regex patterns from config into a dict of {prefix: Pattern}."""
    raw = cfg.get("namespace_patterns", {})
    return {k: re.compile(v) for k, v in raw.items()}


def validate_id_format(value: str,
                       patterns: dict[str, re.Pattern]) -> bool:
    """
    Return True if value matches any configured namespace pattern.
    Also accepts combined 'ID label' values after extracting the ID part.
    """
    id_part = extract_id(value)
    return any(p.fullmatch(id_part) for p in patterns.values())


def detect_id_prefix(value: str) -> str | None:
    """Return the namespace prefix of an ID (e.g. 'GO' from 'GO:0007165')."""
    id_part = extract_id(value)
    if ":" in id_part:
        return id_part.split(":")[0]
    if id_part.startswith("IPR"):
        return "IPR"
    if id_part.startswith("MOD"):
        return "MOD"
    return None


# ══════════════════════════════════════════════════════════════════════════
# Controlled-vocabulary lookup
# ══════════════════════════════════════════════════════════════════════════

def build_cv_lookup(ann: pd.DataFrame) -> dict[str, frozenset[str]]:
    """
    Build controlled-vocabulary lookup from annotation table.
    Returns: {term_type: frozenset of valid tokens}

    Accepts ALL of these formats for a term (e.g. DOP:000001 NearTransmembraneHelix):
      - Plain name:           "NearTransmembraneHelix"
      - Canonical name:       "near transmembrane helix"
      - Synonym:              "NearTransmembraneHelix"
      - ID alone:             "DOP:000001"
      - ID + name combined:   "DOP:000001 NearTransmembraneHelix"
      - Any of the above lowercased
    """
    from collections import defaultdict
    cv: dict[str, set] = defaultdict(set)
    for _, row in ann.iterrows():
        ttype = str(row.get("term_type", "")).strip()
        if not ttype or ttype == "nan":
            continue
        term_id   = str(row.get("ID",   "")).strip()
        term_name = str(row.get("Name", "")).strip()

        for field in ("Name", "canonical_name", "synonym"):
            val = str(row.get(field, "")).strip()
            if val and val != "nan":
                cv[ttype].add(val)
                cv[ttype].add(val.lower())

        # Also register the bare ID  (e.g. "DOP:000001")
        if term_id and term_id != "nan" and ":" in term_id:
            cv[ttype].add(term_id)
            cv[ttype].add(term_id.lower())

        # Also register the combined "ID Name" form (e.g. "DOP:000001 NearTransmembraneHelix")
        if (term_id and term_id != "nan" and ":" in term_id
                and term_name and term_name != "nan"):
            combined = f"{term_id} {term_name}"
            cv[ttype].add(combined)
            cv[ttype].add(combined.lower())

    return {k: frozenset(v) for k, v in cv.items()}


def is_cv_backed(token: str, valid_set: frozenset[str]) -> bool:
    """
    Check if a token is in the CV valid set.

    Accepts ALL of these formats:
      "NearTransmembraneHelix"            — plain name
      "DOP:000001"                        — bare ID
      "DOP:000001 NearTransmembraneHelix" — combined ID + name (both parts checked)
      Any of the above in any case

    The combined format is split into:
      id_part  = "DOP:000001"
      lbl_part = "NearTransmembraneHelix"
    and EITHER part alone is sufficient to match.
    """
    token = str(token).strip()
    if not token:
        return False
    id_part, lbl_part = normalise_id_label(token)
    candidates = {
        token,
        token.lower(),
        id_part,
        id_part.lower(),
    }
    if lbl_part:
        candidates.add(lbl_part)
        candidates.add(lbl_part.lower())
    return bool(candidates & valid_set)


# ══════════════════════════════════════════════════════════════════════════
# Miscellaneous
# ══════════════════════════════════════════════════════════════════════════

def safe_int(val: Any) -> int | None:
    """Convert value to int; return None on failure."""
    try:
        return int(float(val))
    except (ValueError, TypeError):
        return None


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    """
    Apply Benjamini-Hochberg FDR correction.
    Returns corrected p-values in the same order as input.
    Implemented without statsmodels dependency for portability;
    will use statsmodels if available.
    """
    try:
        from statsmodels.stats.multitest import multipletests
        _, corrected, _, _ = multipletests(p_values, method="fdr_bh")
        return list(corrected)
    except ImportError:
        # Manual BH
        import numpy as np
        n = len(p_values)
        if n == 0:
            return []
        order = sorted(range(n), key=lambda i: p_values[i])
        corrected = [0.0] * n
        for rank, idx in enumerate(order):
            corrected[idx] = min(1.0, p_values[idx] * n / (rank + 1))
        # Enforce monotonicity
        for i in range(n - 2, -1, -1):
            corrected[order[i]] = min(corrected[order[i]],
                                      corrected[order[i + 1]])
        return corrected


def stage_banner(logger: logging.Logger, stage_id: str,
                 title: str) -> None:
    """Print a clear stage header to the log."""
    logger.info("=" * 66)
    logger.info(f"  {stage_id}: {title}")
    logger.info("=" * 66)
