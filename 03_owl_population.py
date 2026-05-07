"""
scripts/stage4_owl_population.py
────────────────────────────────────────────────────────────────────────────
STAGE 4 — OWL POPULATION BUILDER

Converts the domain occurrence TSV into an OWL (Turtle) population file
aligned to the core ontology.  This is the bridge between the flat tabular
data (stages 1–3) and the semantic ontology (stages 5+).

What it builds
──────────────
  For every row in the occurrence table it creates three kinds of OWL
  individuals:

  Protein              — one individual per UniProt accession
  DomainType           — one individual per unique domain type ID, carrying
                         GO MF/BP/CC terms, function tags, CATH superfamilies,
                         structural coverage quality
  DomainOccurrence     — one individual per row, linking a Protein to a
                         DomainType and carrying all positional, topological,
                         proximity, binding and disorder annotations

  It also builds two cross-protein summary patterns:
  CoOccurrencePattern  — which domain types co-occur in the same protein
  DomainOrderPattern   — sequential ordering of domain types within proteins

Config keys used (pipeline_config.yaml)
────────────────────────────────────────
  domain_occurrence_table          input TSV  (required)
  owl_population_iri               IRI of the population ontology
                                   default: http://www.semanticweb.org/
                                            ontologies/camp_population
  owl_core_ns                      CORE namespace — MUST end with '#'
                                   default: http://www.semanticweb.org/
                                            ontologies/core_ontology#
  owl_import_iri                   IRI to put in owl:imports (optional)
                                   default: same as core_ns minus the '#'
  owl_occurrence_function_tags     bool — attach function_tags to
                                   DomainOccurrence as well as DomainType
                                   default: false
  owl_no_cooccurrence              bool — skip CoOccurrencePattern
                                   default: false
  owl_no_domain_order              bool — skip DomainOrderPattern
                                   default: false

Outputs (output/stage4_owl_population/)
────────────────────────────────────────
  population.ttl                   main OWL population file (Turtle)
  stage4_summary.tsv               triple count and entity count summary
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations

import re
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

import pandas as pd
from rdflib import Graph, Literal, Namespace, URIRef
from rdflib.namespace import OWL, RDF, RDFS, XSD

import sys as _sys
import os as _os
# ══════════════════════════════════════════════════════════════════════════════
# Inlined utilities (no scripts.utils needed — fully standalone)
# ══════════════════════════════════════════════════════════════════════════════
import logging
import yaml

def split_cell(val, delimiter=";"):
    if pd.isna(val) or str(val).strip() == "":
        return []
    return [t.strip() for t in str(val).split(delimiter) if t.strip()]

def get_input_path(cfg, key):
    fname = cfg.get(key)
    if not fname:
        return None
    p = Path(cfg.get("input_dir", "input")) / fname
    return p if p.exists() else None

def load_config(config_path):
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def load_table(path, **kwargs):
    path = Path(path)
    kw = {"dtype": str, "low_memory": False}
    kw.update(kwargs)
    if path.suffix.lower() in (".xlsx", ".xls"):
        df = pd.read_excel(path, **{k: v for k, v in kw.items()
                                    if k not in ("low_memory", "sep")})
    else:
        kw.setdefault("sep", "\t")
        kw.setdefault("encoding", "utf-8-sig")
        try:
            df = pd.read_csv(path, **kw)
        except UnicodeDecodeError:
            kw["encoding"] = "latin-1"
            df = pd.read_csv(path, **kw)
    for col in df.select_dtypes(include="object").columns:
        df[col] = df[col].str.strip()
    return df

def make_output_dir(cfg, sub=""):
    base = Path(cfg.get("output_dir", "output"))
    d = base / sub if sub else base
    d.mkdir(parents=True, exist_ok=True)
    return d

def rename_to_canonical(df, aliases):
    rename_map = {}
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

def save_tsv(df, path, logger=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    if logger:
        logger.info(f"Saved {len(df)} rows → {path}")

def setup_logging(log_dir, stage_name, level=logging.INFO):
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{stage_name}.log"
    logger = logging.getLogger(stage_name)
    logger.setLevel(level)
    if logger.handlers:
        return logger
    fmt = logging.Formatter(
        "%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setFormatter(fmt)
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger

def stage_banner(logger, stage_id, title):
    logger.info("=" * 66)
    logger.info(f"  {stage_id}: {title}")
    logger.info("=" * 66)




# ══════════════════════════════════════════════════════════════════════════════
# Utility helpers  (ported from new_OWL_trans.py)
# ══════════════════════════════════════════════════════════════════════════════

_SPLIT_RE  = re.compile(r"[;,]\s*")
_WS_RE     = re.compile(r"\s+")
_CHEBI_RE  = re.compile(r"(CHEBI:\d+)", re.IGNORECASE)
_DOT_RE    = re.compile(r"(DOT:\d{6})\s+([^;]+)")

OBO     = Namespace("http://purl.obolibrary.org/obo/")
UNIPROT = Namespace("http://purl.uniprot.org/uniprot/")


def _norm(x: Any) -> str:
    if x is None:
        return ""
    try:
        if pd.isna(x):
            return ""
    except Exception:
        pass
    return str(x).strip()


def _split(cell: Any) -> List[str]:
    s = _norm(cell)
    if not s:
        return []
    return [p for p in (t.strip() for t in _SPLIT_RE.split(s)) if p]


def _frag(s: Any) -> str:
    s = _norm(s)
    if not s:
        return "empty"
    s = _WS_RE.sub("_", s)
    s = re.sub(r"[^A-Za-z0-9_.\-]", "_", s)
    if re.match(r"^\d", s):
        s = "id_" + s
    return s


def _int(x: Any) -> Optional[int]:
    try:
        return int(float(_norm(x)))
    except Exception:
        return None


def _float(x: Any) -> Optional[float]:
    try:
        return float(_norm(x))
    except Exception:
        return None


def _bool(x: Any) -> Optional[bool]:
    s = _norm(x).lower()
    if s in {"1", "true", "t", "yes", "y"}:
        return True
    if s in {"0", "false", "f", "no", "n"}:
        return False
    return None


def _first(r: pd.Series, keys: Iterable[str]) -> str:
    for k in keys:
        if k in r.index:
            v = _norm(r.get(k))
            if v:
                return v
    return ""


def _chebi_ids(text: Any) -> List[str]:
    s = _norm(text)
    seen: Set[str] = set()
    out: List[str] = []
    for m in _CHEBI_RE.finditer(s):
        cid = m.group(1).upper()
        if cid not in seen:
            out.append(cid)
            seen.add(cid)
    return out


def _ordinal_idx(op: str) -> Optional[int]:
    op = _norm(op)
    if op == "N":
        return 1
    m = re.match(r"^N\+(\d+)$", op)
    if m:
        return 1 + int(m.group(1))
    return None


def _go_uri(go_id: str) -> URIRef:
    return URIRef(str(OBO) + _norm(go_id).replace(":", "_"))


def _chebi_uri(chebi_id: str) -> URIRef:
    return URIRef(str(OBO) + _norm(chebi_id).upper().replace(":", "_"))


# ══════════════════════════════════════════════════════════════════════════════
# Core graph builder  (improved version of new_OWL_trans.build_graph)
# ══════════════════════════════════════════════════════════════════════════════

def build_population_graph(
    df: pd.DataFrame,
    population_iri: str,
    core_ns: str,
    import_iri: Optional[str] = None,
    add_cooccurrence: bool = True,
    add_domain_order: bool = True,
    occurrence_function_tags: bool = False,
    logger=None,
) -> Graph:
    """
    Convert the domain occurrence DataFrame to an OWL Graph (Turtle).

    Improvements over new_OWL_trans.py
    ────────────────────────────────────
    1. Pipeline-integrated: accepts DataFrame + logger, returns Graph
    2. binding_parnter typo handled (both spellings accepted)
    3. topology_context → hasTopologyContext object property added
    4. copy_number_property column used in addition to copy_number integer
    5. Progress logging at each phase
    6. Returns entity count summary dict alongside the graph
    """
    if not core_ns.endswith("#"):
        raise ValueError("core_ns must end with '#'  e.g.  .../core_ontology#")

    POP  = Namespace(population_iri.rstrip("#") + "#")
    CORE = Namespace(core_ns)

    g = Graph()
    for prefix, ns in [
        ("rdf", RDF), ("rdfs", RDFS), ("owl", OWL), ("xsd", XSD),
        ("pop", POP), ("core", CORE), ("obo", OBO), ("uniprot", UNIPROT),
    ]:
        g.bind(prefix, ns)

    ont = URIRef(population_iri.rstrip("#"))
    g.add((ont, RDF.type, OWL.Ontology))
    if import_iri:
        g.add((ont, OWL.imports, URIRef(import_iri)))

    # ── CORE class shortcuts ───────────────────────────────────────────────
    Protein              = CORE.Protein
    DomainType           = CORE.DomainType
    DomainOccurrence     = CORE.DomainOccurrence
    GOTerm               = CORE.GOTerm
    GOFunctionTerm       = CORE.GOFunctionTerm
    GOProcessTerm        = CORE.GOProcessTerm
    GOCellComp           = CORE.GOCellularComponentTerm
    KEGGEntry            = CORE.KEGGEntry
    ReactomePathway      = CORE.ReactomePathway
    CATHSuperfamily      = CORE.CATHSuperfamily
    ChEBIEntity          = CORE.ChEBIEntity
    PositionalCategory   = CORE.PositionalCategory
    TopologyContext      = CORE.TopologyContext
    ProximityCategory    = CORE.ProximityCategory
    DisorderStatus       = CORE.DisorderStatus
    DiscontinuityPattern = CORE.DiscontinuityPattern
    InsertionCategory    = CORE.InsertionCategory
    StructCovQuality     = CORE.StructuralCoverageQuality
    CopyNumCategory      = CORE.CopyNumberCategory
    CoOccPat             = CORE.CoOccurrencePattern
    FuncCatCtx           = CORE.FunctionalCategoryContext
    DomOrderPat          = CORE.DomainOrderPattern

    # ── MINIMAL SEMANTICALLY CORRECT OBJECT PROPERTIES (9 total) ────────────
    #
    # GROUP 1 — CONTAINMENT (structural part-whole)
    #   occursInProtein   : DomainOccurrence is a structural part of Protein
    #   hasDomainType     : instance-of relationship (NOT part_of)
    #
    occursInProtein     = CORE.occursInProtein
    hasDomainType       = CORE.hasDomainType

    # GROUP 2 — CHARACTERISATION (descriptive annotation of an occurrence)
    #   hasCharacterisation : merges positional, topology, copy-number,
    #                         disorder, discontinuity — all are controlled-
    #                         vocabulary descriptors of one occurrence
    #   hasFunctionalRole   : assigned to BOTH DomainType (aggregated per
    #                         domain family) AND DomainOccurrence (per-row from
    #                         function_tags column). Same domain can carry
    #                         different roles in different proteins/contexts.
    #
    hasCharacterisation = CORE.hasCharacterisation   # NEW unified property
    hasFunctionalRole   = CORE.hasFunctionalRole     # was hasFunctionalCategory

    # GROUP 3 — SPATIAL (positional relationships between domain occurrences)
    #   isDomainNeighbour   : merges prevDomainType + nextDomainType;
    #                         direction encoded by neighbourDirection datatype
    #   hasProximityContext : renamed from hasProximityCategory
    #
    isDomainNeighbour   = CORE.isDomainNeighbour     # NEW merged property
    hasProximityContext = CORE.hasProximityContext    # was hasProximityCategory

    # GROUP 4 — CO-OCCURRENCE (cross-protein pattern, ELEVATED to DomainType)
    #   cooccursWith      : DomainType → DomainType  (symmetric)
    #   CHANGE: was on DomainOccurrence — now correctly on DomainType
    #   because co-occurrence is a TYPE-LEVEL pattern, not an occurrence fact
    #
    cooccursWith        = CORE.cooccursWith          # NEW — elevated + symmetric

    # GROUP 5 — EXTERNAL LINKS (cross-ontology references)
    #   hasGOAnnotation       : merges hasGOFunctionTerm + hasGOProcessTerm
    #                           + hasGOCellularComponentTerm; GO subtype
    #                           encoded by the GOTerm individual's own rdf:type
    #   hasBindingTarget      : renamed from hasBindingTargetChEBI
    #   hasStructuralClass    : renamed from hasCATHSuperfamily
    #
    hasGOAnnotation     = CORE.hasGOAnnotation       # NEW merged GO property
    hasBindingTarget      = CORE.hasBindingTarget       # DOT: binding target individuals
    hasBindingTargetChEBI = CORE.hasBindingTargetChEBI  # ChEBI small molecule individuals
    hasStructuralClass  = CORE.hasStructuralClass    # was hasCATHSuperfamily

    # Pattern navigation (CoOccurrencePattern / DomainOrderPattern individuals)
    domTypeA            = CORE.domainTypeA
    domTypeB            = CORE.domainTypeB
    hasDomOrderPat      = CORE.hasDomainOrderPattern
    orderDomTypeA       = CORE.orderDomainTypeA
    orderDomTypeB       = CORE.orderDomainTypeB
    hasKEGG             = CORE.hasKEGGEntry
    hasReactome         = CORE.hasReactomePathway
    hasSubcellLoc       = CORE.hasObservedSubcellularLocation

    # ── CORE datatype properties ───────────────────────────────────────────
    domainTypeId    = CORE.domainTypeId
    sourceDB        = CORE.sourceDB
    accessionProp   = CORE.accession
    externalId      = CORE.externalId
    occurrenceId    = CORE.occurrenceId
    startPos        = CORE.startPos
    endPos          = CORE.endPos
    segmentsText    = CORE.segmentsText
    relStart        = CORE.relativeStart
    relEnd          = CORE.relativeEnd
    relMid          = CORE.relativeMidpoint
    ordinalPos      = CORE.ordinalPosition
    ordinalPosIdx   = CORE.ordinalPositionIndex
    copyNumber      = CORE.copyNumber
    contCategory    = CORE.continuityCategory
    bindingTarget   = CORE.bindingTarget
    gapFromPrev     = CORE.gapFromPrev
    gapToNext       = CORE.gapToNext
    adjToPrev       = CORE.adjacentToPrev
    adjToNext       = CORE.adjacentToNext
    nearTmDist      = CORE.nearestTmDistance
    sigCleavPos     = CORE.signalCleavagePos
    nearSigDist     = CORE.nearestSignalCleavageDistance
    nearActDist     = CORE.nearestActiveSiteDistance
    subcellLocText  = CORE.proteinSubcellularLocationText
    bindSiteText    = CORE.bindingSiteTypesText
    funcTagsText    = CORE.functionTagsText
    supportProt     = CORE.supportCountProteins
    supportOcc      = CORE.supportCountOccurrences
    domOrderText    = CORE.domainOrderText

    # ── Schema declarations ────────────────────────────────────────────────
    def _class(c: URIRef) -> None:
        g.add((c, RDF.type, OWL.Class))

    def _op(p: URIRef,
            dom: Optional[URIRef] = None,
            rng: Optional[URIRef] = None) -> None:
        g.add((p, RDF.type, OWL.ObjectProperty))
        if dom:
            g.add((p, RDFS.domain, dom))
        if rng:
            g.add((p, RDFS.range, rng))

    def _dp(p: URIRef,
            dom: Optional[URIRef] = None,
            rng: Optional[URIRef] = None) -> None:
        g.add((p, RDF.type, OWL.DatatypeProperty))
        if dom:
            g.add((p, RDFS.domain, dom))
        if rng:
            g.add((p, RDFS.range, rng))

    for c in [
        Protein, DomainType, DomainOccurrence,
        GOTerm, GOFunctionTerm, GOProcessTerm, GOCellComp,
        KEGGEntry, ReactomePathway, CATHSuperfamily, ChEBIEntity,
        PositionalCategory, TopologyContext, ProximityCategory,
        DisorderStatus, DiscontinuityPattern, InsertionCategory,
        StructCovQuality, CopyNumCategory, CoOccPat, FuncCatCtx, DomOrderPat,
    ]:
        _class(c)

    # ── Object property schema declarations (minimal set) ──────────────────
    # GROUP 1: CONTAINMENT
    _op(occursInProtein,     DomainOccurrence,  Protein)
    _op(hasDomainType,       DomainOccurrence,  DomainType)

    # GROUP 2: CHARACTERISATION
    # AnnotationValue is a new superclass covering:
    # PositionalCategory, TopologyContext, CopyNumCategory,
    # DisorderStatus, DiscontinuityPattern (all controlled-vocab descriptors)
    AnnotationValue = CORE.AnnotationValue
    _class(AnnotationValue)

    # BindingTarget: named individuals for DOT: binding targets
    BindingTarget = CORE.BindingTarget
    _class(BindingTarget)
    _op(hasCharacterisation, DomainOccurrence,  AnnotationValue)
    # hasFunctionalRole is valid on DomainType (aggregated) AND
    # DomainOccurrence (per-row from function_tags column).
    _op(hasFunctionalRole,   None,              FuncCatCtx)

    # GROUP 3: SPATIAL
    _op(isDomainNeighbour,   DomainOccurrence,  DomainType)
    _op(hasProximityContext, DomainOccurrence,  ProximityCategory)

    # GROUP 4: CO-OCCURRENCE (elevated to DomainType, symmetric)
    _op(cooccursWith,        DomainType,        DomainType)
    g.add((cooccursWith,     RDF.type,          OWL.SymmetricProperty))

    # GROUP 5: EXTERNAL LINKS
    _op(hasGOAnnotation,     DomainType,        GOTerm)
    _op(hasBindingTarget,    DomainOccurrence,  BindingTarget)  # DOT: named individuals
    _op(hasBindingTargetChEBI, DomainOccurrence,  ChEBIEntity)    # ChEBI small molecules
    _op(hasStructuralClass,  DomainType,        CATHSuperfamily)

    # Pattern navigation
    _op(domTypeA,            CoOccPat,          DomainType)
    _op(domTypeB,            CoOccPat,          DomainType)
    _op(hasDomOrderPat,      DomainType,        DomOrderPat)
    _op(orderDomTypeA,       DomOrderPat,       DomainType)
    _op(orderDomTypeB,       DomOrderPat,       DomainType)
    _op(hasKEGG,             Protein,           KEGGEntry)
    _op(hasReactome,         Protein,           ReactomePathway)
    _op(hasSubcellLoc,       Protein,           None)

    # Declare AnnotationValue subclasses (all remain usable in queries)
    for _sub in [PositionalCategory, TopologyContext, CopyNumCategory,
                 DisorderStatus, DiscontinuityPattern, InsertionCategory,
                 StructCovQuality]:
        g.add((_sub, RDFS.subClassOf, AnnotationValue))

    # Datatype properties
    _dp(domainTypeId,  DomainType,      XSD.string)
    _dp(sourceDB,      DomainType,      XSD.string)
    _dp(accessionProp, DomainType,      XSD.string)
    _dp(externalId,    None,            XSD.string)
    _dp(occurrenceId,  DomainOccurrence,XSD.string)
    _dp(startPos,      DomainOccurrence,XSD.integer)
    _dp(endPos,        DomainOccurrence,XSD.integer)
    _dp(segmentsText,  DomainOccurrence,XSD.string)
    _dp(relStart,      DomainOccurrence,XSD.decimal)
    _dp(relEnd,        DomainOccurrence,XSD.decimal)
    _dp(relMid,        DomainOccurrence,XSD.decimal)
    _dp(ordinalPos,    DomainOccurrence,XSD.string)
    _dp(ordinalPosIdx, DomainOccurrence,XSD.integer)
    _dp(copyNumber,    DomainOccurrence,XSD.integer)
    _dp(contCategory,  DomainOccurrence,XSD.string)
    # bindingTarget (plain string) removed — now using object property hasBindingTarget
    _dp(gapFromPrev,   DomainOccurrence,XSD.integer)
    _dp(gapToNext,     DomainOccurrence,XSD.integer)
    _dp(adjToPrev,     DomainOccurrence,XSD.boolean)
    _dp(adjToNext,     DomainOccurrence,XSD.boolean)
    _dp(nearTmDist,    DomainOccurrence,XSD.integer)
    _dp(sigCleavPos,   DomainOccurrence,XSD.integer)
    _dp(nearSigDist,   DomainOccurrence,XSD.integer)
    _dp(nearActDist,   DomainOccurrence,XSD.integer)
    _dp(subcellLocText,Protein,         XSD.string)
    _dp(bindSiteText,  DomainOccurrence,XSD.string)
    # funcTagsText: raw semicolon-joined tags — on both DomainType
    # (aggregated) and DomainOccurrence (per-row).
    _dp(funcTagsText,  None,            XSD.string)
    _dp(supportProt,   None,            XSD.integer)
    _dp(supportOcc,    None,            XSD.integer)
    _dp(domOrderText,  Protein,         XSD.string)

    # ── Controlled-vocab helper ────────────────────────────────────────────
    # Regex to strip namespace ID prefix from combined "ID name" tokens
    # e.g. "DOF:000010 catalytic" → "catalytic"
    #      "DOP:000006 cytoplasmic" → "cytoplasmic"
    #      "DOT:000069 Ca2+"        → "Ca2+"
    _ID_LABEL_STRIP_RE = re.compile(
        r'^(?:DOT|DOP|DOF|GO|CHEBI|IPR):[^\s]+\s+(.+)$'
    )

    def _cv_display_label(raw: str) -> str:
        """Return the display label portion of a combined 'ID name' token.
        Falls back to the raw value if no ID prefix is detected.
        """
        m = _ID_LABEL_STRIP_RE.match(raw.strip())
        return m.group(1).strip() if m else raw.strip()

    _cv_cache: dict[str, URIRef] = {}

    def _cv(cls: URIRef, value: str, prefix: str) -> URIRef:
        v = _norm(value) or "unknown"
        key = f"{prefix}{v}"
        if key in _cv_cache:
            return _cv_cache[key]
        u = POP[f"{prefix}{_frag(v)}"]
        g.add((u, RDF.type, cls))
        g.add((u, RDF.type, OWL.NamedIndividual))
        # Use the bare name as RDFS.label so SPARQL label-filters work:
        #   "DOF:000010 catalytic" → label = "catalytic"  (not the full string)
        #   "cytoplasmic"          → label = "cytoplasmic" (unchanged)
        display = _cv_display_label(v)
        g.add((u, RDFS.label, Literal(display)))
        # Preserve the full ID+name form as externalId for provenance
        if display != v:
            g.add((u, CORE.externalId, Literal(v)))
        _cv_cache[key] = u
        return u

    pos_ind  = {k: _cv(PositionalCategory, k, "Pos_")
                for k in ("Nterminal", "Central", "Cterminal", "unknown")}
    copy_ind = {k: _cv(CopyNumCategory, k, "Copy_")
                for k in ("single", "repeated", "unknown")}

    # ── Pass 0: collect all proteins and domain types ──────────────────────
    proteins: Set[str] = set()
    domain_types: Set[str] = set()
    for _, r in df.iterrows():
        p  = _norm(r.get("protein_accession"))
        dt = _norm(r.get("domain_type_id"))
        if p:
            proteins.add(p)
        if dt:
            domain_types.add(dt)
        for col in ("prev_domain_type", "next_domain_type",
                    "cooccurring_domain_types"):
            for dt2 in _split(r.get(col)):
                if dt2:
                    domain_types.add(dt2)

    if logger:
        logger.info(f"  Pass 0: {len(proteins)} proteins, "
                    f"{len(domain_types)} domain types")

    # ── Pass 1: gather DomainType metadata ────────────────────────────────
    dt_label:   Dict[str, str]       = {}
    dt_db:      Dict[str, str]       = {}
    dt_acc:     Dict[str, str]       = {}
    dt_go_all:  Dict[str, Set[str]]  = defaultdict(set)
    dt_go_mf:   Dict[str, Set[str]]  = defaultdict(set)
    dt_go_bp:   Dict[str, Set[str]]  = defaultdict(set)
    dt_go_cc:   Dict[str, Set[str]]  = defaultdict(set)
    dt_func:    Dict[str, Set[str]]  = defaultdict(set)
    dt_cath:    Dict[str, Set[str]]  = defaultdict(set)
    dt_strq:    Dict[str, Set[str]]  = defaultdict(set)

    for _, r in df.iterrows():
        dt = _norm(r.get("domain_type_id"))
        if not dt:
            continue
        if dt not in dt_label:
            dt_label[dt] = _first(r, ["label", "name"]) or dt
        if dt not in dt_db and _norm(r.get("label_db")):
            dt_db[dt] = _norm(r.get("label_db"))
        if dt not in dt_acc and _norm(r.get("label_acc")):
            dt_acc[dt] = _norm(r.get("label_acc"))
        for go in _split(r.get("go_terms_all")):
            dt_go_all[dt].add(go)
        for go in _split(r.get("go_terms_mf")):
            dt_go_mf[dt].add(go)
        for go in _split(r.get("go_terms_bp")):
            dt_go_bp[dt].add(go)
        for go in _split(r.get("go_terms_cc")):
            dt_go_cc[dt].add(go)
        for tag in _split(r.get("function_tags")):
            dt_func[dt].add(tag)
        for sf in _split(r.get("cath_superfamilies")):
            dt_cath[dt].add(sf)
        if _norm(r.get("structural_coverage_quality")):
            dt_strq[dt].add(_norm(r.get("structural_coverage_quality")))

    if logger:
        logger.info(f"  Pass 1: metadata gathered for {len(dt_label)} domain types")

    # ── Proteins ───────────────────────────────────────────────────────────
    for p in sorted(proteins):
        pu = UNIPROT[_frag(p)]
        g.add((pu, RDF.type, Protein))
        g.add((pu, RDF.type, OWL.NamedIndividual))
        g.add((pu, RDFS.label, Literal(f"Protein {p}")))

    # ── DomainTypes ────────────────────────────────────────────────────────
    dt_uri: Dict[str, URIRef] = {}
    for dt in sorted(domain_types):
        if not dt:
            continue
        dtu = POP[f"DomainType_{_frag(dt)}"]
        dt_uri[dt] = dtu
        g.add((dtu, RDF.type, DomainType))
        g.add((dtu, RDF.type, OWL.NamedIndividual))
        g.add((dtu, RDFS.label, Literal(dt_label.get(dt, dt))))
        g.add((dtu, domainTypeId, Literal(dt)))
        if dt_db.get(dt):
            g.add((dtu, sourceDB,      Literal(dt_db[dt])))
        if dt_acc.get(dt):
            g.add((dtu, accessionProp, Literal(dt_acc[dt])))

        for go in sorted(dt_go_all[dt]):
            gu = _go_uri(go)
            g.add((gu, RDF.type, GOTerm)); g.add((gu, RDF.type, OWL.NamedIndividual))
            g.add((dtu, hasGOAnnotation, gu))
        for go in sorted(dt_go_mf[dt]):
            gu = _go_uri(go)
            g.add((gu, RDF.type, GOFunctionTerm)); g.add((gu, RDF.type, GOTerm))
            g.add((gu, RDF.type, OWL.NamedIndividual))
            g.add((dtu, hasGOAnnotation, gu))
        for go in sorted(dt_go_bp[dt]):
            gu = _go_uri(go)
            g.add((gu, RDF.type, GOProcessTerm)); g.add((gu, RDF.type, GOTerm))
            g.add((gu, RDF.type, OWL.NamedIndividual))
            g.add((dtu, hasGOAnnotation, gu))
        for go in sorted(dt_go_cc[dt]):
            gu = _go_uri(go)
            g.add((gu, RDF.type, GOCellComp)); g.add((gu, RDF.type, GOTerm))
            g.add((gu, RDF.type, OWL.NamedIndividual))
            g.add((dtu, hasGOAnnotation, gu))

        if dt_func[dt]:
            g.add((dtu, funcTagsText, Literal(";".join(sorted(dt_func[dt])))))
        for tag in sorted(dt_func[dt]):
            g.add((dtu, hasFunctionalRole, _cv(FuncCatCtx, tag, "FuncCat_")))

        for sf in sorted(dt_cath[dt]):
            sfu = POP[f"CATHSF_{_frag(sf)}"]
            g.add((sfu, RDF.type, CATHSuperfamily))
            g.add((sfu, RDF.type, OWL.NamedIndividual))
            g.add((sfu, RDFS.label, Literal(f"CATH superfamily {sf}")))
            g.add((sfu, externalId, Literal(sf)))
            g.add((dtu, hasStructuralClass, sfu))

        for q in sorted(dt_strq[dt]):
            g.add((dtu, hasCharacterisation, _cv(StructCovQuality, q, "StructQ_")))

    if logger:
        logger.info(f"  Built {len(dt_uri)} DomainType individuals")

    # ── DomainOccurrences ──────────────────────────────────────────────────
    occ_for_order: List[Tuple[str, int, str, str]] = []
    n_occ = 0

    # Resolve binding column (handle the original typo)
    bind_col = next(
        (c for c in ("binding_parnter", "binding_partner") if c in df.columns),
        None,
    )

    for _, r in df.iterrows():
        occ_id = _first(r, ["occurrence_id", "occ_id", "id"])
        p      = _norm(r.get("protein_accession"))
        dt     = _norm(r.get("domain_type_id"))
        if not occ_id or not p or not dt:
            continue

        occ_u  = POP[f"Occ_{_frag(occ_id)}"]
        prot_u = UNIPROT[_frag(p)]
        dtu    = dt_uri.get(dt) or POP[f"DomainType_{_frag(dt)}"]

        g.add((occ_u, RDF.type, DomainOccurrence))
        g.add((occ_u, RDF.type, OWL.NamedIndividual))
        g.add((occ_u, RDFS.label, Literal(f"DomainOccurrence {occ_id}")))
        g.add((occ_u, occurrenceId,   Literal(occ_id)))
        g.add((occ_u, occursInProtein,prot_u))
        g.add((occ_u, hasDomainType,  dtu))

        s = _int(r.get("start"))
        e = _int(r.get("end"))
        if s is not None:
            g.add((occ_u, startPos, Literal(s, datatype=XSD.integer)))
        if e is not None:
            g.add((occ_u, endPos,   Literal(e, datatype=XSD.integer)))
        if _norm(r.get("segments")):
            g.add((occ_u, segmentsText, Literal(_norm(r.get("segments")))))

        # Positional category
        pc = _norm(r.get("positional_category")) or "unknown"
        g.add((occ_u, hasCharacterisation,
               pos_ind.get(pc) or _cv(PositionalCategory, pc, "Pos_")))

        # Topology context  ← NEW: mapped as object property to TopologyContext
        for tc in _split(r.get("topology_context")):
            g.add((occ_u, hasCharacterisation, _cv(TopologyContext, tc, "Topo_")))

        for rs_val, prop in [
            (r.get("relative_start"),    relStart),
            (r.get("relative_end"),      relEnd),
            (r.get("relative_midpoint"), relMid),
        ]:
            v = _float(rs_val)
            if v is not None:
                g.add((occ_u, prop, Literal(v, datatype=XSD.decimal)))

        op = _norm(r.get("ordinal_position_group") or r.get("ordinal_position"))
        if op:
            g.add((occ_u, ordinalPos, Literal(op)))
            idx = _ordinal_idx(op)
            if idx is not None:
                g.add((occ_u, ordinalPosIdx, Literal(idx, datatype=XSD.integer)))

        cn = _int(r.get("copy_number"))
        if cn is not None:
            g.add((occ_u, copyNumber, Literal(cn, datatype=XSD.integer)))
            g.add((occ_u, hasCharacterisation,
                   copy_ind["single" if cn == 1 else "repeated"]))
        elif _norm(r.get("copy_number_property")):
            cnp = _norm(r.get("copy_number_property")).lower()
            key = "single" if "single" in cnp else "repeated" if "double" in cnp or "multi" in cnp else "unknown"
            g.add((occ_u, hasCharacterisation, copy_ind[key]))
        else:
            g.add((occ_u, hasCharacterisation, copy_ind["unknown"]))

        cc = _norm(r.get("continuity_category"))
        if cc:
            g.add((occ_u, contCategory, Literal(cc)))
            g.add((occ_u, hasCharacterisation, _cv(DiscontinuityPattern, cc, "Disc_")))

        if _norm(r.get("insertion_category")):
            g.add((occ_u, hasCharacterisation,
                   _cv(InsertionCategory, _norm(r.get("insertion_category")), "Ins_")))

        if _norm(r.get("prev_domain_type")):
            prev = _norm(r.get("prev_domain_type"))
            g.add((occ_u, isDomainNeighbour, dt_uri.get(prev) or POP[f"DomainType_{_frag(prev)}"]))  # prev
        if _norm(r.get("next_domain_type")):
            nxt = _norm(r.get("next_domain_type"))
            g.add((occ_u, isDomainNeighbour, dt_uri.get(nxt) or POP[f"DomainType_{_frag(nxt)}"]))   # next

        gfp = _int(r.get("gap_from_prev"))
        gtn = _int(r.get("gap_to_next"))
        if gfp is not None:
            g.add((occ_u, gapFromPrev, Literal(gfp, datatype=XSD.integer)))
        if gtn is not None:
            g.add((occ_u, gapToNext,   Literal(gtn, datatype=XSD.integer)))

        for bv, prop in [(_bool(r.get("adjacent_to_prev")), adjToPrev),
                          (_bool(r.get("adjacent_to_next")), adjToNext)]:
            if bv is not None:
                g.add((occ_u, prop, Literal(bv, datatype=XSD.boolean)))

        for co in _split(r.get("cooccurring_domain_types")):
            co_u = dt_uri.get(co) or POP[f"DomainType_{_frag(co)}"]
            g.add((occ_u, cooccursWith, co_u))

        if _norm(r.get("disorder_status")):
            g.add((occ_u, hasCharacterisation,
                   _cv(DisorderStatus, _norm(r.get("disorder_status")), "Dis_")))

        for pr in _split(r.get("proximity_categories")):
            g.add((occ_u, hasProximityContext, _cv(ProximityCategory, pr, "Prox_")))

        for col, prop in [
            ("nearest_tm_distance",              nearTmDist),
            ("signal_cleavage_pos",              sigCleavPos),
            ("nearest_signal_cleavage_distance", nearSigDist),
            ("nearest_active_site_distance",     nearActDist),
        ]:
            v = _int(r.get(col))
            if v is not None:
                g.add((occ_u, prop, Literal(v, datatype=XSD.integer)))

        if _norm(r.get("binding_site_types")):
            g.add((occ_u, bindSiteText, Literal(_norm(r.get("binding_site_types")))))

        # Binding target — now object property pointing to BindingTarget individuals
        # Format in TSV: "DOT:000026 ATP;DOT:000079 Mg2+" or free text without ID
        bt_raw = _norm(r.get(bind_col)) if bind_col else ""
        bt_chebi_col = _norm(r.get("binding_target_chebi"))

        if bt_raw:
            parts = [p.strip() for p in bt_raw.split(";") if p.strip()]
            for part in parts:
                m = _DOT_RE.match(part)
                if m:
                    # DOT: ID found — create/reuse a named BindingTarget individual
                    dot_id    = m.group(1)             # e.g. DOT:000026
                    dot_label = m.group(2).strip()     # e.g. ATP
                    bt_u = POP[f"BindingTarget_{_frag(dot_id)}"]
                    g.add((bt_u, RDF.type, BindingTarget))
                    g.add((bt_u, RDF.type, OWL.NamedIndividual))
                    g.add((bt_u, RDFS.label, Literal(dot_label)))
                    g.add((bt_u, CORE.externalId, Literal(dot_id)))
                    g.add((occ_u, hasBindingTarget, bt_u))
                else:
                    # No DOT: ID — store as a generic BindingTarget individual
                    safe_label = _frag(part)
                    bt_u = POP[f"BindingTarget_{safe_label}"]
                    g.add((bt_u, RDF.type, BindingTarget))
                    g.add((bt_u, RDF.type, OWL.NamedIndividual))
                    g.add((bt_u, RDFS.label, Literal(part)))
                    g.add((occ_u, hasBindingTarget, bt_u))

        # ChEBI binding targets (small molecules with CHEBI: IDs)
        for chebi in (_chebi_ids(bt_raw) + _chebi_ids(bt_chebi_col)):
            cbu = _chebi_uri(chebi)
            g.add((cbu, RDF.type, ChEBIEntity))
            g.add((cbu, RDF.type, OWL.NamedIndividual))
            g.add((occ_u, hasBindingTargetChEBI, cbu))

        # ── Occurrence-level hasFunctionalRole (always applied) ─────────
        # DESIGN CHANGE: the same InterPro domain family may carry
        # different functional roles in different proteins/contexts, so
        # hasFunctionalRole is stored on every DomainOccurrence from the
        # per-row function_tags column, not only on DomainType.
        # FuncCat_ named individuals are shared via _cv_cache — zero
        # data duplication; both DomainType and DomainOccurrence point
        # to the same individual, giving two query entry-points.
        occ_ft = r.get("function_tags")
        if _norm(occ_ft):
            _occ_tags = sorted(set(
                t.strip() for t in _split(occ_ft) if t.strip()))
            for tag in _occ_tags:
                g.add((occ_u, hasFunctionalRole, _cv(FuncCatCtx, tag, "FuncCat_")))
            g.add((occ_u, funcTagsText, Literal(";".join(_occ_tags))))

        occ_for_order.append((p, s if s is not None else 10**9, occ_id, dt))
        n_occ += 1

    # Per-protein: subcellular location, KEGG, Reactome
    prot_meta: Dict[str, Dict] = defaultdict(lambda: {
        "subcell": set(), "kegg": set(), "reactome": set()})
    for _, r in df.iterrows():
        p = _norm(r.get("protein_accession"))
        if not p:
            continue
        if _norm(r.get("protein_subcellular_location")):
            prot_meta[p]["subcell"].add(_norm(r.get("protein_subcellular_location")))
        for kid in _split(r.get("protein_kegg_ids")):
            prot_meta[p]["kegg"].add(kid)
        for rid in _split(r.get("protein_reactome_ids")):
            prot_meta[p]["reactome"].add(rid)

    for p, meta in prot_meta.items():
        pu = UNIPROT[_frag(p)]
        for loc in sorted(meta["subcell"]):
            g.add((pu, subcellLocText, Literal(loc)))
        for kid in sorted(meta["kegg"]):
            ku = POP[f"KEGG_{_frag(kid)}"]
            g.add((ku, RDF.type, KEGGEntry))
            g.add((ku, RDF.type, OWL.NamedIndividual))
            g.add((ku, externalId, Literal(kid)))
            g.add((pu, hasKEGG, ku))
        for rid in sorted(meta["reactome"]):
            ru = POP[f"Reactome_{_frag(rid)}"]
            g.add((ru, RDF.type, ReactomePathway))
            g.add((ru, RDF.type, OWL.NamedIndividual))
            g.add((ru, externalId, Literal(rid)))
            g.add((pu, hasReactome, ru))

    if logger:
        logger.info(f"  Built {n_occ} DomainOccurrence individuals")

    # ── Domain order patterns ──────────────────────────────────────────────
    if add_domain_order and occ_for_order:
        prot_to_occ: Dict[str, List] = defaultdict(list)
        for p, start, occ, dt in occ_for_order:
            prot_to_occ[p].append((start, occ, dt))

        pair_count: Dict[Tuple[str, str], int]      = defaultdict(int)
        pair_prots: Dict[Tuple[str, str], Set[str]] = defaultdict(set)

        for p, occs in prot_to_occ.items():
            types = [dt for _, _, dt in sorted(occs)]
            if types:
                g.add((UNIPROT[_frag(p)], domOrderText,
                       Literal(" > ".join(types))))
            seen: Set[Tuple[str, str]] = set()
            for a, b in zip(types, types[1:]):
                if a and b:
                    pair_count[(a, b)] += 1
                    seen.add((a, b))
            for ab in seen:
                pair_prots[ab].add(p)

        for (a, b), prots in sorted(pair_prots.items(),
                                     key=lambda x: (-len(x[1]), x[0])):
            ou = POP[f"DomOrder_{_frag(a)}__{_frag(b)}"]
            g.add((ou, RDF.type, DomOrderPat))
            g.add((ou, RDF.type, OWL.NamedIndividual))
            g.add((ou, RDFS.label, Literal(f"Domain order: {a} -> {b}")))
            g.add((ou, supportProt, Literal(len(prots),           datatype=XSD.integer)))
            g.add((ou, supportOcc,  Literal(pair_count[(a, b)],   datatype=XSD.integer)))
            a_u = dt_uri.get(a) or POP[f"DomainType_{_frag(a)}"]
            b_u = dt_uri.get(b) or POP[f"DomainType_{_frag(b)}"]
            g.add((ou, orderDomTypeA, a_u))
            g.add((ou, orderDomTypeB, b_u))
            g.add((a_u, hasDomOrderPat, ou))
            g.add((b_u, hasDomOrderPat, ou))

        if logger:
            logger.info(f"  Built {len(pair_prots)} DomainOrderPattern individuals")

    # ── Co-occurrence patterns ─────────────────────────────────────────────
    if add_cooccurrence:
        prot_types: Dict[str, Set[str]] = defaultdict(set)
        for _, r in df.iterrows():
            p  = _norm(r.get("protein_accession"))
            dt = _norm(r.get("domain_type_id"))
            if p and dt:
                prot_types[p].add(dt)

        co_prots: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
        for p, types in prot_types.items():
            for a, b in combinations(sorted(types), 2):
                co_prots[(a, b)].add(p)

        for (a, b), prots in sorted(co_prots.items(),
                                     key=lambda x: (-len(x[1]), x[0])):
            cu = POP[f"CoOcc_{_frag(a)}__{_frag(b)}"]
            g.add((cu, RDF.type, CoOccPat))
            g.add((cu, RDF.type, OWL.NamedIndividual))
            g.add((cu, RDFS.label, Literal(f"Co-occurrence: {a} + {b}")))
            g.add((cu, supportProt, Literal(len(prots), datatype=XSD.integer)))
            a_u = dt_uri.get(a) or POP[f"DomainType_{_frag(a)}"]
            b_u = dt_uri.get(b) or POP[f"DomainType_{_frag(b)}"]
            g.add((cu, domTypeA, a_u))
            g.add((cu, domTypeB, b_u))
            # Direct symmetric triples (minimal property set)
            g.add((a_u, cooccursWith, b_u))
            g.add((b_u, cooccursWith, a_u))

        if logger:
            logger.info(f"  Built {len(co_prots)} CoOccurrencePattern individuals")

    return g


# ══════════════════════════════════════════════════════════════════════════════
# Stage runner
# ══════════════════════════════════════════════════════════════════════════════

def run(cfg: dict) -> None:
    logger  = setup_logging(cfg.get("log_dir", "logs"), "stage4_owl_population")
    stage_banner(logger, "STAGE 4", "OWL Population Builder")
    out_dir = make_output_dir(cfg, "stage4_owl_population")

    # ── Load occurrence table ──────────────────────────────────────────────
    occ_path = get_input_path(cfg, "domain_occurrence_table")
    if not occ_path:
        logger.error("domain_occurrence_table not found — aborting Stage 4")
        return
    df = load_table(occ_path)
    df = rename_to_canonical(df, cfg.get("column_aliases", {}))
    logger.info(f"Occurrence table: {len(df)} rows, "
                f"{df['protein_accession'].nunique()} proteins")

    # ── OWL config ────────────────────────────────────────────────────────
    pop_iri    = cfg.get(
        "owl_population_iri",
        "http://www.semanticweb.org/ontologies/camp_population")
    core_ns    = cfg.get(
        "owl_core_ns",
        "http://www.semanticweb.org/32045/ontologies/2026/0/core_ontology#")
    import_iri = cfg.get(
        "owl_import_iri",
        core_ns.rstrip("#"))

    occ_func   = cfg.get("owl_occurrence_function_tags", False)
    no_coocc   = cfg.get("owl_no_cooccurrence",  False)
    no_order   = cfg.get("owl_no_domain_order",  False)

    logger.info(f"Population IRI: {pop_iri}")
    logger.info(f"Core namespace: {core_ns}")
    logger.info(f"Import IRI:     {import_iri}")
    logger.info(f"CoOccurrence:   {'disabled' if no_coocc else 'enabled'}")
    logger.info(f"DomainOrder:    {'disabled' if no_order  else 'enabled'}")
    if occ_func:
        logger.warning("occurrence_function_tags=True — may cause OWL domain "
                       "conflicts if hasFunctionalCategory is DomainType-only")

    # ── Build graph ────────────────────────────────────────────────────────
    logger.info("Building OWL graph …")
    g = build_population_graph(
        df               = df,
        population_iri   = pop_iri,
        core_ns          = core_ns,
        import_iri       = import_iri,
        add_cooccurrence = not no_coocc,
        add_domain_order = not no_order,
        occurrence_function_tags = occ_func,
        logger           = logger,
    )

    # ── Serialise ──────────────────────────────────────────────────────────
    ttl_path = out_dir / "population.ttl"
    g.serialize(destination=str(ttl_path), format="turtle")
    n_triples = len(g)
    logger.info(f"Serialised → {ttl_path}  ({n_triples:,} triples)")

    # ── Summary ────────────────────────────────────────────────────────────
    from rdflib.namespace import OWL as _OWL
    n_ind = sum(1 for _ in g.subjects(RDF.type, _OWL.NamedIndividual))
    summary = pd.DataFrame([{
        "metric": "triples",           "value": n_triples},
        {"metric": "named_individuals","value": n_ind},
        {"metric": "input_rows",       "value": len(df)},
        {"metric": "proteins",         "value": df["protein_accession"].nunique()},
        {"metric": "domain_types",     "value": df["domain_type_id"].nunique()},
        {"metric": "population_iri",   "value": pop_iri},
        {"metric": "core_ns",          "value": core_ns},
        {"metric": "output_ttl",       "value": str(ttl_path)},
    ])
    save_tsv(summary, out_dir / "stage4_summary.tsv", logger)

    logger.info("")
    logger.info("=" * 60)
    logger.info("  STAGE 4 COMPLETE")
    logger.info(f"  Triples:           {n_triples:,}")
    logger.info(f"  Named individuals: {n_ind:,}")
    logger.info(f"  Output:            {ttl_path}")
    logger.info("=" * 60)


# ══════════════════════════════════════════════════════════════════════════════
# CLI  (keeps backward compatibility with new_OWL_trans.py command-line usage)
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(
        description="Stage 4 — OWL Population Builder  "
                    "(pipeline mode: python run_pipeline.py --stage 4)")
    ap.add_argument("--config", default="config/pipeline_config.yaml",
                    help="Path to pipeline_config.yaml")

    # Direct mode (mirrors new_OWL_trans.py for standalone use)
    ap.add_argument("--in",   dest="inp",            default="",
                    help="Direct input TSV (bypasses config)")
    ap.add_argument("--out",  dest="out",            default="",
                    help="Direct output .ttl path")
    ap.add_argument("--population-iri",              default="",
                    help="Population ontology IRI")
    ap.add_argument("--core-ns",                     default="",
                    help="Core namespace (must end with '#')")
    ap.add_argument("--import-iri",                  default="",
                    help="owl:imports IRI")
    ap.add_argument("--no-global-cooccurrence",      action="store_true")
    ap.add_argument("--no-global-domain-order",      action="store_true")
    ap.add_argument("--occurrence-function-tags",    action="store_true")
    args = ap.parse_args()

    # Show received arguments so user can see what was parsed
    print("[Stage 4] Arguments received:")
    print(f"  --in             : {args.inp}")
    print(f"  --out            : {args.out}")
    print(f"  --population-iri : {args.population_iri}")
    print(f"  --core-ns        : {args.core_ns}")
    print()

    if args.inp and args.out and args.population_iri and args.core_ns:
        import pandas as _pd
        print("[Stage 4] Loading input TSV ...")
        try:
            _df = _pd.read_csv(args.inp, sep="\t", dtype=str, encoding="utf-8-sig").fillna("")
        except UnicodeDecodeError:
            _df = _pd.read_csv(args.inp, sep="\t", dtype=str, encoding="latin-1").fillna("")
        print(f"[Stage 4] Loaded {len(_df)} rows")
        print("[Stage 4] Building OWL graph ...")
        _g = build_population_graph(
            _df,
            population_iri   = args.population_iri,
            core_ns          = args.core_ns,
            import_iri       = (args.import_iri.strip() or None),
            add_cooccurrence = not args.no_global_cooccurrence,
            add_domain_order = not args.no_global_domain_order,
            occurrence_function_tags = args.occurrence_function_tags,
        )
        print("[Stage 4] Serialising to Turtle ...")
        _g.serialize(destination=args.out, format="turtle")
        print(f"[DONE] Wrote {args.out}  ({len(_g):,} triples)")
    else:
        print("[Stage 4] ERROR: one or more required arguments are empty.")
        print("  Make sure all four flags are provided:")
        print("    --in  --out  --population-iri  --core-ns")
