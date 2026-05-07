"""
scripts/stage3_shacl.py  — SHACL Validation (pure-Python, no rdflib needed)
Parses instance_data.ttl by splitting on blank lines, then validates each
subject block against the 8 sh:property constraints in shapes.ttl.
Falls back to pyshacl/rdflib if available.
"""
from __future__ import annotations
import re
from pathlib import Path
import pandas as pd
import sys as _sys
import os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from scripts.utils import (get_input_path, load_config, make_output_dir,
                            save_tsv, setup_logging, stage_banner)

# Mirror of shapes.ttl property constraints
# (shape, target_class, pred_local, min_count, max_count, datatype, min_val, max_val, message)
# Expanded to 18 constraints across 5 shapes matching shapes.ttl v2.
SHAPES = [
    # ── DomainOccurrence (S1–S12) ──────────────────────────────────────────
    ("DomainOccurrenceShape","DomainOccurrence","occursInProtein",    1, 1,    None,      None, None,
     "S1: must link to exactly one Protein via occursInProtein"),
    ("DomainOccurrenceShape","DomainOccurrence","hasDomainType",      1, 1,    None,      None, None,
     "S2: must link to exactly one DomainType via hasDomainType"),
    ("DomainOccurrenceShape","DomainOccurrence","startPos",           1, 1,    "integer", 1,    None,
     "S3: startPos required, must be an integer >= 1"),
    ("DomainOccurrenceShape","DomainOccurrence","endPos",             1, 1,    "integer", 1,    None,
     "S4: endPos required, must be an integer >= 1"),
    ("DomainOccurrenceShape","DomainOccurrence","relativeStart",      0, 1,    "decimal", 0.0,  1.0,
     "S5: relativeStart must be a decimal in [0.0, 1.0]"),
    ("DomainOccurrenceShape","DomainOccurrence","relativeEnd",        0, 1,    "decimal", 0.0,  1.0,
     "S6: relativeEnd must be a decimal in [0.0, 1.0]"),
    ("DomainOccurrenceShape","DomainOccurrence","isRepresentative",   0, 1,    "boolean", None, None,
     "S7: isRepresentative must be a boolean"),
    ("DomainOccurrenceShape","DomainOccurrence","copyNumberProperty", 1, 1,    None,      None, None,
     "S8: copyNumberProperty must be present exactly once"),
    ("DomainOccurrenceShape","DomainOccurrence","hasPositionalCategory",1,1,   None,      None, None,
     "S9: must have exactly one PositionalCategory"),
    ("DomainOccurrenceShape","DomainOccurrence","domainEvidenceLevel",1, 1,    None,      None, None,
     "S10: domainEvidenceLevel must be present"),
    ("DomainOccurrenceShape","DomainOccurrence","hasTopologyContext",  1, 1,   None,      None, None,
     "S11: must have exactly one TopologyContext"),
    ("DomainOccurrenceShape","DomainOccurrence","externalId",         1, 1,    None,      None, None,
     "S12: externalId (occurrence key) must be present"),
    # ── Protein (S13–S14) ──────────────────────────────────────────────────
    ("ProteinShape",          "Protein",         "externalId",         1, None, None,      None, None,
     "S13: Protein must have at least one externalId"),
    ("ProteinShape",          "Protein",         "label",              1, None, None,      None, None,
     "S14: Protein must have a rdfs:label"),
    # ── DomainType (S15) ──────────────────────────────────────────────────
    # Note: DomainTypes use rdfs:label as their primary identifier (not externalId)
    ("DomainTypeShape",       "DomainType",      "label",              1, None, None,      None, None,
     "S15: DomainType must have a rdfs:label"),
    # ── FunctionalTag (S17) ────────────────────────────────────────────────
    ("FunctionalTagShape",    "FunctionalTag",   "label",              1, None, None,      None, None,
     "S17: FunctionalTag must have a rdfs:label"),
    # ── ReactomePathway (S18) ──────────────────────────────────────────────
    ("ReactomePathwayShape",  "ReactomePathway", "externalId",         1, None, None,      None, None,
     "S18: ReactomePathway must have an externalId"),
]

# ── Turtle block parser ───────────────────────────────────────────────────────
_PRED_RE = re.compile(
    r'^\s+((?:ex|rdfs|owl|rdf|xsd):[\w\-]+|a)\s+(.+?)\s*[;.]?\s*$'
)

def _parse_ttl(path: Path) -> dict[str, dict[str, list[str]]]:
    """Split on blank lines → one block per subject. Parse pred-obj per line."""
    text = path.read_text(encoding="utf-8", errors="replace")
    # Remove comment lines
    text = re.sub(r'(?m)^#[^\n]*\n', '', text)
    # Replace \r\n → \n
    text = text.replace('\r\n', '\n').replace('\r', '\n')

    subjects: dict[str, dict[str, list[str]]] = {}

    for block in re.split(r'\n{2,}', text):
        block = block.strip()
        if not block:
            continue
        lines = block.splitlines()
        if not lines:
            continue
        # Subject = first token of first line
        first = lines[0].strip()
        subj_m = re.match(r'^(ex:\S+|<[^>]+>)', first)
        if not subj_m:
            continue
        subj = subj_m.group(1)
        props: dict[str, list[str]] = {}
        subjects[subj] = props

        # Parse all lines in block for pred-obj
        for line in lines:
            line_s = line.strip().rstrip(';. ')
            # Match "pred obj"  (indented or first line after subject)
            # Handle first line: "subject pred obj"
            m = re.match(
                r'^(?:' + re.escape(subj) + r'\s+)?'
                r'((?:ex|rdfs|owl|rdf|xsd):[\w\-]+|a)\s+(.+)$',
                line_s
            )
            if not m:
                continue
            pred_full, obj = m.group(1).strip(), m.group(2).strip().rstrip(';.')
            # Local name of predicate
            pred_local = pred_full.split(':')[-1] if ':' in pred_full else pred_full
            # Clean triple-quoted strings: """value""" → value
            obj = re.sub(r'"{3}(.*?)"{3}', r'\1', obj, flags=re.DOTALL)
            obj = obj.strip('"\'').strip()
            props.setdefault(pred_local, []).append(obj)

    return subjects


def _get_class(props: dict[str, list[str]]) -> str | None:
    for v in props.get("a", []):
        if v.startswith("ex:"):
            return v[3:]
        if v.startswith("owl:") or v == "owl:Ontology":
            return None
    return None


def _is_int(v: str) -> bool:
    try:
        int(v.split("^^")[0].strip('" '))
        return True
    except (ValueError, TypeError):
        return False


def _is_bool(v: str) -> bool:
    return v.strip().lower().split("^^")[0].strip('" ') in ("true", "false")


def _validate(subjects: dict[str, dict[str, list[str]]]) -> list[dict]:
    viols: list[dict] = []
    for subj, props in subjects.items():
        cls = _get_class(props)
        if not cls:
            continue
        for (shape, target, pred, min_c, max_c, dtype, min_val, max_val, msg) in SHAPES:
            if target != cls:
                continue
            vals = props.get(pred, [])
            n = len(vals)
            base = dict(shape=shape, subject=subj, property=pred, message=msg)
            # MinCount
            if min_c is not None and n < min_c:
                viols.append({**base, "violation_type": "MinCountViolation",
                               "detail": f"count={n} < minCount={min_c}"})
                continue
            # MaxCount
            if max_c is not None and n > max_c:
                viols.append({**base, "violation_type": "MaxCountViolation",
                               "detail": f"count={n} > maxCount={max_c}"})
            # Datatype
            if dtype and vals:
                for v in vals:
                    if dtype == "integer" and not _is_int(v):
                        viols.append({**base, "violation_type": "DatatypeViolation",
                                       "detail": f"not xsd:integer: {v!r}"})
                    elif dtype == "boolean" and not _is_bool(v):
                        viols.append({**base, "violation_type": "DatatypeViolation",
                                       "detail": f"not xsd:boolean: {v!r}"})
                    elif dtype == "decimal":
                        try:
                            float(v.split("^^")[0].strip('" '))
                        except (ValueError, TypeError):
                            viols.append({**base, "violation_type": "DatatypeViolation",
                                           "detail": f"not xsd:decimal: {v!r}"})
            # MinInclusive / MaxInclusive (integer and decimal)
            if (min_val is not None or max_val is not None) and dtype in ("integer","decimal") and vals:
                for v in vals:
                    try:
                        x = float(v.split("^^")[0].strip('" '))
                        if min_val is not None and x < min_val:
                            viols.append({**base, "violation_type": "MinInclusiveViolation",
                                           "detail": f"value={x} < minInclusive={min_val}"})
                        if max_val is not None and x > max_val:
                            viols.append({**base, "violation_type": "MaxInclusiveViolation",
                                           "detail": f"value={x} > maxInclusive={max_val}"})
                    except (ValueError, TypeError):
                        pass
    return viols


def run(cfg: dict) -> None:
    logger = setup_logging(cfg.get("log_dir", "logs"), "stage3_shacl")
    stage_banner(logger, "STAGE 3", "SHACL Validation")
    out_dir = make_output_dir(cfg, "stage3_shacl")

    data_path   = get_input_path(cfg, "instance_data")
    shapes_path = get_input_path(cfg, "shacl_shapes")
    if not data_path:
        logger.warning("instance_data.ttl not found — skipping Stage 3"); return
    if not shapes_path:
        logger.warning("shapes.ttl not found — skipping Stage 3"); return

    # Try pyshacl first
    try:
        import pyshacl
        from rdflib import Graph, Namespace, RDF
        SH = Namespace("http://www.w3.org/ns/shacl#")

        logger.info("rdflib + pyshacl available — running full SHACL validation")
        dg = Graph().parse(str(data_path))
        sg = Graph().parse(str(shapes_path))
        conforms, rg, rt = pyshacl.validate(dg, shacl_graph=sg,
                                              inference="rdfs", abort_on_first=False)
        # Write the human-readable report
        (out_dir / "shacl_report.txt").write_text(rt, encoding="utf-8")
        status = "CONFORMS" if conforms else "VIOLATIONS_FOUND"
        logger.info(f"pyshacl: {status}")

        # ── Parse violations from the results graph ──────────────────────
        viol_rows = []
        shape_viols: dict[str, int] = {}
        for result in rg.subjects(RDF.type, SH.ValidationResult):
            focus     = str(rg.value(result, SH.focusNode)   or "")
            path      = str(rg.value(result, SH.resultPath)  or "")
            message   = str(rg.value(result, SH.resultMessage) or "")
            severity  = str(rg.value(result, SH.resultSeverity) or "")
            src_shape = str(rg.value(result, SH.sourceShape) or "")
            src_constraint = str(rg.value(result, SH.sourceConstraintComponent) or "")
            # Shorten URIs for readability
            def _short(uri: str) -> str:
                for prefix in ["http://www.w3.org/ns/shacl#",
                                "http://example.org/domain-ontology-core#",
                                "http://www.w3.org/2001/XMLSchema#"]:
                    uri = uri.replace(prefix, "")
                return uri
            shape_name = _short(src_shape).split("/")[-1]
            vtype      = _short(src_constraint).replace("ConstraintComponent", "Violation")
            viol_rows.append({
                "shape":          shape_name,
                "subject":        _short(focus),
                "property":       _short(path),
                "violation_type": vtype,
                "detail":         message,
                "message":        message,
                "severity":       _short(severity),
            })
            shape_viols[shape_name] = shape_viols.get(shape_name, 0) + 1

        vdf = (pd.DataFrame(viol_rows) if viol_rows else
               pd.DataFrame(columns=["shape","subject","property",
                                      "violation_type","detail","message","severity"]))
        save_tsv(vdf, out_dir / "shacl_violations.tsv", logger)

        # violation summary by shape
        vsummary = pd.DataFrame([{"shape": s, "n_violations": n}
                                  for s, n in sorted(shape_viols.items())])
        save_tsv(vsummary, out_dir / "shacl_violation_summary.tsv", logger)

        # shape summary (PASS / FAIL per shape)
        all_shapes = sorted({s[0] for s in SHAPES})
        shape_summary_rows = []
        for sn in all_shapes:
            n = shape_viols.get(sn, 0)
            shape_summary_rows.append({
                "shape":  sn,
                "status": "FAIL" if n > 0 else "PASS",
                "n_violations": n,
            })
        sdf = pd.DataFrame(shape_summary_rows)
        save_tsv(sdf, out_dir / "shacl_shape_summary.tsv", logger)

        n_pass = sum(1 for r in shape_summary_rows if r["status"] == "PASS")
        n_fail = sum(1 for r in shape_summary_rows if r["status"] == "FAIL")

        # Write rich text report
        lines = ["=" * 66,
                 "  STAGE 3 — SHACL VALIDATION REPORT (pyshacl)",
                 "=" * 66,
                 f"  Source:     {data_path}",
                 f"  Shapes:     {shapes_path}",
                 f"  Violations: {len(viol_rows)}   "
                 f"Shapes PASS={n_pass}  FAIL={n_fail}", ""]
        for r in shape_summary_rows:
            tag = "[PASS]" if r["status"] == "PASS" else "[FAIL]"
            lines.append(f"  {tag} {r['shape']:<40} violations={r['n_violations']}")
        lines += ["", rt]
        (out_dir / "shacl_report.txt").write_text("\n".join(lines), encoding="utf-8")

        for r in shape_summary_rows:
            tag = r["status"]
            fn  = logger.warning if tag == "FAIL" else logger.info
            fn(f"  [{tag}] {r['shape']}  violations={r['n_violations']}")

        logger.info(f"Stage 3 complete — {len(viol_rows)} violations  "
                    f"({n_pass} shapes pass, {n_fail} fail)")
        return

    except ImportError:
        logger.info("pyshacl/rdflib not installed — using pure-Python validator")
    except Exception as exc:
        logger.warning(f"pyshacl error ({exc}) — falling back to pure-Python")

    # Pure-Python path
    logger.info(f"Parsing {data_path} …")
    subjects = _parse_ttl(Path(data_path))

    cls_counts: dict[str, int] = {}
    for p in subjects.values():
        c = _get_class(p)
        if c:
            cls_counts[c] = cls_counts.get(c, 0) + 1
    logger.info("Instance counts: " +
                ", ".join(f"{k}={v}" for k, v in sorted(cls_counts.items())))

    viols = _validate(subjects)
    vdf = (pd.DataFrame(viols) if viols else
           pd.DataFrame(columns=["shape","subject","property",
                                  "violation_type","detail","message"]))
    save_tsv(vdf, out_dir / "shacl_violations.tsv", logger)

    summ = (vdf.groupby(["shape","property","violation_type"])
               .size().reset_index(name="count")
            if not vdf.empty else
            pd.DataFrame(columns=["shape","property","violation_type","count"]))
    save_tsv(summ, out_dir / "shacl_violation_summary.tsv", logger)

    tmap = {s[0]: s[1] for s in SHAPES}
    shape_rows = []
    for sn in sorted({s[0] for s in SHAPES}):
        n = int((vdf["shape"] == sn).sum()) if not vdf.empty else 0
        shape_rows.append({"shape": sn, "target_class": tmap.get(sn,""),
                            "n_violations": n, "status": "PASS" if n == 0 else "FAIL"})
    save_tsv(pd.DataFrame(shape_rows), out_dir / "shacl_shape_summary.tsv", logger)

    n_pass = sum(1 for r in shape_rows if r["status"] == "PASS")
    n_fail = sum(1 for r in shape_rows if r["status"] == "FAIL")
    report = ["=" * 66,
              "  STAGE 3 — SHACL VALIDATION REPORT (pure-Python)", "=" * 66,
              f"  Source:    {data_path}",
              f"  Shapes:    {shapes_path}",
              f"  Instances: " + ", ".join(f"{k}={v}" for k,v in sorted(cls_counts.items())),
              f"  Violations: {len(vdf)}   Shapes PASS={n_pass}  FAIL={n_fail}", ""]
    for r in shape_rows:
        report.append(f"  [{r['status']:4s}] {r['shape']:30s}  violations={r['n_violations']}")
    if not vdf.empty:
        report += ["", "  Sample violations (first 20):"]
        for _, r in vdf.head(20).iterrows():
            report.append(f"    {r.get('shape','')}  {r.get('property','')}  "
                          f"{r.get('violation_type','')}  {r.get('detail','')}")
    (out_dir / "shacl_report.txt").write_text("\n".join(report) + "\n", encoding="utf-8")

    for r in shape_rows:
        fn = logger.warning if r["status"] == "FAIL" else logger.info
        fn(f"  [{r['status']}] {r['shape']}  violations={r['n_violations']}")
    logger.info(f"Stage 3 complete — {len(vdf)} violations  "
                f"({n_pass} shapes pass, {n_fail} fail)")


if __name__ == "__main__":
    import sys, os
    # Always run from the ontology_pipeline project root so that
    # relative paths (input/, output/, config/) resolve correctly
    # regardless of where PyCharm / the terminal launched the script from.
    _project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(_project_root)
    _cfg_path = sys.argv[1] if len(sys.argv) > 1 else "config/pipeline_config.yaml"
    run(load_config(_cfg_path))
