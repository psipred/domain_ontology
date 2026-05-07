#!/usr/bin/env python3
"""
two_pathway_transfer_validation_pipeline.py.
Two-Pathway Ontology Transfer Validation Pipeline

Purpose
Validates two ontology-aligned protein domain pathway tables against the same
annotation registry, measures how well the ontology framework transfers from
the reference (cAMP) pathway to a second pathway, and computes biological and
structural overlap between the two datasets.

Designed to be extended later for:
  • human cAMP  →  human other pathway
  • human cAMP  →  other species cAMP
  • human cAMP  →  other pathway in another species

Usage
  python two_pathway_transfer_validation_pipeline.py \
      --camp   cAMP_human_mapped.tsv \
      --other  otherPathway_human_mapped.tsv \
      --ann    annotation_table.xlsx \
      --out    output_dir \
      [--plots]

Required inputs
  cAMP_human_mapped.tsv       — reference pathway occurrence table
  otherPathway_human_mapped.tsv — second pathway occurrence table
  annotation_table.xlsx       — ontology annotation registry (DOP/DOF/DOT)

Outputs
───────
  validation_report_cAMP.tsv       per-check results for pathway 1
  validation_report_otherPathway.tsv per-check results for pathway 2
  invalid_rows.tsv                 every flagged row from both pathways
  unknown_terms_report.tsv         tokens not found in any term_type
  missing_term_report.tsv          unmatched CV tokens with closest-match hint
  pathway_summary_metrics.tsv      per-pathway counts and coverage stats
  transfer_metrics.tsv             schema/vocab/rule portability scores
  overlap_metrics.tsv              Jaccard and set overlap per category
  comparison_report.tsv            side-by-side biological comparison

Optional (--plots)
  term_reuse_barplot.png
  jaccard_heatmap.png
  validation_summary_plot.png
"""

from __future__ import annotations
import argparse, re, sys
from collections import Counter, defaultdict
from difflib import get_close_matches
from pathlib import Path
from typing import Any
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

plt.rcParams.update({
    "font.family":"DejaVu Sans","font.size":8,"axes.titlesize":9,"axes.labelsize":8,
    "xtick.labelsize":7,"ytick.labelsize":7,"legend.fontsize":7,
    "axes.linewidth":0.7,"axes.spines.top":False,"axes.spines.right":False,
    "xtick.major.width":0.7,"ytick.major.width":0.7,"xtick.major.size":3,
    "ytick.major.size":3,"xtick.direction":"out","ytick.direction":"out",
    "legend.frameon":False,"figure.dpi":300,"savefig.dpi":300,"savefig.bbox":"tight",
    "axes.grid":False,"figure.facecolor":"white","axes.facecolor":"white",
    "pdf.fonttype":42,"ps.fonttype":42,
})
PAL = {
    "blue":"#0077BB","orange":"#EE7733","green":"#009988","red":"#CC3311",
    "purple":"#AA3377","yellow":"#CCBB44","cyan":"#33BBEE","grey":"#BBBBBB",
    "darkgrey":"#555555","lightgrey":"#EEEEEE","pass":"#44AA99","warn":"#DDAA33",
    "fail":"#BB5566","shared":"#44AA99","ref_only":"#0077BB","tgt_only":"#EE7733",
    "unused":"#BBBBBB",
}

RE_GO        = re.compile(r"^GO:\d{7}$")
RE_REACTOME  = re.compile(r"^R-[A-Z]{3}-\d+$")
RE_KEGG      = re.compile(r"^[a-z]{2,5}:\d+$")
RE_CHEBI_OK  = re.compile(r"^CHEBI:\d+$")
RE_CHEBI_DBL = re.compile(r"^CHEBI:CHEBI:\d+$")
RE_UNIPROT   = re.compile(r"^[A-Z][0-9][A-Z0-9]{3}[0-9]([A-Z][0-9][A-Z0-9]{3}[0-9])?$")
RE_DT_ID     = re.compile(r"^(CDD:cd\d+|INTERPRO:IPR\d{6}|NCBIFAM:TIGR\d{5}|PFAM:PF\d{5}|PRINTS:PR\d{5}|PROFILE:PS\d{5}|SMART:SM\d{5}|SUPERFAMILY:SSF\d+|PANTHER:PTHR\d+|GENE3D:G3DSA:[\d.]+)$")
RE_ID_LABEL  = re.compile(r"^((DOP|DOF|DOT|GO|IPR|CHEBI|PR|MOD|R-[A-Z]{3}):[^\s]+)\s+(.+)$")

NAMESPACE_RULES = [
    ("positional_category","DOP","positional_category must use DOP terms"),
    ("topology_context","DOP","topology_context must use DOP terms"),
    ("proximity_categories","DOP","proximity_categories must use DOP terms"),
    ("copy_number_property","DOP","copy_number_property must use DOP terms"),
    ("function_tags","DOF","function_tags must use DOF terms"),
    ("binding_partner","DOT","binding_partner must use DOT terms"),
    ("binding_parnter","DOT","binding_parnter must use DOT terms"),
]
CV_COL_TO_TERM_TYPE = {
    "positional_category":"positional_or_topology_term",
    "topology_context":"positional_or_topology_term",
    "proximity_categories":"positional_or_topology_term",
    "copy_number_property":"CopyNumberCategory",
    "function_tags":"functional_term",
    "binding_partner":"binding_target_term",
    "binding_parnter":"binding_target_term",
}
# FIX-2: tiered required columns
REQUIRED_CORE       = ["occurrence_id","protein_accession","domain_type_id","start","end"]
REQUIRED_ANNOTATION = ["positional_category","topology_context","function_tags"]
REQUIRED_OPTIONAL   = ["copy_number_property"]  # not yet in all pipeline versions
REQUIRED_COLUMNS    = REQUIRED_CORE + REQUIRED_ANNOTATION + REQUIRED_OPTIONAL
GO_COLS = ["go_terms_all","go_terms_mf","go_terms_bp","go_terms_cc"]


# ── I/O ──────────────────────────────────────────────────────────────────────
def load_table(path):
    path = Path(path)
    if not path.exists(): raise FileNotFoundError(f"Not found: {path}")
    kw = {"dtype":str}
    df = pd.read_excel(path,**kw) if path.suffix.lower() in (".xlsx",".xls") else pd.read_csv(path,sep="\t",low_memory=False,**kw)
    for col in df.select_dtypes(include=["object","string"]).columns:
        df[col] = df[col].str.strip()
    return df

def save_tsv(df,path,label=""):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    df.to_csv(path,sep="\t",index=False)
    print(f"  ✓  Saved {len(df):>5} rows → {path}{' ['+label+']' if label else ''}")

def print_section(t): print(f"\n{'─'*68}\n  {t}\n{'─'*68}")
def print_check(cid,st,n,desc):
    icon={"PASS":"✓","FAIL":"✗","WARN":"⚠"}.get(st,"?")
    suf={"FAIL":"  ← ACTION REQUIRED","WARN":"  ← review"}.get(st,"")
    print(f"  {icon} [{cid}] {st:4}  n={n:<5}  {desc}{suf}")
def _panel_label(ax,letter,x=-0.10,y=1.06):
    ax.text(x,y,letter,transform=ax.transAxes,fontsize=11,fontweight="bold",va="top",ha="left")


# ── Token helpers ─────────────────────────────────────────────────────────────
def split_cell(val,delim=";"):
    if pd.isna(val) or str(val).strip()=="": return []
    return [t.strip() for t in str(val).split(delim) if t.strip()]

def extract_id(tok):
    m=RE_ID_LABEL.match(tok.strip()); return m.group(1) if m else tok.strip()

def extract_label(tok):
    m=RE_ID_LABEL.match(tok.strip()); return m.group(3) if m else tok.strip()

def get_id_namespace(tok):
    ip=extract_id(tok)
    for p in ("DOP","DOF","DOT"):
        if ip.startswith(p+":"): return p
    return None

def build_cv_lookup(ann):
    cv=defaultdict(set)
    for _,row in ann.iterrows():
        tt=str(row.get("term_type","")).strip()
        if not tt or tt=="nan": continue
        for f in ("ID","Name","canonical_name","synonym"):
            v=str(row.get(f,"")).strip()
            if v and v!="nan": cv[tt].add(v); cv[tt].add(v.lower())
    return {k:frozenset(v) for k,v in cv.items()}

def is_cv_backed(tok,valid):
    t=tok.strip(); s=re.sub(r":\d+$","",t)
    return any(v in valid for v in (t,t.lower(),extract_id(t),extract_id(t).lower(),extract_label(t),extract_label(t).lower(),s,s.lower()))

def all_valid_tokens(cv):
    return [v for s in cv.values() for v in s if not v.islower()]

def _normalise_chebi(val):
    if not isinstance(val,str) or "CHEBI:CHEBI:" not in val: return val,0
    return re.sub(r"CHEBI:CHEBI:","CHEBI:",val), val.count("CHEBI:CHEBI:")


# ── Validation ────────────────────────────────────────────────────────────────
def _row_dict(df,idx):
    return {k:(str(v)[:300] if isinstance(v,str) else v) for k,v in df.loc[idx].to_dict().items()}

def _detect_binding_col(df,override=None):
    if override and override in df.columns: return override
    for c in ("binding_partner","binding_parnter"):
        if c in df.columns: return c
    return ""

def validate_table(df,pathway_label,cv,binding_col_override=None):
    results,invalid,unknown=[],[],[]
    shape=df.shape

    def _add(cid,desc,st,n,detail=""):
        results.append({"pathway":pathway_label,"check_id":cid,"description":desc,
                        "table_rows":shape[0],"table_cols":shape[1],"status":st,
                        "n_issues":n,"detail":detail})
        print_check(cid,st,n,desc)

    def _flag(idx,cid,reason):
        invalid.append({"pathway":pathway_label,"check_id":cid,"reason":reason,
                        "source_row":int(idx)+2,**_row_dict(df,idx)})

    bc=_detect_binding_col(df,binding_col_override)
    print(f"\n  Pathway: {pathway_label}  ({shape[0]} rows × {shape[1]} cols)")
    print(f"  Binding column: '{bc}'" if bc else "  ⚠  No binding partner column found")

    # C1 — tiered (FIX-2)
    ma=[c for c in REQUIRED_ANNOTATION if c not in df.columns]
    mo=[c for c in REQUIRED_OPTIONAL if c not in df.columns]
    mc = [c for c in REQUIRED_CORE if c not in df.columns]
    ma = [c for c in REQUIRED_ANNOTATION if c not in df.columns]
    mo = [c for c in REQUIRED_OPTIONAL if c not in df.columns]
    _add("C1","Required structural columns","FAIL" if mc else "PASS",len(mc),", ".join(mc))
    _add("C1b","Required annotation columns","WARN" if ma else "PASS",len(ma),", ".join(ma))
    _add("C1c","Optional annotation columns","WARN" if mo else "PASS",len(mo),f"Not yet generated: {', '.join(mo)}")

    # C2
    oid=next((c for c in ("occurrence_id","domain_occurrence_id") if c in df.columns),"occurrence_id")
    if oid in df.columns:
        blank_oid=df[oid].isna()|(df[oid].str.strip()=="")
        n=int(blank_oid.sum()); _add("C2",f"Blank {oid}","FAIL" if n else "PASS",n)
        for idx in df[blank_oid].index: _flag(idx,"C2",f"blank {oid}")
    else:
        blank_oid=pd.Series(False,index=df.index); _add("C2",f"Blank {oid}","WARN",0,"Column absent")

    # C3a/b
    if oid in df.columns:
        dup=df.duplicated(subset=oid,keep=False)&~blank_oid; n=int(dup.sum())
        _add("C3a",f"Duplicate {oid}","FAIL" if n else "PASS",n)
        for idx in df[dup].index: _flag(idx,"C3a",f"dup {oid}: '{df.loc[idx,oid]}'")
    else: _add("C3a",f"Duplicate {oid}","WARN",0,"absent")
    kc=[c for c in df.columns if c!="batch_gene"]
    dr=df.duplicated(subset=kc,keep=False)&~blank_oid
    _add("C3b","Full-row duplicate excl batch_gene","WARN" if dr.sum() else "PASS",int(dr.sum()))

    # C4
    if "protein_accession" in df.columns:
        fil=df["protein_accession"].notna()&(df["protein_accession"].str.strip()!="")
        bad=fil&~df["protein_accession"].str.strip().str.fullmatch(RE_UNIPROT.pattern)
        n=int(bad.sum()); _add("C4","protein_accession UniProtKB format","FAIL" if n else "PASS",n)
        for idx in df[bad].index: _flag(idx,"C4",f"bad: '{df.loc[idx,'protein_accession']}'")
    else: _add("C4","protein_accession format","WARN",0,"absent")

    # C5
    if "domain_type_id" in df.columns:
        fil=df["domain_type_id"].notna()&(df["domain_type_id"].str.strip()!="")
        bad=fil&~df["domain_type_id"].str.strip().str.fullmatch(RE_DT_ID.pattern)
        n=int(bad.sum()); _add("C5","domain_type_id format","FAIL" if n else "PASS",n)
        for idx in df[bad].index: _flag(idx,"C5",f"bad: '{df.loc[idx,'domain_type_id']}'")
    else: _add("C5","domain_type_id format","WARN",0,"absent")

    # C6
    if "start" in df.columns and "end" in df.columns:
        cb=[]
        for idx,row in df[["start","end"]].iterrows():
            parts=[]
            try: si=int(float(row["start"]))
            except: si=None; parts.append(f"start not int: '{row['start']}'") if pd.notna(row["start"]) else None
            try: ei=int(float(row["end"]))
            except: ei=None; parts.append(f"end not int: '{row['end']}'") if pd.notna(row["end"]) else None
            if si and si<=0: parts.append(f"start≤0:{si}")
            if ei and ei<=0: parts.append(f"end≤0:{ei}")
            if si and ei and si>=ei: parts.append(f"start({si})≥end({ei})")
            if parts: cb.append((idx,"; ".join(parts)))
        n=len(cb); _add("C6","start/end valid (positive, start<end)","FAIL" if n else "PASS",n)
        for idx,r in cb: _flag(idx,"C6",r)
    else: _add("C6","start/end coordinates","WARN",0,"absent")

    # C7
    go_bad=[(idx,f"{col}='{tok}'") for col in GO_COLS if col in df.columns
            for idx,val in df[col].items() for tok in split_cell(val) if not RE_GO.fullmatch(tok)]
    n=len(go_bad); _add("C7","GO term format","FAIL" if n else "PASS",n)
    for idx,r in go_bad: _flag(idx,"C7",f"bad GO: {r}")

    # C8 — auto-fix double prefix (FIX-3)
    total_fixes=0; chebi_residual=[]
    if "binding_target_chebi" in df.columns:
        fixed_vals=[]
        for idx,val in df["binding_target_chebi"].items():
            if pd.isna(val): fixed_vals.append(val); continue
            fv,nf=_normalise_chebi(str(val)); total_fixes+=nf; fixed_vals.append(fv)
            for tok in split_cell(fv):
                if not RE_CHEBI_OK.fullmatch(tok): chebi_residual.append((idx,f"bad after fix: '{tok}'"))
        df=df.copy(); df["binding_target_chebi"]=fixed_vals
    if total_fixes>0 and not chebi_residual:
        _add("C8","binding_target_chebi format (CHEBI:nnnnn)","WARN",total_fixes,f"Auto-fixed {total_fixes} CHEBI:CHEBI: double-prefix values in-memory")
    elif chebi_residual:
        n=len(chebi_residual); _add("C8","binding_target_chebi format","FAIL",n,f"Also auto-fixed {total_fixes} double-prefix")
        for idx,r in chebi_residual: _flag(idx,"C8",r)
    else: _add("C8","binding_target_chebi format","PASS",0)

    # C9
    if "protein_reactome_ids" in df.columns:
        rb=[(idx,tok) for idx,val in df["protein_reactome_ids"].items() for tok in split_cell(val) if not RE_REACTOME.fullmatch(tok)]
        n=len(rb); _add("C9","Reactome ID format","FAIL" if n else "PASS",n)
        for idx,tok in rb: _flag(idx,"C9",f"bad: '{tok}'")
    else: _add("C9","Reactome ID format","WARN",0,"absent")

    # C10
    if "protein_kegg_ids" in df.columns:
        kb=[(idx,tok) for idx,val in df["protein_kegg_ids"].items() for tok in split_cell(val) if not RE_KEGG.fullmatch(tok)]
        n=len(kb); _add("C10","KEGG ID format","FAIL" if n else "PASS",n)
        for idx,tok in kb: _flag(idx,"C10",f"bad: '{tok}'")
    else: _add("C10","KEGG ID format","WARN",0,"absent")

    # C11-C16
    _bc16=bc if bc else "binding_partner"
    cv_checks=[("C11","positional_category","positional_or_topology_term"),
               ("C12","copy_number_property","CopyNumberCategory"),
               ("C13","topology_context","positional_or_topology_term"),
               ("C14","proximity_categories","positional_or_topology_term"),
               ("C15","function_tags","functional_term"),
               ("C16",_bc16,"binding_target_term")]
    all_valid=frozenset(tok for s in cv.values() for tok in s)
    for cid,col,ttype in cv_checks:
        if col not in df.columns:
            _add(cid,f"{col} CV-backed ({ttype})","WARN",0,f"Column '{col}' absent"); continue
        valid=cv.get(ttype,frozenset()); br=defaultdict(list)
        for idx,val in df[col].items():
            for tok in split_cell(val):
                if not is_cv_backed(tok,valid):
                    br[idx].append(tok)
                    if not is_cv_backed(tok,all_valid):
                        unknown.append({"pathway":pathway_label,"column":col,"token":tok,"term_type":ttype,"source_row":int(idx)+2})
        nt=sum(len(v) for v in br.values())
        _add(cid,f"{col} CV-backed ({ttype})","FAIL" if nt else "PASS",nt,f"{len(br)} rows")
        for idx,toks in br.items(): _flag(idx,cid,f"unmatched in '{col}': "+"; ".join(toks))

    # C17
    nsv=[(idx,col,tok,ns,rn,rd) for col,rn,rd in NAMESPACE_RULES if col in df.columns
         for idx,val in df[col].items() for tok in split_cell(val)
         if (ns:=get_id_namespace(tok)) and ns!=rn]
    n=len(nsv); _add("C17","Relation-to-namespace rules","FAIL" if n else "PASS",n)
    for idx,col,tok,ns,rn,rd in nsv: _flag(idx,"C17",f"{rd}: '{tok}' ns={ns} in '{col}'")

    return results,invalid,unknown


# ── Missing term report ───────────────────────────────────────────────────────
def build_missing_term_report(dfs,cv):
    counter=defaultdict(int)
    for pl,df in dfs.items():
        for col,ttype in CV_COL_TO_TERM_TYPE.items():
            if col not in df.columns: continue
            valid=cv.get(ttype,frozenset())
            for val in df[col].dropna():
                for tok in split_cell(val):
                    if not is_cv_backed(tok,valid): counter[(pl,col,tok)]+=1
    av=all_valid_tokens(cv); rows=[]
    for (pw,col,tok),cnt in sorted(counter.items(),key=lambda x:(-x[1],x[0])):
        m=get_close_matches(tok,av,n=1,cutoff=0.6)
        cl=m[0] if m else ""
        rows.append({"pathway":pw,"column":col,"unmatched_token":tok,"occurrence_count":cnt,
                     "closest_match":cl,"suggested_action":f"Closest: '{cl}' — check spelling" if cl else "Add to annotation_table.xlsx"})
    return pd.DataFrame(rows)


# ── Pathway metrics ───────────────────────────────────────────────────────────
def compute_pathway_metrics(df,label,cv):
    m={"pathway":label}
    m["n_rows"]=len(df)
    m["n_unique_proteins"]=df["protein_accession"].nunique() if "protein_accession" in df.columns else 0
    m["n_unique_domain_types"]=df["domain_type_id"].nunique() if "domain_type_id" in df.columns else 0
    for col in CV_COL_TO_TERM_TYPE:
        if col in df.columns:
            fil=df[col].notna()&(df[col].str.strip()!="")
            m[f"pct_filled_{col}"]=round(100*fil.sum()/max(len(df),1),1)
            toks=set()
            for val in df[col].dropna():
                for tok in split_cell(val): toks.add(extract_label(tok).lower())
            m[f"n_unique_terms_{col}"]=len(toks)
            valid=cv.get(CV_COL_TO_TERM_TYPE[col],frozenset())
            tot=bk=0
            for val in df[col].dropna():
                for tok in split_cell(val):
                    tot+=1
                    if is_cv_backed(tok,valid): bk+=1
            m[f"cv_backed_rate_{col}"]=round(100*bk/max(tot,1),1)
        else:
            m[f"pct_filled_{col}"]=m[f"n_unique_terms_{col}"]=m[f"cv_backed_rate_{col}"]=None
    m["schema_coverage_pct"]=round(100*sum(1 for c in REQUIRED_COLUMNS if c in df.columns)/len(REQUIRED_COLUMNS),1)
    return m


# ── Transfer metrics ──────────────────────────────────────────────────────────
def compute_transfer_metrics(df_ref,df_tgt,cv,lr="cAMP",lt="other"):
    m={"reference_pathway":lr,"target_pathway":lt}
    rc,tc=set(df_ref.columns),set(df_tgt.columns)
    sh=rc&tc
    m["schema_reuse_rate"]=round(len(sh)/max(len(rc),1),4)
    m["n_shared_columns"]=len(sh); m["n_ref_only_columns"]=len(rc-tc); m["n_target_only_columns"]=len(tc-rc)
    def _ct(df):
        t=set()
        for col in CV_COL_TO_TERM_TYPE:
            if col not in df.columns: continue
            for val in df[col].dropna():
                for tok in split_cell(val): t.add(extract_label(tok).lower())
        return t
    rt,tt=_ct(df_ref),_ct(df_tgt); st=rt&tt
    m["vocab_reuse_rate"]=round(len(st)/max(len(rt),1),4)
    m["n_ref_cv_tokens"]=len(rt); m["n_target_cv_tokens"]=len(tt); m["n_shared_cv_tokens"]=len(st); m["n_target_only_tokens"]=len(tt-rt)
    test=sum(1 for col,_,_ in NAMESPACE_RULES if col in df_tgt.columns)
    m["rule_portability_rate"]=round(test/max(len(NAMESPACE_RULES),1),4); m["n_testable_rules"]=test; m["n_total_rules"]=len(NAMESPACE_RULES)
    av=frozenset(tok for s in cv.values() for tok in s)
    tot=cov=0
    for col,ttype in CV_COL_TO_TERM_TYPE.items():
        if col not in df_tgt.columns: continue
        valid=cv.get(ttype,frozenset())
        for val in df_tgt[col].dropna():
            for tok in split_cell(val):
                tot+=1
                if is_cv_backed(tok,valid): cov+=1
    m["mapping_completeness"]=round(cov/max(tot,1),4); m["n_target_tokens_total"]=tot; m["n_target_tokens_mapped"]=cov
    reg=set(tok.lower() for s in cv.values() for tok in s if not tok.islower())
    m["n_reused_terms"]=len(rt&tt&reg); m["n_new_terms_required"]=len(tt-reg); m["n_registry_terms_unused"]=len(reg-rt-tt)
    return m


# ── Jaccard ───────────────────────────────────────────────────────────────────
def jaccard(a,b): u=a|b; return round(len(a&b)/max(len(u),1),4)

def compute_overlap_metrics(df1,df2,l1="cAMP",l2="other"):
    def _toks(df,col):
        t=set()
        if col not in df.columns: return t
        for val in df[col].dropna():
            for tok in split_cell(val): t.add(extract_label(tok).lower())
        return t
    def _cv(df,col):
        if col not in df.columns: return set()
        return set(df[col].dropna().str.strip().str.lower().unique())-{"","nan"}
    cats=[("domain_types",_cv(df1,"domain_type_id"),_cv(df2,"domain_type_id")),
          ("proteins",_cv(df1,"protein_accession"),_cv(df2,"protein_accession")),
          ("function_tags",_toks(df1,"function_tags"),_toks(df2,"function_tags")),
          ("binding_targets",_toks(df1,_detect_binding_col(df1)),_toks(df2,_detect_binding_col(df2))),
          ("topology_context",_toks(df1,"topology_context"),_toks(df2,"topology_context")),
          ("positional_category",_toks(df1,"positional_category"),_toks(df2,"positional_category")),
          ("proximity_categories",_toks(df1,"proximity_categories"),_toks(df2,"proximity_categories")),
          ("copy_number_property",_toks(df1,"copy_number_property"),_toks(df2,"copy_number_property")),
          ("go_terms_mf",_toks(df1,"go_terms_mf"),_toks(df2,"go_terms_mf")),
          ("go_terms_bp",_toks(df1,"go_terms_bp"),_toks(df2,"go_terms_bp"))]
    rows=[]
    for cat,s1,s2 in cats:
        rows.append({"category":cat,f"{l1}_size":len(s1),f"{l2}_size":len(s2),
                     "shared":len(s1&s2),f"unique_to_{l1}":len(s1-s2),f"unique_to_{l2}":len(s2-s1),
                     "jaccard":jaccard(s1,s2),f"{l1}_items":"; ".join(sorted(s1)),
                     f"{l2}_items":"; ".join(sorted(s2)),"shared_items":"; ".join(sorted(s1&s2)),
                     f"unique_to_{l1}_items":"; ".join(sorted(s1-s2)),f"unique_to_{l2}_items":"; ".join(sorted(s2-s1))})
    return pd.DataFrame(rows)


# ── Comparison report ─────────────────────────────────────────────────────────
def build_comparison_report(df1,df2,l1,l2,cv):
    def _cv_rate(df,col,ttype):
        if col not in df.columns: return 0.0
        valid=cv.get(ttype,frozenset()); tot=bk=0
        for val in df[col].dropna():
            for tok in split_cell(val):
                tot+=1
                if is_cv_backed(tok,valid): bk+=1
        return round(100*bk/max(tot,1),1)
    def _ut(df,col):
        if col not in df.columns: return 0
        t=set()
        for val in df[col].dropna():
            for tok in split_cell(val): t.add(extract_label(tok).lower())
        return len(t)
    bc1=_detect_binding_col(df1); bc2=_detect_binding_col(df2)
    metrics=[("Total occurrence rows",len(df1),len(df2)),
             ("Unique proteins",df1["protein_accession"].nunique() if "protein_accession" in df1.columns else 0,df2["protein_accession"].nunique() if "protein_accession" in df2.columns else 0),
             ("Unique domain types",df1["domain_type_id"].nunique() if "domain_type_id" in df1.columns else 0,df2["domain_type_id"].nunique() if "domain_type_id" in df2.columns else 0),
             ("Schema columns present",len(df1.columns),len(df2.columns)),
             ("Unique function tags",_ut(df1,"function_tags"),_ut(df2,"function_tags")),
             ("Unique topology terms",_ut(df1,"topology_context"),_ut(df2,"topology_context")),
             ("Unique positional terms",_ut(df1,"positional_category"),_ut(df2,"positional_category")),
             ("Unique proximity terms",_ut(df1,"proximity_categories"),_ut(df2,"proximity_categories")),
             ("Unique binding targets",_ut(df1,bc1),_ut(df2,bc2)),
             ("Unique copy-number terms",_ut(df1,"copy_number_property"),_ut(df2,"copy_number_property")),
             ("CV-backed rate: function_tags (%)",_cv_rate(df1,"function_tags","functional_term"),_cv_rate(df2,"function_tags","functional_term")),
             ("CV-backed rate: topology_context (%)",_cv_rate(df1,"topology_context","positional_or_topology_term"),_cv_rate(df2,"topology_context","positional_or_topology_term")),
             ("CV-backed rate: positional_category (%)",_cv_rate(df1,"positional_category","positional_or_topology_term"),_cv_rate(df2,"positional_category","positional_or_topology_term")),
             ("CV-backed rate: binding_partner (%)",_cv_rate(df1,bc1,"binding_target_term"),_cv_rate(df2,bc2,"binding_target_term")),
             ("CV-backed rate: copy_number_property (%)",_cv_rate(df1,"copy_number_property","CopyNumberCategory"),_cv_rate(df2,"copy_number_property","CopyNumberCategory"))]
    rows=[]
    for name,v1,v2 in metrics:
        try: diff=round(float(v2)-float(v1),2)
        except: diff=None
        rows.append({"metric":name,l1:v1,l2:v2,"difference":diff})
    return pd.DataFrame(rows)


# ════════════════════════════════════════════════════════════════════════════
# Figures
# ════════════════════════════════════════════════════════════════════════════

def plot_validation_summary(r1,r2,out,l1,l2):
    sp={"PASS":PAL["pass"],"WARN":PAL["warn"],"FAIL":PAL["fail"]}
    fig=plt.figure(figsize=(12,5.5))
    gs=gridspec.GridSpec(1,3,figure=fig,wspace=0.38,left=0.06,right=0.97,top=0.88,bottom=0.08,width_ratios=[2,3,3])
    ax_b,ax_g,ax_t=fig.add_subplot(gs[0]),fig.add_subplot(gs[1]),fig.add_subplot(gs[2])

    for yi,(lbl,res) in enumerate([(l2,r2),(l1,r1)]):
        cnt=Counter(r["status"] for r in res); tot=sum(cnt.values()); left=0
        for st in ["PASS","WARN","FAIL"]:
            v=cnt.get(st,0)
            if v==0: continue
            ax_b.barh(yi,v,height=0.5,left=left,color=sp[st],edgecolor="white",lw=0.7)
            col="white" if st!="WARN" else "#2c3e50"
            ax_b.text(left+v/2,yi,f"{v}\n({100*v/max(tot,1):.0f}%)",ha="center",va="center",fontsize=6.5,fontweight="bold",color=col)
            left+=v
    ax_b.set_yticks([0,1]); ax_b.set_yticklabels([l2,l1],fontsize=8)
    ax_b.set_xlabel("Number of checks"); ax_b.set_title("Validation Summary",fontweight="bold",pad=6)
    ax_b.legend(handles=[mpatches.Patch(facecolor=sp[s],label=s) for s in ["PASS","WARN","FAIL"]],fontsize=6.5,loc="lower right",handlelength=1.2)
    _panel_label(ax_b,"A")

    all_cids=sorted({r["check_id"] for r in r1+r2})
    snum={"PASS":2,"WARN":1,"FAIL":0}
    gd=np.full((2,len(all_cids)),np.nan)
    for pi,res in enumerate([r1,r2]):
        lk={r["check_id"]:r["status"] for r in res}
        for ci,cid in enumerate(all_cids):
            if cid in lk: gd[pi,ci]=snum.get(lk[cid],-1)
    cmap3=LinearSegmentedColormap.from_list("pf",[PAL["fail"],PAL["warn"],PAL["pass"]],N=3)
    ax_g.imshow(gd,aspect="auto",cmap=cmap3,vmin=-0.5,vmax=2.5,interpolation="nearest")
    ax_g.set_xticks(range(len(all_cids))); ax_g.set_xticklabels(all_cids,rotation=60,ha="right",fontsize=6.5)
    ax_g.set_yticks([0,1]); ax_g.set_yticklabels([l1,l2],fontsize=8)
    ax_g.set_title("Per-Check Status Grid",fontweight="bold",pad=6)
    for pi,res in enumerate([r1,r2]):
        lk={r["check_id"]:r["status"] for r in res}
        for ci,cid in enumerate(all_cids):
            st=lk.get(cid,"")
            ax_g.text(ci,pi,st[:1],ha="center",va="center",fontsize=5.5,
                      color="white" if st=="FAIL" else "#2c3e50",fontweight="bold")
    _panel_label(ax_g,"B")

    ax_t.axis("off")
    rows_t=[[r["check_id"],r["status"],lbl[:10],r["description"][:32]+"…" if len(r["description"])>32 else r["description"]]
             for lbl,res in [(l1,r1),(l2,r2)] for r in res if r["status"] in ("FAIL","WARN")]
    if rows_t:
        tbl=ax_t.table(cellText=rows_t,colLabels=["Check","Status","Pathway","Description"],loc="center",cellLoc="left",bbox=[0,0,1,1])
        tbl.auto_set_font_size(False); tbl.set_fontsize(6.5); tbl.auto_set_column_width([0,1,2,3])
        for (r,c),cell in tbl.get_celld().items():
            cell.set_edgecolor("#dddddd"); cell.set_linewidth(0.4)
            if r==0: cell.set_facecolor("#2c3e50"); cell.set_text_props(color="white",fontweight="bold")
            elif rows_t[r-1][1]=="FAIL": cell.set_facecolor("#ffeaea")
            elif rows_t[r-1][1]=="WARN": cell.set_facecolor("#fff8e1")
            else: cell.set_facecolor("white")
    else: ax_t.text(0.5,0.5,"All checks PASS",ha="center",va="center",color=PAL["pass"],fontsize=10,fontweight="bold")
    ax_t.set_title("FAIL / WARN Detail",fontweight="bold",pad=6); _panel_label(ax_t,"C")
    fig.suptitle("Ontology Validation Report",fontsize=10,fontweight="bold",y=0.97)
    fig.savefig(out); plt.close(fig); print(f"  ✓  Plot saved → {out}")


def plot_cv_coverage_heatmap(m1,m2,l1,l2,out):
    cv_c=["positional_category","topology_context","proximity_categories","function_tags","binding_partner","copy_number_property"]
    present=[c for c in cv_c if m1.get(f"pct_filled_{c}") is not None or m2.get(f"pct_filled_{c}") is not None]
    if not present: print("  ⚠  No CV data — skipping"); return
    nr=len(present); data=np.zeros((nr,4))
    for i,col in enumerate(present):
        data[i,0]=float(m1.get(f"pct_filled_{col}") or 0)
        data[i,1]=float(m1.get(f"cv_backed_rate_{col}") or 0)
        data[i,2]=float(m2.get(f"pct_filled_{col}") or 0)
        data[i,3]=float(m2.get(f"cv_backed_rate_{col}") or 0)
    fig,(ax,ax_l)=plt.subplots(1,2,figsize=(7.5,max(2.8,nr*0.6+1.2)),gridspec_kw={"width_ratios":[8,1]})
    im=ax.imshow(data,aspect="auto",cmap="YlGn",vmin=0,vmax=100,interpolation="nearest")
    ax.set_xticks(range(4)); ax.set_xticklabels([f"{l1[:10]}\nFill %",f"{l1[:10]}\nCV-backed %",f"{l2[:10]}\nFill %",f"{l2[:10]}\nCV-backed %"],fontsize=7,fontweight="bold")
    ax.set_yticks(range(nr)); ax.set_yticklabels(present,fontsize=7.5); ax.axvline(1.5,color="white",lw=3)
    for i in range(nr):
        for j in range(4):
            v=data[i,j]; fc="white" if v>60 else "#2c3e50"
            ax.text(j,i,f"{v:.0f}%",ha="center",va="center",fontsize=7,fontweight="bold",color=fc)
    cb=fig.colorbar(im,cax=ax_l); cb.ax.tick_params(labelsize=6); cb.set_label("Percentage (%)",fontsize=7)
    ax.set_title(f"CV Column Coverage  ·  {l1}  |  {l2}\nFill rate = rows with value   CV-backed = tokens in registry",fontsize=8,fontweight="bold",pad=10)
    fig.tight_layout(); fig.savefig(out); plt.close(fig); print(f"  ✓  Plot saved → {out}")


def plot_term_reuse(transfer,out,l1,l2):
    nr=transfer.get("n_ref_cv_tokens",0); nt=transfer.get("n_target_cv_tokens",0)
    ns=transfer.get("n_shared_cv_tokens",0); nu=transfer.get("n_registry_terms_unused",0)
    nn=transfer.get("n_new_terms_required",0)
    fig,axes=plt.subplots(1,3,figsize=(11,3.2),gridspec_kw={"width_ratios":[3,1.8,1.5]})
    ax,ax2,ax3=axes
    bh=0.42
    for yi,(lbl,ns_,nu_,cu) in enumerate([(l1,ns,nr-ns,PAL["ref_only"]),(l2,ns,nt-ns,PAL["tgt_only"])]):
        ax.barh(yi,ns_,height=bh,color=PAL["shared"],edgecolor="white",lw=0.6,label="Shared" if yi==0 else "_")
        ax.barh(yi,nu_,height=bh,left=ns_,color=cu,edgecolor="white",lw=0.6,label=f"Unique to {lbl}")
        for w,left in [(ns_,0),(nu_,ns_)]:
            if w>0: ax.text(left+w/2,yi,str(w),ha="center",va="center",fontsize=7,fontweight="bold",color="white")
    ax.set_yticks([0,1]); ax.set_yticklabels([l1,l2]); ax.set_xlabel("Unique CV tokens")
    ax.set_title("Vocabulary Usage per Pathway",fontweight="bold",pad=6); ax.spines["left"].set_visible(False); ax.tick_params(axis="y",length=0)
    hl,ll=ax.get_legend_handles_labels(); ax.legend(dict(zip(ll,hl)).values(),dict(zip(ll,hl)).keys(),fontsize=6.5,handlelength=1.2,loc="lower right")
    _panel_label(ax,"A")
    cats=["Reused\n(both)",f"{l1[:8]}\nonly","New terms\nneeded","Registry\nunused"]
    vals=[ns,nr-ns,nn,nu]; cols=[PAL["pass"],PAL["ref_only"],PAL["warn"],PAL["grey"]]
    yp=range(len(cats)); bars=ax2.barh(list(yp),vals,height=0.55,color=cols,edgecolor="white",lw=0.6)
    for bar,v in zip(bars,vals):
        if v>0: ax2.text(v+max(vals)*0.02,bar.get_y()+bar.get_height()/2,str(v),va="center",fontsize=7,color="#2c3e50")
    ax2.set_yticks(list(yp)); ax2.set_yticklabels(cats,fontsize=7); ax2.set_xlabel("Count"); ax2.set_title("Term Transfer Status",fontweight="bold",pad=6); _panel_label(ax2,"B")
    used=nr-ns+nt-ns+ns; sizes=[used,max(nu,0)]
    wedges,_=ax3.pie(sizes,colors=[PAL["green"],PAL["grey"]],startangle=90,wedgeprops=dict(width=0.52,edgecolor="white",lw=1.5))
    pct=100*used/max(used+max(nu,0),1)
    ax3.text(0,0,f"{pct:.0f}%\nused",ha="center",va="center",fontsize=9,fontweight="bold",color="#2c3e50")
    ax3.set_title("Registry\nCoverage",fontweight="bold",pad=6)
    ax3.legend(handles=[mpatches.Patch(color=PAL["green"],label=f"Used ({used})"),mpatches.Patch(color=PAL["grey"],label=f"Unused ({max(nu,0)})")],fontsize=6.5,loc="lower center",bbox_to_anchor=(0.5,-0.18))
    _panel_label(ax3,"C",x=-0.06)
    fig.suptitle(f"Ontology Vocabulary Transfer: {l1} → {l2}",fontsize=9,fontweight="bold",y=1.02)
    fig.tight_layout(); fig.savefig(out); plt.close(fig); print(f"  ✓  Plot saved → {out}")


def plot_jaccard_heatmap(ov,out,l1,l2):
    cats=ov["category"].tolist(); jacs=ov["jaccard"].astype(float).tolist(); sh=ov["shared"].astype(int).tolist()
    u1c=f"unique_to_{l1}"; u2c=f"unique_to_{l2}"
    u1=ov[u1c].astype(int).tolist() if u1c in ov.columns else (ov[f"{l1}_size"].astype(int)-ov["shared"].astype(int)).tolist()
    u2=ov[u2c].astype(int).tolist() if u2c in ov.columns else (ov[f"{l2}_size"].astype(int)-ov["shared"].astype(int)).tolist()
    n=len(cats); y=np.arange(n)
    fig,(ax_l,ax_c)=plt.subplots(1,2,figsize=(9,max(3.5,n*0.52+1.2)),gridspec_kw={"width_ratios":[2.2,3]})
    cmap=plt.cm.RdYlGn
    for i,(j,cat) in enumerate(zip(jacs,cats)):
        ax_l.plot([0,j],[i,i],color="#d5d5d5",lw=1.2,zorder=1)
        ax_l.scatter(j,i,s=65,color=cmap(j),zorder=3,linewidths=0.6,edgecolors="white")
        ax_l.text(j+0.025,i,f"{j:.3f}",va="center",fontsize=6.5,color="#2c3e50")
    ax_l.axvline(0.5,color="#aaaaaa",lw=0.8,ls="--",alpha=0.7)
    ax_l.set_xlim(-0.04,1.2); ax_l.set_yticks(y); ax_l.set_yticklabels(cats,fontsize=7.5)
    ax_l.set_xlabel("Jaccard similarity index"); ax_l.set_title(f"Biological Overlap\n{l1}  vs  {l2}",fontweight="bold",pad=6)
    sm=plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,1)); sm.set_array([])
    cb=fig.colorbar(sm,ax=ax_l,shrink=0.5,pad=0.03,aspect=16); cb.set_label("Jaccard",fontsize=6.5); cb.ax.tick_params(labelsize=6)
    _panel_label(ax_l,"A")
    bh=0.22
    ax_c.barh(y+bh,sh,height=bh,color=PAL["shared"],label="Shared")
    ax_c.barh(y,u1,height=bh,color=PAL["ref_only"],label=f"Unique {l1}")
    ax_c.barh(y-bh,u2,height=bh,color=PAL["tgt_only"],label=f"Unique {l2}")
    for i,(sv,u1v,u2v) in enumerate(zip(sh,u1,u2)):
        for val,dy in [(sv,bh),(u1v,0),(u2v,-bh)]:
            if val>0: ax_c.text(val+0.3,i+dy,str(val),va="center",fontsize=6,color="#2c3e50")
    ax_c.set_yticks(y); ax_c.set_yticklabels([])
    ax_c.set_xlabel("Item count"); ax_c.set_title("Set Breakdown\n(shared / unique per pathway)",fontweight="bold",pad=6)
    ax_c.legend(fontsize=6.5,loc="lower right",handlelength=1.2); _panel_label(ax_c,"B")
    fig.tight_layout(); fig.savefig(out); plt.close(fig); print(f"  ✓  Plot saved → {out}")


def plot_transfer_radar(transfer,out,l1,l2):
    mls=["Schema\nreuse","Vocabulary\nreuse","Rule\nportability","Mapping\ncompleteness"]
    vals=[transfer.get(k,0.0) for k in ["schema_reuse_rate","vocab_reuse_rate","rule_portability_rate","mapping_completeness"]]
    N=len(mls); angles=np.linspace(0,2*np.pi,N,endpoint=False).tolist(); angles+=angles[:1]; vs=vals+vals[:1]
    fig,ax=plt.subplots(figsize=(4.5,4.5),subplot_kw={"polar":True})
    for lv,ls in [(0.25,":"),(0.5,"--"),(0.75,":"),(1.0,"-")]:
        ax.plot(angles,[lv]*(N+1),color="#cccccc",lw=0.6,ls=ls,zorder=1)
        ax.text(np.pi/2,lv+0.04,f"{lv:.0%}",ha="center",va="bottom",fontsize=5.5,color="#999999")
    ax.fill(angles,vs,color=PAL["blue"],alpha=0.18,zorder=2)
    ax.plot(angles,vs,color=PAL["blue"],lw=2.0,marker="o",ms=6,markerfacecolor=PAL["blue"],markeredgecolor="white",markeredgewidth=1.0,zorder=3)
    for angle,val in zip(angles[:-1],vals):
        ax.text(angle,val+0.12,f"{val:.1%}",ha="center",va="center",fontsize=7.5,fontweight="bold",color=PAL["blue"])
    ax.set_xticks(angles[:-1]); ax.set_xticklabels(mls,fontsize=8,color="#2c3e50"); ax.set_ylim(0,1.18); ax.set_yticks([])
    ax.spines["polar"].set_linewidth(0.7); ax.spines["polar"].set_color("#cccccc")
    ax.set_title(f"Ontology Transfer Metrics\n{l1} → {l2}",fontsize=9,fontweight="bold",pad=22)
    nr=transfer.get("n_new_terms_required",0); nrs=transfer.get("n_reused_terms",0)
    fig.text(0.5,0.02,f"Reused terms: {nrs}   ·   New terms needed: {nr}",ha="center",fontsize=7,color="#555555",bbox=dict(boxstyle="round,pad=0.3",fc="#f5f5f5",ec="#cccccc",lw=0.5))
    fig.tight_layout(rect=[0,0.06,1,1]); fig.savefig(out); plt.close(fig); print(f"  ✓  Plot saved → {out}")


def plot_pathway_comparison(comp,ov,transfer,l1,l2,out):
    fig=plt.figure(figsize=(11,11))
    gs=fig.add_gridspec(3,2,hspace=0.55,wspace=0.40,left=0.09,right=0.97,top=0.93,bottom=0.04)
    axes=[fig.add_subplot(gs[r,c]) for r in range(3) for c in range(2)]
    ax_A,ax_B,ax_C,ax_D,ax_E,ax_F=axes
    c1,c2=PAL["blue"],PAL["orange"]
    def _cv(name,lbl):
        r=comp[comp["metric"]==name]; return float(r[lbl].values[0]) if not r.empty else 0.0
    # A
    lA=["Occurrences","Unique\nproteins","Domain\ntypes"]; nA=["Total occurrence rows","Unique proteins","Unique domain types"]
    x=np.arange(3); w=0.33; v1=[_cv(n,l1) for n in nA]; v2=[_cv(n,l2) for n in nA]
    ax_A.bar(x-w/2,v1,w,color=c1,edgecolor="white",lw=0.5,label=l1); ax_A.bar(x+w/2,v2,w,color=c2,edgecolor="white",lw=0.5,label=l2)
    m=max(v1+v2)
    for xi,(a,b) in enumerate(zip(v1,v2)):
        ax_A.text(xi-w/2,a+m*0.02,f"{int(a)}",ha="center",va="bottom",fontsize=6.5)
        ax_A.text(xi+w/2,b+m*0.02,f"{int(b)}",ha="center",va="bottom",fontsize=6.5)
    ax_A.set_xticks(x); ax_A.set_xticklabels(lA); ax_A.set_ylabel("Count"); ax_A.set_title("Dataset Size",fontweight="bold",pad=6)
    ax_A.legend(fontsize=6.5,handlelength=1.2,loc="upper right"); _panel_label(ax_A,"A")
    # B
    lB=["Function\ntags","Topology","Positional","Proximity","Binding\ntargets"]
    nB=["Unique function tags","Unique topology terms","Unique positional terms","Unique proximity terms","Unique binding targets"]
    xB=np.arange(5); wB=0.35; v1B=[_cv(n,l1) for n in nB]; v2B=[_cv(n,l2) for n in nB]
    ax_B.bar(xB-wB/2,v1B,wB,color=c1,edgecolor="white",lw=0.5); ax_B.bar(xB+wB/2,v2B,wB,color=c2,edgecolor="white",lw=0.5)
    mB=max(v1B+v2B+[1])
    for xi,(a,b) in enumerate(zip(v1B,v2B)):
        ax_B.text(xi-wB/2,a+mB*0.02,f"{int(a)}",ha="center",va="bottom",fontsize=6); ax_B.text(xi+wB/2,b+mB*0.02,f"{int(b)}",ha="center",va="bottom",fontsize=6)
    ax_B.set_xticks(xB); ax_B.set_xticklabels(lB); ax_B.set_ylabel("Unique terms"); ax_B.set_title("Annotation Vocabulary Richness",fontweight="bold",pad=6); _panel_label(ax_B,"B")
    # C
    rc=[("CV-backed rate: function_tags (%)","Function tags"),("CV-backed rate: topology_context (%)","Topology"),("CV-backed rate: positional_category (%)","Positional"),("CV-backed rate: binding_partner (%)","Binding partner")]
    yC=np.arange(len(rc)); r1C=[_cv(m,l1) for m,_ in rc]; r2C=[_cv(m,l2) for m,_ in rc]
    for i,(r1,r2) in enumerate(zip(r1C,r2C)):
        ax_C.plot([r1],[i-0.15],"o",color=c1,ms=7,markeredgecolor="white",markeredgewidth=0.6,zorder=3)
        ax_C.plot([r2],[i+0.15],"D",color=c2,ms=7,markeredgecolor="white",markeredgewidth=0.6,zorder=3)
        ax_C.plot([r1,r2],[i-0.15,i+0.15],color="#cccccc",lw=0.8,zorder=1)
        ax_C.text(r1+1,i-0.15,f"{r1:.0f}%",va="center",fontsize=6,color=c1); ax_C.text(r2+1,i+0.15,f"{r2:.0f}%",va="center",fontsize=6,color=c2)
    ax_C.set_yticks(yC); ax_C.set_yticklabels([s for _,s in rc]); ax_C.set_xlim(-5,115); ax_C.axvline(100,color="#cccccc",lw=0.7,ls="--")
    ax_C.set_xlabel("CV-backed rate (%)"); ax_C.set_title("Ontology Coverage per Column",fontweight="bold",pad=6)
    ax_C.legend(handles=[mpatches.Patch(color=c1,label=l1),mpatches.Patch(color=c2,label=l2)],fontsize=6.5,loc="lower right",handlelength=1.2); _panel_label(ax_C,"C")
    # D
    cats_D=ov["category"].tolist(); jacs_D=ov["jaccard"].astype(float).tolist(); yD=np.arange(len(cats_D))
    cmap_D=plt.cm.RdYlGn
    bars=ax_D.barh(yD,jacs_D,height=0.6,color=[cmap_D(j) for j in jacs_D],edgecolor="white",lw=0.5)
    for bar,j in zip(bars,jacs_D):
        fc="white" if j<0.3 or j>0.7 else "#2c3e50"
        ax_D.text(min(j-0.01,0.98),bar.get_y()+bar.get_height()/2,f"{j:.3f}",ha="right",va="center",fontsize=6.5,fontweight="bold",color=fc,clip_on=True)
    ax_D.axvline(0.5,color="#aaaaaa",lw=0.7,ls="--",alpha=0.7); ax_D.set_xlim(0,1.1)
    ax_D.set_yticks(yD); ax_D.set_yticklabels(cats_D); ax_D.set_xlabel("Jaccard similarity index"); ax_D.set_title("Set Overlap by Category",fontweight="bold",pad=6)
    sm_D=plt.cm.ScalarMappable(cmap=cmap_D,norm=plt.Normalize(0,1)); sm_D.set_array([])
    cb_D=fig.colorbar(sm_D,ax=ax_D,shrink=0.65,pad=0.02,aspect=18); cb_D.ax.tick_params(labelsize=5.5); _panel_label(ax_D,"D")
    # E
    tm=["Schema reuse","Vocabulary reuse","Rule portability","Mapping completeness"]
    tv=[transfer.get(k,0.0) for k in ["schema_reuse_rate","vocab_reuse_rate","rule_portability_rate","mapping_completeness"]]
    yE=np.arange(4); cmap_E=plt.cm.RdYlGn
    bE=ax_E.barh(yE,tv,height=0.55,color=[cmap_E(v) for v in tv],edgecolor="white",lw=0.5)
    for bar,v in zip(bE,tv):
        fc="white" if v<0.35 or v>0.75 else "#2c3e50"
        ax_E.text(min(v-0.01,1.08),bar.get_y()+bar.get_height()/2,f"{v:.1%}",ha="right",va="center",fontsize=7,fontweight="bold",color=fc,clip_on=True)
    ax_E.axvline(0.5,color="#aaaaaa",lw=0.7,ls="--",alpha=0.7); ax_E.set_xlim(0,1.12)
    ax_E.set_yticks(yE); ax_E.set_yticklabels(tm); ax_E.set_xlabel("Score (0→1)"); ax_E.set_title(f"Transfer Performance\n{l1} → {l2}",fontweight="bold",pad=6)
    sm_E=plt.cm.ScalarMappable(cmap=cmap_E,norm=plt.Normalize(0,1)); sm_E.set_array([])
    fig.colorbar(sm_E,ax=ax_E,shrink=0.65,pad=0.02,aspect=18).ax.tick_params(labelsize=5.5); _panel_label(ax_E,"E")
    # F
    ax_F.axis("off")
    rows_F=[["Metric",l1[:12],l2[:12]],
            ["Occurrence rows",str(int(_cv("Total occurrence rows",l1))),str(int(_cv("Total occurrence rows",l2)))],
            ["Unique proteins",str(int(_cv("Unique proteins",l1))),str(int(_cv("Unique proteins",l2)))],
            ["Domain types",str(int(_cv("Unique domain types",l1))),str(int(_cv("Unique domain types",l2)))],
            ["Function tags",str(int(_cv("Unique function tags",l1))),str(int(_cv("Unique function tags",l2)))],
            ["Binding targets",str(int(_cv("Unique binding targets",l1))),str(int(_cv("Unique binding targets",l2)))],
            ["Schema columns",str(int(_cv("Schema columns present",l1))),str(int(_cv("Schema columns present",l2)))],
            ["Shared CV tokens",str(transfer.get("n_shared_cv_tokens","—")),"↑↑"],
            ["New terms needed","—",str(transfer.get("n_new_terms_required","—"))],
            ["Mapping completeness",f"{transfer.get('mapping_completeness',0):.1%}","same"]]
    tbl=ax_F.table(cellText=rows_F[1:],colLabels=rows_F[0],loc="center",cellLoc="center",bbox=[0,0,1,1])
    tbl.auto_set_font_size(False); tbl.set_fontsize(6.5)
    for (r,c),cell in tbl.get_celld().items():
        cell.set_edgecolor("#dddddd"); cell.set_linewidth(0.4)
        if r==0: cell.set_facecolor("#2c3e50"); cell.set_text_props(color="white",fontweight="bold")
        elif r%2==0: cell.set_facecolor("#f7f7f7")
        else: cell.set_facecolor("white")
    ax_F.set_title("Summary Statistics",fontweight="bold",pad=6); _panel_label(ax_F,"F",x=-0.04)
    fig.suptitle(f"Two-Pathway Comparison: {l1} vs {l2}",fontsize=10,fontweight="bold")
    fig.savefig(out); plt.close(fig); print(f"  ✓  Plot saved → {out}")


def plot_domain_function_profile(df1,df2,l1,l2,out):
    def _tag_counts(df):
        ctr=Counter()
        if "function_tags" not in df.columns: return ctr
        for val in df["function_tags"].dropna():
            for tok in split_cell(val):
                lb=extract_label(tok)
                if lb: ctr[lb.lower()]+=1
        return ctr
    c1t,c2t=_tag_counts(df1),_tag_counts(df2)
    all_tags=sorted((c1t|c2t).keys(),key=lambda t:-(c1t.get(t,0)+c2t.get(t,0)))
    top=all_tags[:15]; n1=max(sum(c1t.values()),1); n2=max(sum(c2t.values()),1)
    fig=plt.figure(figsize=(12,5.5))
    gs=fig.add_gridspec(1,2,wspace=0.42,left=0.07,right=0.97,top=0.88,bottom=0.15,width_ratios=[2,3])
    ax_A,ax_B=fig.add_subplot(gs[0]),fig.add_subplot(gs[1])
    tpal=[PAL["blue"],PAL["orange"],PAL["green"],PAL["red"],PAL["purple"],PAL["yellow"],PAL["cyan"],"#997700","#004488","#882255","#BBBBBB","#44AA99","#EE7733","#009988","#CC3311"]
    for yi,(ctr,total) in enumerate([(c1t,n1),(c2t,n2)]):
        left=0.0
        for tag,col in zip(top,tpal):
            frac=ctr.get(tag,0)/total
            if frac==0: continue
            ax_A.barh(1-yi,frac*100,height=0.38,left=left,color=col,edgecolor="white",lw=0.4)
            if frac*100>4: ax_A.text(left+frac*100/2,1-yi,f"{frac*100:.0f}%",ha="center",va="center",fontsize=5.5,color="white",fontweight="bold")
            left+=frac*100
    ax_A.set_yticks([0,1]); ax_A.set_yticklabels([l2,l1],fontsize=8)
    ax_A.set_xlabel("% of domain occurrences"); ax_A.set_title("Function-Tag Composition\n(top 15 tags)",fontweight="bold",pad=6)
    ax_A.set_xlim(0,105); ax_A.spines["left"].set_visible(False); ax_A.tick_params(axis="y",length=0)
    ax_A.legend(handles=[mpatches.Patch(color=c,label=t[:22]) for t,c in zip(top,tpal)],fontsize=5,ncol=2,loc="lower right",handlelength=1.0,handletextpad=0.4,columnspacing=0.8,bbox_to_anchor=(1.0,-0.25))
    _panel_label(ax_A,"A")
    if "domain_type_id" in df1.columns and "domain_type_id" in df2.columns:
        s1=set(df1["domain_type_id"].dropna().unique()); s2=set(df2["domain_type_id"].dropna().unique())
        cnt1=df1["domain_type_id"].value_counts(); cnt2=df2["domain_type_id"].value_counts()
        def _sh(d): return d.split(":")[-1] if ":" in d else d
        sh_d=[d for d in sorted(s1&s2,key=lambda x:-(cnt1.get(x,0)+cnt2.get(x,0)))[:8]]
        o1=[d for d in sorted(s1-s2,key=lambda x:-cnt1.get(x,0))[:6]]
        o2=[d for d in sorted(s2-s1,key=lambda x:-cnt2.get(x,0))[:6]]
        doms=[_sh(d) for d in o1]+[_sh(d) for d in sh_d]+[_sh(d) for d in o2]
        v1s=[cnt1.get(d,0) for d in o1]+[cnt1.get(d,0) for d in sh_d]+[0]*len(o2)
        v2s=[0]*len(o1)+[cnt2.get(d,0) for d in sh_d]+[cnt2.get(d,0) for d in o2]
        cc=[PAL["ref_only"]]*len(o1)+[PAL["shared"]]*len(sh_d)+[PAL["tgt_only"]]*len(o2)
        yB=np.arange(len(doms)); wB=0.33
        ax_B.barh(yB+wB/2,v1s,height=wB,color=cc,edgecolor="white",lw=0.4,alpha=0.85)
        ax_B.barh(yB-wB/2,v2s,height=wB,color=[PAL["tgt_only"] if c==PAL["ref_only"] else PAL["ref_only"] if c==PAL["tgt_only"] else c for c in cc],edgecolor="white",lw=0.4,alpha=0.65)
        ax_B.set_yticks(yB); ax_B.set_yticklabels(doms,fontsize=6.5); ax_B.set_xlabel("Domain occurrence count")
        ax_B.set_title(f"Top Domain Types\nGreen=shared  Blue={l1}-only  Orange={l2}-only",fontweight="bold",pad=6,fontsize=7.5)
        for lbt,mid_y,col in [(f"{l1}-only",len(o1)/2-0.5,PAL["ref_only"]),("Shared",len(o1)+len(sh_d)/2-0.5,PAL["shared"]),(f"{l2}-only",len(o1)+len(sh_d)+len(o2)/2-0.5,PAL["tgt_only"])]:
            ax_B.text(-max(v1s+v2s+[1])*0.18,mid_y,lbt,ha="right",va="center",fontsize=6,color=col,fontweight="bold")
        ax_B.legend(handles=[mpatches.Patch(color=PAL["blue"],label=l1),mpatches.Patch(color=PAL["orange"],label=l2)],fontsize=6.5,loc="lower right",handlelength=1.2)
    _panel_label(ax_B,"B")
    fig.suptitle(f"Domain Functional Profile: {l1}  vs  {l2}",fontsize=9.5,fontweight="bold",y=0.97)
    fig.savefig(out); plt.close(fig); print(f"  ✓  Plot saved → {out}")


# ── Pipeline orchestrator ─────────────────────────────────────────────────────
def run_pipeline(path_camp,path_other,path_ann,out_dir,make_plots=False,label1="cAMP",label2="otherPathway",binding_col_override=None):
    out_dir=Path(out_dir); out_dir.mkdir(parents=True,exist_ok=True)
    print("="*68+"\n  TWO-PATHWAY ONTOLOGY TRANSFER VALIDATION PIPELINE\n"+"="*68)
    print(f"  Pathway 1 (reference) : {Path(path_camp).name}  [{label1}]")
    print(f"  Pathway 2 (target)    : {Path(path_other).name}  [{label2}]")
    print(f"  Annotation registry   : {Path(path_ann).name}")
    print(f"  Output directory      : {out_dir}")
    print_section("1 / 7   Loading inputs")
    df1=load_table(path_camp); print(f"  ✓  {label1}: {df1.shape[0]} rows × {df1.shape[1]} cols")
    df2=load_table(path_other); print(f"  ✓  {label2}: {df2.shape[0]} rows × {df2.shape[1]} cols")
    ann=load_table(path_ann); print(f"  ✓  Annotation registry: {ann.shape[0]} terms")
    cv=build_cv_lookup(ann); print("  ✓  CV lookup: "+" | ".join(f"{k}={len(v)//2}" for k,v in cv.items()))
    print_section("2 / 7   Validating pathway tables")
    bc1=_detect_binding_col(df1,binding_col_override); bc2=_detect_binding_col(df2,binding_col_override)
    print(f"\n── {label1} ──"); res1,inv1,unk1=validate_table(df1,label1,cv,bc1)
    print(f"\n── {label2} ──"); res2,inv2,unk2=validate_table(df2,label2,cv,bc2)
    save_tsv(pd.DataFrame(res1),out_dir/f"validation_report_{label1}.tsv","checks")
    save_tsv(pd.DataFrame(res2),out_dir/f"validation_report_{label2}.tsv","checks")
    ai=pd.DataFrame(inv1+inv2)
    if not ai.empty:
        meta=[c for c in ["pathway","check_id","reason","source_row"] if c in ai.columns]
        save_tsv(ai[meta+[c for c in ai.columns if c not in meta]],out_dir/"invalid_rows.tsv","flagged rows")
    else: save_tsv(pd.DataFrame(),out_dir/"invalid_rows.tsv")
    au=pd.DataFrame(unk1+unk2)
    save_tsv(au if not au.empty else pd.DataFrame(),out_dir/"unknown_terms_report.tsv","unknown tokens")
    print_section("3 / 7   Building missing term report")
    md=build_missing_term_report({label1:df1,label2:df2},cv)
    save_tsv(md,out_dir/"missing_term_report.tsv","unmatched tokens")
    print_section("4 / 7   Computing per-pathway summary metrics")
    m1=compute_pathway_metrics(df1,label1,cv); m2=compute_pathway_metrics(df2,label2,cv)
    save_tsv(pd.DataFrame([m1,m2]),out_dir/"pathway_summary_metrics.tsv","summary")
    print_section("5 / 7   Computing transfer metrics")
    transfer=compute_transfer_metrics(df1,df2,cv,label1,label2)
    save_tsv(pd.DataFrame([transfer]),out_dir/"transfer_metrics.tsv","transfer")
    for k,v in [("Schema reuse rate",f"{transfer['schema_reuse_rate']:.1%}"),("Vocabulary reuse rate",f"{transfer['vocab_reuse_rate']:.1%}"),("Rule portability rate",f"{transfer['rule_portability_rate']:.1%}"),("Mapping completeness",f"{transfer['mapping_completeness']:.1%}"),("Reused CV terms",transfer['n_reused_terms']),("New terms required",transfer['n_new_terms_required']),("Registry terms unused",transfer['n_registry_terms_unused'])]:
        print(f"  {k:<30}: {v}")
    print_section("6 / 7   Computing overlap and Jaccard metrics")
    ov=compute_overlap_metrics(df1,df2,label1,label2)
    save_tsv(ov,out_dir/"overlap_metrics.tsv","Jaccard")
    for _,row in ov.iterrows(): print(f"  {row['category']:<30}  Jaccard={row['jaccard']:.3f}  (shared {row['shared']})")
    comp=build_comparison_report(df1,df2,label1,label2,cv)
    save_tsv(comp,out_dir/"comparison_report.tsv","comparison")
    if make_plots:
        print_section("7 / 7   Generating publication-quality figures")
        plot_validation_summary(res1,res2,out_dir/"validation_summary_plot.png",label1,label2)
        plot_cv_coverage_heatmap(m1,m2,label1,label2,out_dir/"cv_coverage_heatmap.png")
        plot_term_reuse(transfer,out_dir/"term_reuse_barplot.png",label1,label2)
        plot_jaccard_heatmap(ov,out_dir/"jaccard_heatmap.png",label1,label2)
        plot_transfer_radar(transfer,out_dir/"transfer_metrics_radar.png",label1,label2)
        plot_pathway_comparison(comp,ov,transfer,label1,label2,out_dir/"pathway_comparison_plot.png")
        plot_domain_function_profile(df1,df2,label1,label2,out_dir/"domain_function_profile.png")
    else: print_section("7 / 7   Plots skipped  (pass --plots to enable)")
    t1=len(res1); p1=sum(1 for r in res1 if r["status"]=="PASS"); f1=sum(1 for r in res1 if r["status"]=="FAIL")
    t2=len(res2); p2=sum(1 for r in res2 if r["status"]=="PASS"); f2=sum(1 for r in res2 if r["status"]=="FAIL")
    print("\n"+"="*68+"\n  PIPELINE COMPLETE\n"+"="*68)
    print(f"  {label1:<20} {t1} checks | {p1} PASS | {f1} FAIL")
    print(f"  {label2:<20} {t2} checks | {p2} PASS | {f2} FAIL")
    print(f"  Flagged rows  : {len(inv1)+len(inv2)}")
    print(f"  Missing terms : {len(md)}")
    print(f"  All outputs   → {out_dir}/")
    if f1+f2>0: print("\n  ⚠  Some checks FAILED — review invalid_rows.tsv")
    print("="*68)


def main():
    p=argparse.ArgumentParser(description="Two-Pathway Ontology Transfer Validation Pipeline")
    p.add_argument("--camp",required=True); p.add_argument("--other",required=True)
    p.add_argument("--ann",required=True); p.add_argument("--out",default="output_transfer")
    p.add_argument("--label1",default="cAMP"); p.add_argument("--label2",default="otherPathway")
    p.add_argument("--plots",action="store_true"); p.add_argument("--binding_col",default=None)
    a=p.parse_args()
    run_pipeline(Path(a.camp),Path(a.other),Path(a.ann),Path(a.out),a.plots,a.label1,a.label2,a.binding_col)

if __name__=="__main__": main()
