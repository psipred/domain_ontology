"""
convert_toolkit_output.py
──────────────────────────────────────────────────────────────────────────────
Converts the MMseqs2 cluster output downloaded from the Tübingen Bioinformatics
Toolkit into the format expected by stage 13 and stage 16.

Real file format (NOT a two-column TSV — full FASTA sequences per cluster):

    Number of sequences in the input set: 143
    Number of sequences in the reduced set: 88

    Cluster#    1
    >sp|Q13370|PDE3B_HUMAN cGMP-inhibited 3',5'-cyclic phosphodiesterase 3B ...
    MRRDERDAKAM...

    Cluster#    2
    >sp|Q13464|ROCK1_HUMAN Rho-associated protein kinase 1 ...
    MSTGDSFETR...
    >sp|O75116|ROCK2_HUMAN Rho-associated protein kinase 2 ...
    MSRPPPTGKM...
    ...

The first sequence in each cluster block is treated as the representative.
All sequences in the block (including the representative itself) become rows
in the output files.

Output files:
  camp_family_cluster_ids.tsv
  protein_group_labels_sequencing.tsv   <- for stage 13 PCA colouring
──────────────────────────────────────────────────────────────────────────────
"""

from __future__ import annotations
import csv
import re
import time
from pathlib import Path
from collections import Counter

# ── CONFIGURE THESE PATHS ─────────────────────────────────────────────────────
INPUT_TSV      = r"C:\Users\32045\Desktop\测试\integrated_pipeline _calcium_human\input\mmseqs2_2668130_1.tsv"
OUTPUT_CLUSTER = r"C:\Users\32045\Desktop\测试\integrated_pipeline _calcium_human\input\calcium_family_cluster_ids.tsv"
OUTPUT_LABELS  = r"C:\Users\32045\Desktop\测试\integrated_pipeline _calcium_human\input\calcium_protein_group_labels_sequencing.tsv"

FETCH_GENE_NAMES = True

# ── UniProt accession regex ───────────────────────────────────────────────────
_UNIPROT_RE = re.compile(
    r'[OPQ][0-9][A-Z0-9]{3}[0-9]'
    r'|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}'
)


# ── Helpers ───────────────────────────────────────────────────────────────────

def accession_from_fasta_header(header: str) -> str | None:
    """
    Extract UniProt accession from a FASTA header line, e.g.:
      >sp|Q13370|PDE3B_HUMAN cGMP-inhibited ...   ->  'Q13370'
      >tr|A0A000|GENE_HUMAN description           ->  'A0A000'
      >P07550                                     ->  'P07550'
    """
    s = header.lstrip(">").strip()
    parts = s.split("|")
    if len(parts) >= 2:
        candidate = parts[1].strip()
        if _UNIPROT_RE.fullmatch(candidate):
            return candidate
    first = s.split()[0] if s else ""
    if _UNIPROT_RE.fullmatch(first):
        return first
    m = _UNIPROT_RE.search(s)
    return m.group(0) if m else None


def fetch_gene_name(accession: str) -> str:
    """Fetch gene name from UniProt REST API. Returns accession on failure."""
    import urllib.request, json
    url = (f"https://rest.uniprot.org/uniprotkb/{accession}"
           f"?fields=gene_names,protein_name&format=json")
    try:
        with urllib.request.urlopen(url, timeout=10) as r:
            data = json.loads(r.read())
        for g in data.get("genes", []):
            name = g.get("geneName", {}).get("value", "")
            if name:
                return name
        pname = (data.get("proteinDescription", {})
                     .get("recommendedName", {})
                     .get("fullName", {})
                     .get("value", ""))
        if pname:
            return pname.split()[0]
    except Exception:
        pass
    return accession


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    input_path = Path(INPUT_TSV)
    if not input_path.exists():
        print(f"ERROR: {INPUT_TSV} not found.")
        return

    # ── Parse the FASTA-by-cluster format ─────────────────────────────────
    entries: list[tuple[str, int, str]] = []   # (accession, cluster_id, rep)
    current_cluster_id: int | None = None
    current_rep:        str | None = None

    _CLUSTER_LINE = re.compile(r'^Cluster#\s+(\d+)\s*$')

    with open(input_path, encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")

            # Cluster separator: "Cluster#\t1"
            m = _CLUSTER_LINE.match(line)
            if m:
                current_cluster_id = int(m.group(1))
                current_rep        = None
                continue

            # FASTA header: ">sp|ACC|NAME ..."
            if line.startswith(">"):
                if current_cluster_id is None:
                    continue
                acc = accession_from_fasta_header(line)
                if acc is None:
                    print(f"  WARN: cannot parse accession from: {line!r}")
                    continue
                if current_rep is None:
                    current_rep = acc        # first sequence = representative
                entries.append((acc, current_cluster_id, current_rep))
                continue

            # Sequence lines, blank lines, metadata lines — skip

    if not entries:
        print("ERROR: No FASTA entries found in the input file.")
        print("Check that the file contains Cluster# blocks with >sp|...|... headers.")
        return

    n_clusters = len({cid for _, cid, _ in entries})
    n_proteins = len(entries)
    print(f"Parsed : {n_proteins} proteins in {n_clusters} clusters")

    sizes     = Counter(cid for _, cid, _ in entries)
    size_dist = Counter(sizes.values())
    print(f"Size distribution (size: n_clusters): "
          f"{dict(sorted(size_dist.items()))}")
    print()

    # ── Write camp_family_cluster_ids.tsv ─────────────────────────────────
    Path(OUTPUT_CLUSTER).parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_CLUSTER, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["accession", "cluster_id", "representative"],
            delimiter="\t",
        )
        w.writeheader()
        for acc, cid, rep in entries:
            w.writerow({"accession": acc, "cluster_id": cid, "representative": rep})
    print(f"Written: {OUTPUT_CLUSTER}")

    # ── Fetch gene names per cluster representative ────────────────────────
    rep_by_cluster: dict[int, str] = {}
    for acc, cid, rep in entries:
        if cid not in rep_by_cluster:
            rep_by_cluster[cid] = rep

    cluster_labels: dict[int, str] = {}
    if FETCH_GENE_NAMES:
        print(f"\nFetching gene names for {n_clusters} cluster representatives...")
        for cid in sorted(rep_by_cluster):
            rep_acc = rep_by_cluster[cid]
            label   = fetch_gene_name(rep_acc)
            label   = (label.replace(" ", "_")
                            .replace("/", "_")
                            .replace(",", "")
                            .replace("(", "")
                            .replace(")", ""))[:35]
            cluster_labels[cid] = label
            print(f"  Cluster {cid:3d}  rep={rep_acc}  -> {label}")
            time.sleep(0.25)
    else:
        for cid in rep_by_cluster:
            cluster_labels[cid] = f"Cluster_{cid}"
        print("Skipped gene name fetching — using Cluster_N labels")

    # ── Write protein_group_labels_sequencing.tsv ─────────────────────────
    Path(OUTPUT_LABELS).parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_LABELS, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["protein_accession", "protein_group",
                        "cluster_id", "representative"],
            delimiter="\t",
        )
        w.writeheader()
        for acc, cid, rep in entries:
            w.writerow({
                "protein_accession": acc,
                "protein_group":     cluster_labels.get(cid, f"Cluster_{cid}"),
                "cluster_id":        cid,
                "representative":    rep,
            })
    print(f"Written: {OUTPUT_LABELS}")

    print()
    print("=" * 60)
    print("  CONVERSION COMPLETE")
    print("=" * 60)
    print()
    print("Next steps:")
    print()
    print("  1. Copy output files to your input folder:")
    print("       camp_family_cluster_ids.tsv          ->  input/")
    print("       protein_group_labels_sequencing.tsv  ->  input/")
    print()
    print("  2. Add to pipeline_config.yaml:")
    print("       mmseqs2_clusters:     input/camp_family_cluster_ids.tsv")
    print("       protein_group_labels: input/protein_group_labels_sequencing.tsv")
    print()
    print("  3. Run stages 13 and 16:")
    print("       python run_pipeline.py --stages 13,16")
    print()
    print("  NOTE: pairwise_identity.tsv is NOT available from the web toolkit.")
    print("  Stage 16 will use binary cluster-membership distance proxy.")


if __name__ == "__main__":
    main()
