"""
03_parse_mmseqs2_results.py
────────────────────────────────────────────────────────────────────────────
Step 3 — Parse MMseqs2 clustering results.

Input
─────
  --clusters   two-column TSV: representative TAB member
               (from Tübingen Toolkit or any MMseqs2 output)
  --fasta      combined_camp_proteins.fasta  (to read species= tags)
  --data       data directory (for human/mouse_proteins.tsv fallback)

Output
──────
  mmseqs2_clusters.tsv
    protein_accession | species | cluster_id | representative |
    cluster_size | has_human | has_mouse | cross_species_cluster
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations
import argparse
import re
from collections import Counter, defaultdict
from pathlib import Path
import pandas as pd


def _extract_acc(s: str) -> str:
    s     = s.strip().lstrip(">")
    parts = s.split("|")
    return parts[1].strip() if len(parts) >= 2 else s.split()[0].strip()


def _load_species_map(fasta_path: Path, human_tsv: Path,
                       mouse_tsv: Path) -> dict[str, str]:
    """
    Build {accession: species} and register EVERY possible ID form so the
    mapping works regardless of which format the toolkit uses in the cluster file.

    For each FASTA sequence the following keys are all added to the map:
      O00329              bare UniProt accession
      sp|O00329|PK3CD_HUMAN   full UniProt identifier (no description)
      PK3CD_HUMAN         entry name (third pipe field)
      O00329|PK3CD_HUMAN  no-prefix pipe format

    Species detection order:
      1. species= tag       2. _HUMAN/_MOUSE suffix
      3. OS= organism name  4. OX= NCBI taxonomy ID
      5. Protein TSV fallback
    """
    species_map: dict[str, str] = {}

    if fasta_path.exists():
        with fasta_path.open() as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue

                # Detect species from header
                sp: str | None = None
                m = re.search(r"species=(\w+)", line)
                if m:
                    sp = m.group(1).lower()
                elif "_HUMAN" in line:          sp = "human"
                elif "_MOUSE" in line:          sp = "mouse"
                elif "OS=Homo sapiens" in line: sp = "human"
                elif "OS=Mus musculus" in line: sp = "mouse"
                elif "OX=9606"  in line:        sp = "human"
                elif "OX=10090" in line:        sp = "mouse"

                if sp is None:
                    continue

                # Parse all ID forms from the header
                # Header: >sp|O00329|PK3CD_HUMAN Phosphatidylinositol...
                raw  = line.strip().lstrip(">")
                first_word = raw.split()[0]          # sp|O00329|PK3CD_HUMAN
                parts = first_word.split("|")

                # Register every form the cluster file might use
                if len(parts) == 3:
                    # Standard UniProt: sp|O00329|PK3CD_HUMAN
                    db, acc, entry = parts
                    species_map[acc]                  = sp  # O00329
                    species_map[first_word]           = sp  # sp|O00329|PK3CD_HUMAN
                    species_map[entry]                = sp  # PK3CD_HUMAN
                    species_map[f"{acc}|{entry}"]     = sp  # O00329|PK3CD_HUMAN
                elif len(parts) == 2:
                    # No db prefix: O00329|PK3CD_HUMAN
                    acc, entry = parts
                    species_map[acc]                  = sp
                    species_map[entry]                = sp
                    species_map[first_word]           = sp
                else:
                    # Plain accession
                    species_map[first_word]           = sp

    # Fallback: protein TSV tables from step 01
    for tsv, species in [(human_tsv, "human"), (mouse_tsv, "mouse")]:
        if tsv.exists():
            df    = pd.read_csv(tsv, sep="\t", dtype=str)
            acc_c = next((c for c in ["protein_accession", "accession"]
                          if c in df.columns), None)
            if acc_c:
                for acc in df[acc_c].dropna().str.strip():
                    if acc not in species_map:
                        species_map[acc] = species

    return species_map


def main(args: argparse.Namespace) -> None:
    data_dir   = Path(args.data)
    cluster_file = Path(args.clusters)
    out_dir    = data_dir

    out_dir.mkdir(parents=True, exist_ok=True)
    print("=" * 60)
    print("  03 — Parse MMseqs2 Results")
    print("=" * 60)

    if not cluster_file.exists():
        print(f"\nERROR: {cluster_file} not found.")
        print("Run MMseqs2 on combined_camp_proteins.fasta first:")
        print("  Tübingen Toolkit: https://toolkit.tuebingen.mpg.de/tools/mmseqs2")
        print("  Min seq id = 0.50, Coverage = 0.80")
        print(f"  Download result and save as {cluster_file}")
        return

    print("\nLoading species assignments …")
    species_map = _load_species_map(
        data_dir / "combined_camp_proteins.fasta",
        data_dir / "human_proteins.tsv",
        data_dir / "mouse_proteins.tsv",
    )
    print(f"  {sum(s=='human' for s in species_map.values())} human, "
          f"{sum(s=='mouse' for s in species_map.values())} mouse proteins")

    print(f"\nParsing {cluster_file} …")

    # ── Auto-detect cluster file format ──────────────────────────────────────
    # The Tübingen Toolkit exports two formats depending on which button you use:
    #   Format A — FASTA-per-cluster  (what this file is):
    #     "Number of sequences in the input set: 277"
    #     "Cluster#\t1"
    #     ">sp|P97490|ADCY8_MOUSE ... "
    #     "MELSDVHCLSG..."
    #     ""
    #     "Cluster#\t2"
    #     ...
    #
    #   Format B — Two-column TSV  (representative TAB member):
    #     "O00329\tO00459"
    #     "O00329\tA2ARP1"
    #     ...
    #
    # We detect by checking the first non-empty content line.

    with cluster_file.open() as fh:
        first_content = next(
            (l.strip() for l in fh
             if l.strip() and not l.startswith("Number of")),
            ""
        )

    is_fasta_cluster = first_content.startswith("Cluster#") or first_content.startswith(">")
    rows: list[dict] = []

    if is_fasta_cluster:
        # ── Format A: FASTA-per-cluster ───────────────────────────────────────
        # Each cluster block:
        #   Cluster#\tN
        #   >sp|ACC|ENTRY description
        #   SEQUENCE...
        #   (blank line)
        # Representative = first FASTA entry in each block.
        print("  Detected format: FASTA-per-cluster (Tübingen toolkit default)")
        cluster_id    = 0
        cluster_members: list[str] = []

        def _flush(cid, members, rows):
            if not members:
                return
            rep = members[0]   # first listed = representative
            for mem in members:
                rows.append({"accession": mem, "representative": rep,
                              "cluster_id": cid})

        with cluster_file.open() as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith("Cluster#"):
                    _flush(cluster_id, cluster_members, rows)
                    cluster_id += 1
                    cluster_members = []
                elif line.startswith(">"):
                    # Extract accession from full UniProt header
                    acc = _extract_acc(line)
                    # Also register species directly from this header
                    # (more reliable than looking up in species_map later)
                    sp: str | None = None
                    if "_HUMAN" in line:          sp = "human"
                    elif "_MOUSE" in line:        sp = "mouse"
                    elif "OX=9606" in line:       sp = "human"
                    elif "OX=10090" in line:      sp = "mouse"
                    elif "OS=Homo sapiens" in line: sp = "human"
                    elif "OS=Mus musculus" in line: sp = "mouse"
                    if sp and acc not in species_map:
                        species_map[acc] = sp
                    cluster_members.append(acc)
            # flush last cluster
            _flush(cluster_id, cluster_members, rows)

    else:
        # ── Format B: two-column TSV ──────────────────────────────────────────
        print("  Detected format: two-column TSV (representative TAB member)")
        rep_to_id: dict[str, int] = {}
        with cluster_file.open() as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t") if "\t" in line else line.split()
                if len(parts) < 2:
                    continue
                rep_raw, mem_raw = parts[0], parts[1]
                if rep_raw.lower() in ("representative","cluster","rep"):
                    continue
                rep = _extract_acc(rep_raw)
                mem = _extract_acc(mem_raw)
                if rep not in rep_to_id:
                    rep_to_id[rep] = len(rep_to_id) + 1
                rows.append({"accession": mem, "representative": rep,
                              "cluster_id": rep_to_id[rep]})

    if not rows:
        print("ERROR: no data rows found. Check file format.")
        return

    df = pd.DataFrame(rows)
    df["species"] = df["accession"].map(lambda a: species_map.get(a, "unknown"))

    sizes = df.groupby("cluster_id").size().rename("cluster_size")
    df    = df.merge(sizes, on="cluster_id")

    clust_species = (df.groupby("cluster_id")["species"]
                       .apply(set).rename("species_set"))
    df = df.merge(clust_species, on="cluster_id")
    df["has_human"]             = df["species_set"].apply(lambda s: "human" in s)
    df["has_mouse"]             = df["species_set"].apply(lambda s: "mouse" in s)
    df["cross_species_cluster"] = df["has_human"] & df["has_mouse"]
    df = df.drop(columns=["species_set"])
    df = df.rename(columns={"accession": "protein_accession"})
    df = df[["protein_accession","species","cluster_id","representative",
             "cluster_size","has_human","has_mouse","cross_species_cluster"]]

    out_path = out_dir / "mmseqs2_clusters.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    print(f"\nSaved → {out_path}")

    n_cross = df[df["cross_species_cluster"]]["cluster_id"].nunique()
    n_human = df[df["species"]=="human"]["protein_accession"].nunique()
    n_mouse = df[df["species"]=="mouse"]["protein_accession"].nunique()
    n_unk   = df[df["species"]=="unknown"]["protein_accession"].nunique()

    print(f"\nSummary:")
    print(f"  Total clusters:         {df['cluster_id'].nunique()}")
    print(f"  Cross-species clusters: {n_cross}  ← ortholog candidates")
    print(f"  Human proteins:         {n_human}")
    print(f"  Mouse proteins:         {n_mouse}")
    if n_unk:
        print(f"  Unknown species:        {n_unk}  ← should be 0; check FASTA path")
    # fix np.int64 keys showing as np.int64(1) in display
    sz_dist = {int(k): int(v)
               for k, v in Counter(df.groupby("cluster_id").size().values).items()}
    print(f"  Cluster size dist:      {dict(sorted(sz_dist.items()))}")

    if n_cross == 0:
        print("\n  WARNING: 0 cross-species clusters.")
        print(f"  Species map: {n_human} human, {n_mouse} mouse, {n_unk} unknown.")
        if n_unk > 0:
            print("  The combined FASTA was found but species could not be assigned.")
            print("  Check that your combined FASTA is at:")
            print(f"    {data_dir / 'combined_camp_proteins.fasta'}")
        elif n_human == 0 or n_mouse == 0:
            print("  One species has 0 proteins — cluster file may not contain")
            print("  both species, or the FASTA / protein TSVs are missing.")

    print("\nStep 3 complete.  Next: python 04_build_ortholog_candidates.py")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Parse MMseqs2 clustering results")
    p.add_argument("--clusters", required=True,
                   help="Two-column TSV: representative TAB member")
    p.add_argument("--data", default="transfer_pipeline/data",
                   help="Data directory (for FASTA and protein TSVs)")
    main(p.parse_args())
