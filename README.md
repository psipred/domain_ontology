# Domain-Centric OWL Ontology Pipeline
## A Three-Namespace Controlled Vocabulary for Functional Annotation Transfer in cAMP Signalling

This repository contains the full pipeline for constructing, populating, and validating a domain-centric OWL 2 DL knowledge graph applied to the human cAMP signalling pathway, with cross-pathway (human calcium signalling) and cross-species (mouse cAMP) builds.

---

## Repository structure

```
README.md
01_data_and_merge.py          # Stage 1+2: data collection and occurrence modelling
02_id_injection.py            # Stage 3: controlled vocabulary ID injection
03_owl_population.py          # Stage 4: OWL/Turtle knowledge graph construction
/Annotation_table             # 211-term DOP/DOF/DOT controlled vocabulary
/Gene_lists                   # Input gene symbol lists for each pathway/species build
/blind_test                   # Blind transfer evaluation scripts (Methods A, B, C)
/two_pathway_transfer         # Cross-pathway comparison pipeline (cAMP vs calcium)
/two_species_transfer         # Cross-species comparison pipeline (human vs mouse cAMP)
/validation_pipeline          # Pipeline for the ontology validation
```

---

## Pipeline overview

The pipeline runs in three sequential stages. Each script must be run in order as each stage depends on the output of the previous stage.

```
Stage 1+2:  01_data_and_merge.py
            - output: *_MERGED_domain_occurrences.tsv
            [Manual curation of function_tags, binding_partner, and reference columns]

Stage 3:    02_id_injection.py
              -output: *_MERGED_domain_occurrences_WITH_IDS.tsv
                
Stage 4:    03_owl_population.py
            -output: population.ttl  (OWL/Turtle knowledge graph)
                
Stage 5:    /validation_pipeline
            ├── 01_background_data.py          # builds Swiss-Prot background for enrichment
            ├── 02_make_pipeline_inputs.py     # prepares inputs; produces *_protein_group_labels.tsv which must be manually completed before proceeding
            ├── 03_fetch_sequences.py          # retrieves canonical FASTA sequences
            │     -> [Run MMseqs2 clustering on combined FASTA
            ├── 04_convert_toolkit_output.py   # converts MMseqs2 cluster output to pipeline format
            └── 05_run_pipeline.py             # runs full validation including enrichment,
                                               # Mantel test, clustering comparison,
                                               # and blind transfer evaluation
            -output: validation_results/

```

---

## Script descriptions

### 01_data_and_merge.py — Stage 1 and 2: Data collection and domain occurrence modelling

Queries the UniProt Swiss-Prot REST API and the InterPro REST API to retrieve per-protein domain occurrence data for a given pathway gene list. Consolidates overlapping InterPro signatures into canonical domain occurrences using a longest-span representative strategy with member-database priority (Pfam > SMART > CDD > PROSITE > CATH/Gene3D). Assigns positional category, proximity category, topology context, and copy number to each occurrence. Merges all per-gene outputs into pathway-level TSV occurrence tables.

**Key output:** `{pathway_name}_MERGED_domain_occurrences.tsv`

**Usage examples:**
```bash
# Human cAMP pathway
python 01_data_and_merge.py \
    --gene-file Gene_lists/cAMP_human.txt \
    --pathway-name cAMP_human \
    --organism-id 9606 \
    --fetch-interpro

# Mouse cAMP pathway
python 01_data_and_merge.py \
    --gene-file Gene_lists/cAMP_human.txt \
    --pathway-name cAMP_mouse \
    --organism-id 10090 \
    --fetch-interpro

# Human calcium signalling
python 01_data_and_merge.py \
    --gene-file Gene_lists/calcium_human.txt \
    --pathway-name calcium_human \
    --organism-id 9606 \
    --fetch-interpro

# Single protein (debug)
python 01_data_and_merge.py \
    --name PDE4D \
    --query "gene:PDE4D* AND organism_id:9606" \
    --fetch-interpro
```

> **Note:** After running Stage 1, several columns in the output TSV require manual curation before proceeding to Stage 3: `function_tags` (DOF vocabulary) and `binding_partner` (DOT vocabulary). Assign the most specific applicable term from the annotation table supported by at least one literature or database source.

---

### 02_id_injection.py — Stage 3: Controlled vocabulary identifier injection

Reads the DOP/DOF/DOT annotation table and the merged domain occurrence TSV produced by Stage 1. Replaces every plain-text vocabulary token with its stable prefixed identifier-label string. For example:

- `Cterminal` → `DOP:000003 Cterminal`
- `kinase` → `DOF:000039 kinase`
- `cAMP` → `DOT:000002 cAMP`

Tokens not found in the annotation table are written to `unmapped_tokens.txt` grouped by column. Resolve all unmapped tokens iteratively until the file is empty before proceeding to Stage 4.

**Key output:** `{input_name}_WITH_IDS.tsv`

**Usage examples:**
```bash
# Minimum required arguments
python 02_id_injection.py \
    --annotation Annotation_table/Annotation_table_human_cAMP.xlsx \
    --input cAMP_human_MERGED_domain_occurrences.tsv

# Carry manual annotations over from a previously corrected file
python 02_id_injection.py \
    --annotation Annotation_table/Annotation_table_human_cAMP.xlsx \
    --input cAMP_human_MERGED_domain_occurrences.tsv \
    --carry-over previous_cAMP_occurrences.tsv

# Strict mode: exit with error if any token is unmapped
python 02_id_injection.py \
    --annotation Annotation_table/Annotation_table_human_cAMP.xlsx \
    --input cAMP_human_MERGED_domain_occurrences.tsv \
    --strict
```

---

### 03_owl_population.py — Stage 4: OWL/Turtle knowledge graph construction

Converts the ID-injected domain occurrence TSV into an OWL 2 DL Turtle population file aligned to the core ontology schema. For every row, creates three kinds of OWL individuals: Protein, DomainType, and DomainOccurrence. Also builds CoOccurrencePattern and DomainOrderPattern individuals encoding cross-protein structural relationships. DomainOccurrence individuals are minted with deterministic IRIs of the form `base#{InterPro_accession}_{protein_accession}_{start}_{end}`.

**Key output:** `output/stage4_owl_population/population.ttl`

**Usage:**

Configure the input and output paths in `pipeline_config.yaml`, then run:
```bash
python 03_owl_population.py
```

**Verify logical consistency** of the output population.ttl using the HermiT OWL 2 DL reasoner in Protégé 5.6.8 before downstream analysis.

---

## Database versions

The results reported in the paper used the following database releases. To exactly reproduce the reported results, use these versions.

| Database | Release | Date retrieved |
|---|---|---|
| UniProt Swiss-Prot | 2025_11 | November 2025 |
| InterPro | 109.0 | November 2025 |

---

## Requirements

Python 3.9 or higher is required. Install all dependencies using:

```bash
pip install -r requirements.txt
```

Core dependencies:

| Package | Version | Used in |
|---|---|---|
| pandas | >= 1.5.0 | All scripts |
| requests | >= 2.28.0 | 01_data_and_merge.py |
| rdflib | 6.3.2 | 03_owl_population.py |
| pyyaml | >= 6.0 | 03_owl_population.py |
| openpyxl | >= 3.0.0 | 02_id_injection.py |

---
