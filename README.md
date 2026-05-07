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
**Verify logical consistency** of the output population.ttl using the HermiT OWL 2 DL reasoner in Protégé 5.6.8 before downstream analysis.

---

### /validation_pipeline — Stage 5: Ontology validation and evaluation

The validation pipeline comprises five scripts that must be run sequentially. It performs enrichment analysis, Mantel testing, clustering comparison, and blind transfer evaluation.

#### 01_background_data.py

Downloads and processes the complete reviewed Swiss-Prot human proteome as the enrichment background. Must be run once before any validation analysis.

**Key output:** `background_domain_counts.tsv`

```bash
python validation_pipeline/01_background_data.py
```

---

#### 02_make_pipeline_inputs.py

Prepares all input files required by the validation pipeline from the populated occurrence tables. This script produces a `*_protein_group_labels.tsv` file that assigns each protein to a functional group (e.g. GPCR, PKA, PDE, adenylyl cyclase).

**Key output:** `*_protein_group_labels.tsv` (requires manual completion)

```bash
python validation_pipeline/02_make_pipeline_inputs.py \
    --input cAMP_human_MERGED_domain_occurrences_WITH_IDS.tsv \
    --output-dir pipeline_inputs/
```

> **Manual step required:** The `*_protein_group_labels.tsv` produced by this script must be manually reviewed and completed before proceeding to `03_fetch_sequences.py`. Open the file, verify that each protein has been assigned to the correct functional group, and fill in any rows where the group label is missing or incorrect.

---

#### 03_fetch_sequences.py

Retrieves canonical FASTA sequences for all human and mouse cAMP proteins from the UniProt REST API and concatenates them into a combined FASTA file for MMseqs2 clustering.

**Key output:** `combined_cAMP_proteins.fasta`

```bash
python validation_pipeline/03_fetch_sequences.py \
    --input pipeline_inputs/ \
    --output combined_cAMP_proteins.fasta
```

> **Manual step required:** After running this script, submit `combined_cAMP_proteins.fasta` to MMseqs2 before proceeding to `04_convert_toolkit_output.py`. See the MMseqs2 parameters section above for the exact command.

---

#### 04_convert_toolkit_output.py

Converts the MMseqs2 cluster output file into the tabular format required by the validation pipeline.

**Key output:** `mmseqs2_clusters_converted.tsv`

```bash
python validation_pipeline/04_convert_toolkit_output.py \
    --mmseqs-output mmseqs2_output_cluster.tsv \
    --fasta combined_cAMP_proteins.fasta \
    --output mmseqs2_clusters_converted.tsv
```

---

#### 05_run_pipeline.py

Runs the complete validation pipeline including domain enrichment analysis (Fisher's exact test with Benjamini-Hochberg correction), permutation testing (1,000 permutations), Mantel test (scikit-bio, 9,999 permutations), ontology-based and sequence-based clustering comparison (ARI and NMI), cross-pathway and cross-species Jaccard comparison, and blind transfer evaluation (Methods A, B, C with bootstrap confidence intervals and permutation-based p-values).

**Key output:** `validation_results/` directory containing all figures and statistical output tables.

```bash
python validation_pipeline/05_run_pipeline.py \
    --input pipeline_inputs/ \
    --clusters mmseqs2_clusters_converted.tsv \
    --background background_domain_counts.tsv \
    --output validation_results/
---

## Cross-pathway comparison — /two_pathway_transfer

Computes Jaccard similarity and five transfer metrics comparing the human cAMP build against the human calcium signalling build. Run after both builds have completed Stages 1–4.

```bash
python two_pathway_transfer/run_two_pathway_transfer.py \
    --reference cAMP_human_MERGED_domain_occurrences_WITH_IDS.tsv \
    --target calcium_human_MERGED_domain_occurrences_WITH_IDS.tsv \
    --output two_pathway_transfer/output_transfer/
```

**Key output:** `two_pathway_transfer/output_transfer/overlap_metrics.tsv` and `transfer_metrics.tsv`

---

## Cross-species comparison — /two_species_transfer

Computes Jaccard similarity, five transfer metrics, domain occurrence frequency correlation (Pearson's r), and PCoA with Procrustes superimposition comparing the human and mouse cAMP builds. Run after both builds have completed Stages 1–4.

```bash
python two_species_transfer/run_two_species_transfer.py \
    --reference cAMP_human_MERGED_domain_occurrences_WITH_IDS.tsv \
    --target cAMP_mouse_MERGED_domain_occurrences_WITH_IDS.tsv \
    --output two_species_transfer/output_transfer/
```

**Key output:** `two_species_transfer/output_transfer/` containing overlap metrics, transfer metrics, occurrence frequency correlation, and PCoA coordinates.

---
---

## Blind transfer evaluation — /blind_test

Evaluates whether domain-architecture matching (Method B) generalises more faithfully than direct ortholog mapping (Method A) or nearest-sequence similarity (Method C) for annotation transfer from human to mouse cAMP proteins. The pipeline comprises eight sequential scripts with two manual curation steps and one external MMseqs2 run.

> **Critical:** Step 5 creates an 80/20 holdout split. The 20% test set (`ortholog_pairs_test.tsv`) must NOT be opened or used until Step 8. Using test set results to tune any parameter invalidates the blind evaluation.

---

### Step 1 — Prepare protein tables

Builds clean normalised protein tables for human and mouse from the ID-injected occurrence TSVs. Assigns each protein to a functional group and family module.

**Key output:** `transfer_pipeline/data/human_proteins.tsv`, `mouse_proteins.tsv`

```bash
python blind_test/01_prepare_protein_tables.py \
    --human cAMP_human_MERGED_domain_occurrences_WITH_IDS.tsv \
    --mouse cAMP_mouse_MERGED_domain_occurrences_WITH_IDS.tsv \
    --out transfer_pipeline/data/
```

Optional group label files can be provided with `--human_groups` and `--mouse_groups` (protein_accession | protein_group TSV format). If omitted, all proteins are assigned to group "other".

---

### Step 2 — Fetch sequences

Retrieves canonical FASTA sequences for all human and mouse proteins from the UniProt REST API and concatenates them into a combined FASTA file for MMseqs2.

**Key output:** `transfer_pipeline/data/combined_camp_proteins.fasta`

```bash
python blind_test/02_fetch_sequences.py \
    --data transfer_pipeline/data/ \
    --out  transfer_pipeline/data/
```

> **Manual step required:** After this script completes, submit `combined_camp_proteins.fasta` to MMseqs2 using the parameters in the MMseqs2 section above before proceeding to Step 3.

---

### Step 3 — Parse MMseqs2 results

Parses the MMseqs2 clustering output and assigns human and mouse proteins to shared clusters. Automatically detects both FASTA-per-cluster format (Tübingen Toolkit default) and two-column TSV format.

**Key output:** `transfer_pipeline/data/mmseqs2_clusters.tsv`

```bash
python blind_test/03_parse_mmseqs2_results.py \
    --clusters mmseqs2_output_cluster.tsv \
    --data transfer_pipeline/data/
```

---

### Step 4 — Build ortholog candidates

Scores all human-mouse protein pairs using four independent signals: MMseqs2 cluster membership, gene symbol similarity, shared protein group, and sequence length ratio. Pairs below a score threshold of 1.0 are excluded.

**Key output:** `transfer_pipeline/data/ortholog_candidates.tsv`

```bash
python blind_test/04_build_ortholog_candidates.py \
    --data  transfer_pipeline/data/ \
    --human cAMP_human_MERGED_domain_occurrences_WITH_IDS.tsv \
    --mouse cAMP_mouse_MERGED_domain_occurrences_WITH_IDS.tsv
```

> **Manual step required:** Open `ortholog_candidates.tsv` and fill in the `curation_status` column for each pair using one of three values: `high_confidence` (clear 1:1 ortholog), `ambiguous` (paralog or unclear), or `excluded` (not orthologs). Save as `ortholog_candidates_curated.tsv`. Then proceed to Step 5.

---

### Step 5 — Finalise ortholog pairs and create blind holdout split

Filters to approved pairs and creates a stratified 80/20 holdout split by confidence tier. The 20% test set is locked and must not be used until Step 8.

**Key outputs:**
- `ortholog_pairs_calibration.tsv` (80% — used in Steps 6 and 7)
- `ortholog_pairs_test.tsv` (20% — LOCKED until Step 8)

```bash
python blind_test/05_finalize_ortholog_pairs.py \
    --data      transfer_pipeline/data/ \
    --test_frac 0.20 \
    --seed      42
```

---

### Step 6 — Build annotation tables

Extracts human and mouse annotations into long-form tables. Uses `ortholog_pairs_calibration.tsv` only — the test set is never touched at this stage. Also writes a `transferability_rules.tsv` file that must be edited before proceeding.

**Key outputs:** `human_annotations_long.tsv`, `mouse_annotations_gold.tsv`, `transferability_rules.tsv`

```bash
python blind_test/06_build_annotation_tables.py \
    --human cAMP_human_MERGED_domain_occurrences_WITH_IDS.tsv \
    --mouse cAMP_mouse_MERGED_domain_occurrences_WITH_IDS.tsv \
    --data  transfer_pipeline/data/
```

> **Manual step required:** Open `transferability_rules.tsv` and confirm which annotation categories should be transferred. Decisions must be justified on biological grounds before running Step 7. Do not consult mouse gold labels when making these decisions.

---

### Step 6b — Build mouse structural annotations

Extracts only structural annotations from the mouse occurrence table (domain type, topology, position, proximity). This is the only mouse data that Method B is allowed to see during prediction, enforcing the blind test protocol.

**Key output:** `transfer_pipeline/data/mouse_annotations_structural.tsv`

```bash
python blind_test/06b_build_mouse_structural_annotations.py \
    --mouse cAMP_mouse_MERGED_domain_occurrences_WITH_IDS.tsv \
    --out   transfer_pipeline/data/
```

---

### Step 7 — Generate transfer predictions

Generates predictions for all three transfer methods using the calibration set and transferability rules.

**Key output:** `transfer_pipeline/data/predictions/`

```bash
python blind_test/07_generate_transfer_predictions.py \
    --data transfer_pipeline/data/
```

---

### Step 8 — Evaluate transfer predictions

Evaluates all three methods against the frozen mouse gold standard using the locked 20% test set. Computes precision, recall and F1 per annotation category with bootstrap confidence intervals (1,000 iterations) and permutation-based p-values (1,000 permutations).

**Key output:** `transfer_pipeline/results/`

```bash
python blind_test/08_evaluate_transfer.py \
    --data transfer_pipeline/data/
```

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
