# Cross-Species 5'UTR RNA Structure Conservation Analysis

**Date:** 2026-04-13
**Status:** Draft

## Goal

Identify conserved RNA structural motifs across 12 eukaryotic species by scanning GraphClust CM models from each species against the 5'UTR sequences of all other species. Classify discovered motifs as known (Rfam) or novel, and quantify their phylogenetic conservation.

## Species

| Category | Species |
|----------|---------|
| Mammals | Human, Mouse, Cow, Bat, Platypus |
| Birds | Guineafowl |
| Fish | Zebrafish |
| Insects | FruitFly |
| Plants | Arabidopsis, Corn, Rice |
| Fungi | Yeast |

## Architecture

Three phases, clearly separated:

### Phase 1: Scan (zalkon)

**Script:** `my_module/cross_species_scan.py`

For each of the 12 source species, collect unique CM models (latest round per cluster, with >0 base pairs) and run `cmsearch` against each of the 11 target species' FASTA files.

- **Total:** 12 × 11 = 132 source→target directions
- **CMs per species:** ~80 unique (after deduplication by cluster_num)
- **Total cmsearch calls:** ~960 CMs × 11 targets ≈ 10,560

**cmsearch parameters:**
```
cmsearch --cpu 2 --noali --tblout <output> -T 25 --toponly <cm> <fasta>
```

- `-T 25`: minimum bitscore threshold
- `--toponly`: forward strand only (5'UTR is sense strand)
- `--noali`: no alignment output (speed)
- `--cpu 2`: 2 threads per cmsearch call

**Parallelisation:** 3 source species concurrently × 8 threads each = 24 threads on zalkon.

**Output structure:**
```
data/cross_species/5UTR/{source_species}/{source}_{cluster}_to_{target}.tblout
```

**Estimated runtime:** 4–6 hours on zalkon (24 threads).

### Phase 2: Aggregate (zalkon or local)

**Function:** `aggregate_cross_species_results()` in `graphclust_results.py`

Parse all tblout files, filter by E-value ≤ 0.001 and bitscore ≥ 25, combine into a single Parquet file.

**Output:** `data/cross_species/5UTR/all_hits.parquet`

**Columns:**
- `source_species`, `source_cluster`, `source_cluster_num`
- `target_species`, `target_name`, `target_from`, `target_to`, `strand`
- `score`, `evalue`

### Phase 3: Analysis (local Jupyter notebook)

**Notebook:** `cross_species_5UTR.ipynb`

Three main analyses:

#### A) Conservation Summary

Per source cluster: count how many target species have significant hits. Rank by conservation degree.

**Output:** Polars DataFrame with `source_species`, `source_cluster`, `n_species_with_hits`, `species_list`, `total_hits`, `mean_score`.

#### B) Bidirectional Motif Clustering

Group clusters from different species that detect each other bidirectionally into motif families. Algorithm:

1. Build a graph: nodes = (species, cluster) pairs; edges = bidirectional hit relationships
2. Two nodes are connected if cluster A from species X has hits in species Y, AND any cluster B from species Y has overlapping hits in species X at the same genomic locus
3. Connected components = motif families
4. Implementation: Union-Find (disjoint set) over bidirectional hit pairs

**Output:** DataFrame mapping each (species, cluster) to a `motif_family_id`, with family-level statistics (n_species, species list, mean SCI, mean bitscore).

#### C) Novel vs. Rfam Classification

Cross-reference motif families against existing Rfam cmscan results (already available in `data/{species}/protein_coding/cmscan_results_specific_library/mRNA/`).

- For each motif family, check if any member cluster's hit regions overlap with Rfam annotations
- Families with Rfam overlap → annotated (report which Rfam family)
- Families without Rfam overlap → novel candidates

**Output:** DataFrame with `motif_family_id`, `rfam_family` (or "novel"), `n_species`, `conservation_score`.

## Visualisations

All plots as PDF, predTED style, no titles, Okabe-Ito palette for categories.

### A) Presence/Absence Heatmap

- Rows: 12 species, ordered by phylogenetic relationship (Mammals → Birds → Fish → Insects → Plants → Fungi)
- Columns: motif families, ordered by conservation degree (most conserved left)
- Colour: best bitscore per cell (viridis scale), grey = absent
- Annotation: Rfam families marked with symbols

### B) Conservation Distribution

- X-axis: number of species with hits (1–12)
- Y-axis: number of motif families
- Stacked bars: Rfam-annotated vs. novel

### C) Taxonomic Conservation

- Group motif families by phylogenetic pattern:
  - Universal (≥10 species)
  - Eukaryotic-wide (6–9 species across kingdoms)
  - Kingdom-specific (e.g. only mammals, only plants)
  - Species-specific (1–2 species)
- Bar chart or Euler diagram showing distribution

### D) Novel Motif Candidates

- Table of top novel motif families sorted by conservation × structural quality
- Columns: motif_family_id, n_species, species list, mean bitscore, mean SCI, representative structure (dot-bracket)

## File Locations

| File | Purpose |
|------|---------|
| `my_module/cross_species_scan.py` | Scan script for zalkon |
| `my_module/graphclust_results.py` | Aggregate + analysis functions |
| `data/cross_species/5UTR/` | Results directory |
| `data/cross_species/5UTR/all_hits.parquet` | Aggregated hits |
| `cross_species_5UTR.ipynb` | Analysis notebook |

## Dependencies

- Infernal (`cmsearch`) from rnaclust conda environment
- Polars, tqdm, matplotlib (from agat conda environment)
- Existing: `my_module.graphclust_results`, `my_module.rfam_overlap`
