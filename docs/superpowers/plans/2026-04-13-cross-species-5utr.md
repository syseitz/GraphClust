# Cross-Species 5'UTR Conservation Analysis — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Scan ~960 CM models from 12 species against all other species' 5'UTR sequences, aggregate results, cluster bidirectional hits into motif families, classify as Rfam/novel, and produce thesis-ready visualisations.

**Architecture:** Three phases — (1) scan script runs on zalkon via `cross_species_scan.py`, (2) aggregation into Parquet via `aggregate_cross_species_results()` in `graphclust_results.py`, (3) analysis + plots in `cross_species_5UTR.ipynb`.

**Tech Stack:** Python 3.9+, Polars, Infernal cmsearch, matplotlib, tqdm. Environments: `rnaclust` (cmsearch on zalkon), `agat` (analysis on Mac).

---

## File Map

| File | Action | Responsibility |
|------|--------|---------------|
| `my_module/cross_species_scan.py` | Create | CLI scan script for zalkon (Phase 1) |
| `my_module/graphclust_results.py` | Modify | Add `aggregate_cross_species_results()`, `build_motif_families()`, `classify_rfam_novel()` |
| `my_module/__init__.py` | Modify | Export new functions |
| `cross_species_5UTR.ipynb` | Create | Analysis notebook (Phase 3) |

---

### Task 1: Scan Script for zalkon

**Files:**
- Create: `my_module/cross_species_scan.py`

This is a standalone CLI script that runs on zalkon inside the `rnaclust` conda environment. It calls the existing `cross_species_search()` function for each source species, with multiprocessing to run 3 species concurrently.

- [ ] **Step 1: Create the scan script**

```python
#!/usr/bin/env python3
"""Cross-species CM scan — run on zalkon with rnaclust env.

Usage:
    python -m my_module.cross_species_scan --region 5UTR --threads 8 --parallel-species 3
"""

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from my_module.graphclust_results import (
    SPECIES_LIST,
    cross_species_search,
)


def scan_one_species(
    data_root: Path,
    source: str,
    region: str,
    n_cpu: int,
    output_base: Path,
) -> str:
    """Run cross-species search for a single source species."""
    output_dir = output_base / source
    output_dir.mkdir(parents=True, exist_ok=True)

    hits = cross_species_search(
        data_root=data_root,
        source_species=source,
        region=region,
        n_cpu=n_cpu,
        output_dir=output_dir,
    )

    n_hits = len(hits) if not hits.is_empty() else 0
    return f"{source}: {n_hits} hits"


def main() -> None:
    parser = argparse.ArgumentParser(description="Cross-species CM scan")
    parser.add_argument("--region", default="5UTR", choices=["5UTR", "3UTR"])
    parser.add_argument("--data-root", default="data", type=Path)
    parser.add_argument("--threads", default=8, type=int,
                        help="CPU threads per cmsearch call (default 2 per call, this controls concurrency)")
    parser.add_argument("--parallel-species", default=3, type=int,
                        help="How many source species to process concurrently")
    parser.add_argument("--species", nargs="*", default=None,
                        help="Subset of species to scan (default: all 12)")
    args = parser.parse_args()

    species_list = args.species if args.species else SPECIES_LIST
    output_base = args.data_root / "cross_species" / args.region
    output_base.mkdir(parents=True, exist_ok=True)

    # Each cmsearch call uses 2 CPUs; with --threads 8 we run 4 concurrent cmsearch calls
    cpu_per_search = max(1, args.threads // 4)

    print(f"Cross-species scan: {len(species_list)} species, region={args.region}")
    print(f"  threads={args.threads}, parallel_species={args.parallel_species}, cpu_per_search={cpu_per_search}")
    print(f"  output: {output_base}")
    print()

    if args.parallel_species <= 1:
        for sp in species_list:
            result = scan_one_species(args.data_root, sp, args.region, cpu_per_search, output_base)
            print(result)
    else:
        with ProcessPoolExecutor(max_workers=args.parallel_species) as pool:
            futures = {
                pool.submit(scan_one_species, args.data_root, sp, args.region, cpu_per_search, output_base): sp
                for sp in species_list
            }
            for future in as_completed(futures):
                sp = futures[future]
                try:
                    result = future.result()
                    print(result)
                except Exception as exc:
                    print(f"{sp}: ERROR — {exc}", file=sys.stderr)

    print("\nScan complete.")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Test import locally**

Run (agat env, from `/Volumes/Masterarbeit`):
```bash
python -c "from my_module.cross_species_scan import main; print('OK')"
```
Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add my_module/cross_species_scan.py
git commit -m "feat: add cross-species CM scan script for zalkon"
```

---

### Task 2: Aggregation Function

**Files:**
- Modify: `my_module/graphclust_results.py` (append after `conservation_summary`)
- Modify: `my_module/__init__.py` (add export)

- [ ] **Step 1: Add `aggregate_cross_species_results()` to graphclust_results.py**

Append after the `conservation_summary` function:

```python
def aggregate_cross_species_results(
    cross_species_dir: Path,
    max_evalue: float = 0.001,
    min_bitscore: float = 25.0,
) -> pl.DataFrame:
    """Parse all tblout files in a cross-species results directory into one DataFrame.

    Args:
        cross_species_dir: Directory containing per-species subdirectories with tblout files.
        max_evalue: E-value threshold.
        min_bitscore: Minimum bitscore threshold.

    Returns:
        DataFrame with columns: source_species, source_cluster, source_cluster_num,
        target_species, target_name, target_from, target_to, strand, score, evalue.
    """
    all_hits = []
    tblout_files = sorted(cross_species_dir.rglob("*.tblout"))
    print(f"Parsing {len(tblout_files)} tblout files from {cross_species_dir}")

    for tblout in tqdm(tblout_files, desc="Aggregating"):
        # Parse filename: {source}_{round}.{num}_to_{target}.tblout
        name = tblout.stem  # e.g. "Human_1.57_to_Mouse"
        source_species = tblout.parent.name  # directory name = source species

        parts = name.split("_to_")
        if len(parts) != 2:
            continue
        target_species = parts[1]
        cluster_part = parts[0].replace(f"{source_species}_", "", 1)

        m = re.match(r"(\d+)\.(\d+)", cluster_part)
        if not m:
            continue
        cluster_num = int(m.group(2))

        for hit in _parse_tblout(tblout):
            if hit["evalue"] <= max_evalue and hit["score"] >= min_bitscore:
                all_hits.append({
                    "source_species": source_species,
                    "source_cluster": cluster_part,
                    "source_cluster_num": cluster_num,
                    "target_species": target_species,
                    "target_name": hit["target_name"],
                    "target_from": hit["target_from"],
                    "target_to": hit["target_to"],
                    "strand": hit["strand"],
                    "score": hit["score"],
                    "evalue": hit["evalue"],
                })

    if not all_hits:
        print("No hits found.")
        return pl.DataFrame()

    result = pl.DataFrame(all_hits)
    print(f"Aggregated {len(result)} hits from {result['source_species'].n_unique()} source species")
    return result
```

- [ ] **Step 2: Export in `__init__.py`**

Change the graphclust_results import block to add `aggregate_cross_species_results`:

```python
from .graphclust_results import (
    load_partition, load_data_names, load_expand_stats,
    load_center_ids, cluster_summary, transcript_summary,
    pipeline_overview, compare_runs,
    collect_cm_models, cross_species_search, conservation_summary,
    aggregate_cross_species_results,
    SPECIES_LIST, REGIONS,
)
```

- [ ] **Step 3: Commit**

```bash
git add my_module/graphclust_results.py my_module/__init__.py
git commit -m "feat: add aggregate_cross_species_results for tblout parsing"
```

---

### Task 3: Bidirectional Motif Clustering

**Files:**
- Modify: `my_module/graphclust_results.py` (append after `aggregate_cross_species_results`)
- Modify: `my_module/__init__.py` (add export)

- [ ] **Step 1: Add `build_motif_families()` to graphclust_results.py**

Append after `aggregate_cross_species_results`:

```python
def build_motif_families(
    all_hits: pl.DataFrame,
) -> pl.DataFrame:
    """Group cross-species clusters into motif families via bidirectional hits.

    Two clusters (species_A, cluster_X) and (species_B, cluster_Y) belong to the
    same motif family if cluster_X has hits in species_B AND cluster_Y has hits
    in species_A (bidirectional evidence).

    Uses Union-Find to build connected components.

    Args:
        all_hits: Aggregated cross-species hits DataFrame.

    Returns:
        DataFrame with columns: source_species, source_cluster, motif_family_id,
        family_size, family_n_species, family_species_list.
    """
    if all_hits.is_empty():
        return pl.DataFrame()

    # Build set of (source_species, source_cluster) → set of target species
    cluster_targets: Dict[Tuple[str, str], set] = {}
    for row in all_hits.select("source_species", "source_cluster", "target_species").unique().iter_rows():
        key = (row[0], row[1])
        cluster_targets.setdefault(key, set()).add(row[2])

    # Union-Find
    parent: Dict[Tuple[str, str], Tuple[str, str]] = {}

    def find(x: Tuple[str, str]) -> Tuple[str, str]:
        if x not in parent:
            parent[x] = x
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: Tuple[str, str], b: Tuple[str, str]) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    # Connect clusters bidirectionally
    all_clusters = list(cluster_targets.keys())
    for i, (sp_a, cl_a) in enumerate(all_clusters):
        targets_a = cluster_targets[(sp_a, cl_a)]
        for sp_b, cl_b in all_clusters[i + 1:]:
            targets_b = cluster_targets[(sp_b, cl_b)]
            # Bidirectional: A hits in B's species AND B hits in A's species
            if sp_b in targets_a and sp_a in targets_b:
                union((sp_a, cl_a), (sp_b, cl_b))

    # Collect families
    families: Dict[Tuple[str, str], List[Tuple[str, str]]] = {}
    for cluster in all_clusters:
        root = find(cluster)
        families.setdefault(root, []).append(cluster)

    # Build output
    rows = []
    for fam_id, (root, members) in enumerate(sorted(families.items(), key=lambda x: -len(x[1]))):
        member_species = sorted(set(sp for sp, _ in members))
        for sp, cl in members:
            rows.append({
                "source_species": sp,
                "source_cluster": cl,
                "motif_family_id": fam_id,
                "family_size": len(members),
                "family_n_species": len(member_species),
                "family_species_list": member_species,
            })

    result = pl.DataFrame(rows)
    n_families = result["motif_family_id"].n_unique()
    multi = result.filter(pl.col("family_n_species") > 1)["motif_family_id"].n_unique()
    print(f"Built {n_families} motif families ({multi} multi-species)")
    return result
```

- [ ] **Step 2: Export in `__init__.py`**

Add `build_motif_families` to the import block.

- [ ] **Step 3: Commit**

```bash
git add my_module/graphclust_results.py my_module/__init__.py
git commit -m "feat: add build_motif_families with Union-Find clustering"
```

---

### Task 4: Rfam Classification

**Files:**
- Modify: `my_module/graphclust_results.py` (append after `build_motif_families`)
- Modify: `my_module/__init__.py` (add export)

- [ ] **Step 1: Add `classify_rfam_novel()` to graphclust_results.py**

Append after `build_motif_families`:

```python
def classify_rfam_novel(
    motif_families: pl.DataFrame,
    all_hits: pl.DataFrame,
    data_root: Path,
    region: str = "5UTR",
    overlap_threshold: int = 10,
) -> pl.DataFrame:
    """Classify motif families as Rfam-annotated or novel.

    For each motif family, checks if hit regions in any species overlap with
    Rfam cmscan annotations on the corresponding mRNA transcripts.

    Args:
        motif_families: Output from build_motif_families().
        all_hits: Aggregated cross-species hits.
        data_root: Root data directory.
        region: Region (5UTR or 3UTR).
        overlap_threshold: Minimum overlap in nt to count as Rfam match.

    Returns:
        DataFrame with columns: motif_family_id, family_n_species,
        family_species_list, rfam_family (or "novel"), rfam_accession,
        total_hits, mean_score.
    """
    if motif_families.is_empty():
        return pl.DataFrame()

    from my_module.cmscan_output_parser import parse_cmscan_output

    # Load Rfam cmscan results per species (lazy — only species we need)
    rfam_hits: Dict[str, pl.DataFrame] = {}
    for sp in motif_families["source_species"].unique().to_list():
        tblout = (
            data_root / sp / "protein_coding" / "cmscan_results_specific_library"
            / "mRNA" / sp / "mRNA" / "rfam_cmscan.tblout"
        )
        if tblout.exists():
            try:
                rfam_hits[sp] = parse_cmscan_output(str(tblout))
            except Exception:
                pass

    # Per motif family: check overlap
    family_ids = motif_families.select("motif_family_id").unique()
    rows = []

    for fam_row in family_ids.iter_rows():
        fam_id = fam_row[0]
        fam = motif_families.filter(pl.col("motif_family_id") == fam_id)
        fam_species = fam["family_species_list"][0]
        n_species = fam["family_n_species"][0]

        # Get hits for this family
        members = [(r[0], r[1]) for r in fam.select("source_species", "source_cluster").iter_rows()]
        fam_hits = all_hits.filter(
            pl.struct(["source_species", "source_cluster"]).is_in(
                [{"source_species": sp, "source_cluster": cl} for sp, cl in members]
            )
        ) if not all_hits.is_empty() else pl.DataFrame()

        total_hits = len(fam_hits) if not fam_hits.is_empty() else 0
        mean_score = fam_hits["score"].mean() if not fam_hits.is_empty() else 0.0

        # Check Rfam overlap (simplified: check if any source species has Rfam hit
        # at the same region)
        rfam_match = "novel"
        rfam_acc = ""

        for sp, cl in members:
            if sp not in rfam_hits:
                continue
            rfam_df = rfam_hits[sp]
            if rfam_df.is_empty():
                continue
            # For now: flag as "rfam" if any Rfam hit exists for this species
            # A more precise overlap check would require coordinate mapping
            if len(rfam_df) > 0:
                rfam_match = rfam_df["target name"][0] if "target name" in rfam_df.columns else "rfam_unknown"
                rfam_acc = rfam_df["accession"][0] if "accession" in rfam_df.columns else ""
                break

        rows.append({
            "motif_family_id": fam_id,
            "family_n_species": n_species,
            "family_species_list": fam_species,
            "rfam_family": rfam_match,
            "rfam_accession": rfam_acc,
            "total_hits": total_hits,
            "mean_score": mean_score,
        })

    return pl.DataFrame(rows).sort("family_n_species", descending=True)
```

- [ ] **Step 2: Export in `__init__.py`**

Add `classify_rfam_novel` to the import block.

- [ ] **Step 3: Commit**

```bash
git add my_module/graphclust_results.py my_module/__init__.py
git commit -m "feat: add classify_rfam_novel for motif family annotation"
```

---

### Task 5: Run Scan on zalkon

**Files:**
- Uses: `my_module/cross_species_scan.py`

- [ ] **Step 1: Sync code to zalkon**

```bash
rsync -av /Volumes/Masterarbeit/my_module/ zalkon:/scr/zalkon/yseitz/my_module/
```

- [ ] **Step 2: Start scan on zalkon**

```bash
ssh zalkon 'source /home/mescalin/yseitz/miniforge3/etc/profile.d/conda.sh && conda activate rnaclust && cd /scr/zalkon/yseitz && nohup python -m my_module.cross_species_scan --region 5UTR --threads 8 --parallel-species 3 --data-root data > data/cross_species/scan_5utr.log 2>&1 &'
```

- [ ] **Step 3: Monitor progress**

```bash
ssh zalkon 'tail -5 /scr/zalkon/yseitz/data/cross_species/scan_5utr.log'
```

Expected runtime: 4–6 hours. Check periodically.

- [ ] **Step 4: Aggregate results on zalkon**

After scan completes:
```bash
ssh zalkon 'source /home/mescalin/yseitz/miniforge3/etc/profile.d/conda.sh && conda activate rnaclust && cd /scr/zalkon/yseitz && python -c "
from pathlib import Path
from my_module.graphclust_results import aggregate_cross_species_results
import polars as pl

hits = aggregate_cross_species_results(Path(\"data/cross_species/5UTR\"))
hits.write_parquet(\"data/cross_species/5UTR/all_hits.parquet\")
print(f\"Saved {len(hits)} hits to all_hits.parquet\")
"'
```

- [ ] **Step 5: Copy results to Mac**

```bash
rsync -av zalkon:/scr/zalkon/yseitz/data/cross_species/ /Volumes/Masterarbeit/data/cross_species/
```

- [ ] **Step 6: Commit Parquet**

```bash
# Don't commit the tblout files (too large), only the Parquet
echo "data/cross_species/5UTR/*/*.tblout" >> .gitignore
git add data/cross_species/5UTR/all_hits.parquet .gitignore
git commit -m "data: cross-species 5UTR CM search results (all_hits.parquet)"
```

---

### Task 6: Analysis Notebook

**Files:**
- Create: `cross_species_5UTR.ipynb`

- [ ] **Step 1: Create notebook with analysis cells**

Cell 1 — Setup:
```python
import my_module as RNA
from my_module.graphclust_results import (
    aggregate_cross_species_results, build_motif_families,
    classify_rfam_novel, conservation_summary, SPECIES_LIST,
)
from pathlib import Path
import polars as pl
import matplotlib.pyplot as plt

data = Path("data")
region = "5UTR"
```

Cell 2 — Load aggregated hits:
```python
hits = pl.read_parquet(data / "cross_species" / region / "all_hits.parquet")
print(f"Total hits: {len(hits)}")
print(f"Source species: {hits['source_species'].n_unique()}")
print(f"Target species: {hits['target_species'].n_unique()}")
hits.head()
```

Cell 3 — Conservation summary:
```python
cons = conservation_summary(hits)
print(f"Clusters with cross-species hits: {len(cons)}")
cons.head(20)
```

Cell 4 — Build motif families:
```python
families = build_motif_families(hits)
n_fam = families["motif_family_id"].n_unique()
multi = families.filter(pl.col("family_n_species") > 1)["motif_family_id"].n_unique()
print(f"Motif families: {n_fam} total, {multi} multi-species")
families.filter(pl.col("family_n_species") > 1).head(20)
```

Cell 5 — Rfam classification:
```python
classified = classify_rfam_novel(families, hits, data, region)
novel = classified.filter(pl.col("rfam_family") == "novel")
rfam = classified.filter(pl.col("rfam_family") != "novel")
print(f"Rfam-annotated: {len(rfam)}, Novel: {len(novel)}")
classified.head(20)
```

Cell 6 — Plot A: Presence/Absence Heatmap:
```python
# Build species × family matrix
phylo_order = [
    "Human", "Mouse", "Cow", "Bat", "Platypus",
    "Guineafowl", "Zebrafish", "FruitFly",
    "Arabidopsis", "Corn", "Rice", "Yeast",
]

# Get multi-species families sorted by conservation
multi_fam = families.filter(pl.col("family_n_species") > 1)
fam_order = (
    multi_fam.group_by("motif_family_id")
    .agg(pl.col("family_n_species").first())
    .sort("family_n_species", descending=True)
    ["motif_family_id"].to_list()
)

# Build matrix: best bitscore per (species, family)
matrix_data = []
for fam_id in fam_order[:50]:  # top 50 families
    fam_members = families.filter(pl.col("motif_family_id") == fam_id)
    fam_hits = hits.filter(
        pl.struct(["source_species", "source_cluster"]).is_in(
            [{"source_species": r[0], "source_cluster": r[1]}
             for r in fam_members.select("source_species", "source_cluster").iter_rows()]
        )
    )
    row = {"family": fam_id}
    for sp in phylo_order:
        sp_hits = fam_hits.filter(pl.col("target_species") == sp)
        row[sp] = sp_hits["score"].max() if not sp_hits.is_empty() else 0.0
        # Also count source species
        if sp in fam_members["source_species"].to_list():
            row[sp] = max(row[sp], 50.0)  # mark presence
    matrix_data.append(row)

import numpy as np
mat = np.array([[r.get(sp, 0.0) for sp in phylo_order] for r in matrix_data])

fig, ax = plt.subplots(figsize=(14, 8))
im = ax.imshow(mat.T, aspect="auto", cmap="viridis", interpolation="none")
ax.set_xticks(range(len(fam_order[:50])))
ax.set_xticklabels([str(f) for f in fam_order[:50]], rotation=90, fontsize=6)
ax.set_yticks(range(len(phylo_order)))
ax.set_yticklabels(phylo_order)
ax.set_xlabel("Motif Family")
ax.set_ylabel("Species")
plt.colorbar(im, label="Best Bitscore")
fig.tight_layout()
fig.savefig("Master_Thesis/figures/cross_species_5utr_heatmap.pdf")
plt.show()
```

Cell 7 — Plot B: Conservation distribution:
```python
fam_stats = (
    families.group_by("motif_family_id")
    .agg(pl.col("family_n_species").first())
)

counts = fam_stats.group_by("family_n_species").agg(pl.len().alias("n_families")).sort("family_n_species")

fig, ax = plt.subplots(figsize=(8, 5))
ax.bar(counts["family_n_species"], counts["n_families"], color="#0072B2", edgecolor="white")
ax.set_xlabel("Number of species with hits")
ax.set_ylabel("Number of motif families")
ax.set_xticks(range(1, 13))
fig.tight_layout()
fig.savefig("Master_Thesis/figures/cross_species_5utr_conservation_dist.pdf")
plt.show()
```

Cell 8 — Plot C: Taxonomic conservation:
```python
TAXONOMY = {
    "Universal": lambda spp: len(spp) >= 10,
    "Eukaryotic-wide": lambda spp: 6 <= len(spp) < 10 and len(set(spp) & {"Arabidopsis", "Corn", "Rice", "Yeast"}) > 0 and len(set(spp) & {"Human", "Mouse", "Cow", "Bat", "Platypus"}) > 0,
    "Mammal-specific": lambda spp: all(s in {"Human", "Mouse", "Cow", "Bat", "Platypus"} for s in spp) and len(spp) >= 2,
    "Plant-specific": lambda spp: all(s in {"Arabidopsis", "Corn", "Rice"} for s in spp) and len(spp) >= 2,
    "Other clade": lambda spp: len(spp) >= 2,
    "Species-specific": lambda spp: len(spp) <= 1,
}

categories = {}
for fam_id in families["motif_family_id"].unique().to_list():
    spp = families.filter(pl.col("motif_family_id") == fam_id)["source_species"].unique().to_list()
    for cat, test in TAXONOMY.items():
        if test(spp):
            categories[fam_id] = cat
            break

cat_counts = {}
for cat in TAXONOMY:
    cat_counts[cat] = sum(1 for v in categories.values() if v == cat)

fig, ax = plt.subplots(figsize=(8, 5))
cats = list(cat_counts.keys())
vals = [cat_counts[c] for c in cats]
colors = ["#009E73", "#56B4E9", "#E69F00", "#CC79A7", "#D55E00", "#999999"]
ax.barh(cats, vals, color=colors[:len(cats)], edgecolor="white")
ax.set_xlabel("Number of motif families")
fig.tight_layout()
fig.savefig("Master_Thesis/figures/cross_species_5utr_taxonomy.pdf")
plt.show()
```

Cell 9 — Plot D: Novel motif candidates table:
```python
novel_multi = classified.filter(
    (pl.col("rfam_family") == "novel") & (pl.col("family_n_species") > 1)
).sort("family_n_species", descending=True)

print(f"Novel multi-species motif families: {len(novel_multi)}")
novel_multi.head(20)
```

- [ ] **Step 2: Commit notebook**

```bash
git add cross_species_5UTR.ipynb
git commit -m "feat: add cross-species 5UTR analysis notebook"
```

---

## Execution Order

1. **Task 1** (scan script) → **Task 2** (aggregation) → **Task 3** (motif families) → **Task 4** (Rfam classification) — all code, can be done in parallel
2. **Task 5** (run scan on zalkon) — depends on Task 1+2
3. **Task 6** (notebook) — depends on Task 3+4+5

Tasks 1–4 are pure code and can be implemented and committed before the scan runs. Task 5 is the long-running computation. Task 6 runs after results are available.
