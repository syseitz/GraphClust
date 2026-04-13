# Motif Families + Rfam Classification v2 — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace `build_motif_families()` with reciprocal source-locus validation and `classify_rfam_novel()` with proper coordinate-mapped Rfam overlap, fixing three bugs (transitive merging, ID mismatch, off-by-one).

**Architecture:** (1) Build a cluster-loci index from `load_partition()` + `data.names` with mRNA-absolute coordinates, (2) validate motif family edges by checking reciprocal locus overlap, (3) classify families by Rfam overlap at the locus level using `load_mRNA_cmscan()`.

**Tech Stack:** Python 3.10, Polars, existing `rfam_overlap.py` helpers.

---

## Coordinate Convention (verified empirically)

- `data.names` `orig_id`: `ENST...:START-STOP_IDX` where START-STOP are **1-based inclusive** positions within the region (5UTR/3UTR) on the mRNA
- For 5UTR: region offset = 0, so these ARE mRNA coordinates
- For 3UTR: region offset = BED 3UTR start (from `load_region_offsets()`)
- `_parse_tblout()` cmsearch hits: `target_from`/`target_to` are **1-based inclusive** positions within the fragment
- `load_mRNA_cmscan()` Rfam hits: already normalised to **0-based inclusive** mRNA coordinates
- **Normalisation target:** 0-based inclusive everywhere. Formula: `overlap = min(endA, endB) - max(startA, startB) + 1`

## File Map

| File | Action | Responsibility |
|------|--------|---------------|
| `my_module/graphclust_results.py` | Modify | Replace `build_motif_families()`, `classify_rfam_novel()`, fix `_parse_tblout()` |

No new files. All changes in one file.

---

### Task 1: Fix `_parse_tblout()` coordinate normalisation

**Files:**
- Modify: `my_module/graphclust_results.py:434-452`

- [ ] **Step 1: Normalise cmsearch coordinates to 0-based inclusive**

In `_parse_tblout()`, change the coordinate parsing from raw 1-based to 0-based inclusive:

```python
def _parse_tblout(tblout_path: Path) -> List[Dict]:
    """Parse Infernal tblout format.

    Coordinates are normalised to 0-based inclusive.
    """
    hits = []
    with open(tblout_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 16:
                continue
            raw_from = int(fields[7])
            raw_to = int(fields[8])
            hits.append({
                "target_name": fields[0],
                "target_from": min(raw_from, raw_to) - 1,  # 1-based → 0-based
                "target_to": max(raw_from, raw_to) - 1,    # 1-based → 0-based
                "strand": fields[9],
                "score": float(fields[14]),
                "evalue": float(fields[15]),
            })
    return hits
```

Note: `min/max` handles reverse-strand hits where from > to.

- [ ] **Step 2: Commit**

```bash
git add my_module/graphclust_results.py
git commit -m "fix: normalise _parse_tblout to 0-based inclusive coordinates"
```

---

### Task 2: Add `_build_cluster_loci()` helper

**Files:**
- Modify: `my_module/graphclust_results.py` (add new function before `build_motif_families`)

This helper builds an index of all member loci for each (species, cluster) pair, in mRNA-absolute 0-based inclusive coordinates.

- [ ] **Step 1: Add the helper function**

```python
def _build_cluster_loci(
    data_root: Path,
    species_list: List[str],
    region: str,
) -> Dict[Tuple[str, str], List[Tuple[str, int, int]]]:
    """Build an index of member loci per (species, cluster).

    For each cluster, collects all member fragment loci from the partition
    and maps them to mRNA-absolute 0-based inclusive coordinates.

    Args:
        data_root: Root data directory.
        species_list: Species to index.
        region: Region (5UTR or 3UTR).

    Returns:
        Dict mapping (species, cluster_name) to a list of
        (transcript_id, mRNA_start_0, mRNA_stop_0) tuples.
    """
    try:
        from .rfam_overlap import load_region_offsets
    except ImportError:
        from my_module.rfam_overlap import load_region_offsets

    cluster_loci: Dict[Tuple[str, str], List[Tuple[str, int, int]]] = {}

    for species in tqdm(species_list, desc="Building cluster loci"):
        gc_dir = _find_graphclust_dir(data_root, species, region)
        if gc_dir is None:
            continue

        try:
            partition = load_partition(gc_dir)
            names = load_data_names(gc_dir)
        except Exception:
            continue

        # Load region offsets for coordinate transformation
        try:
            offsets = load_region_offsets(species, base_dir=str(data_root))
        except Exception:
            offsets = {}

        # Join partition with names to get orig_id per hit
        joined = partition.join(
            names.select("seq_id", "orig_id", "transcript_id"),
            on="seq_id",
            how="left",
        )

        for row in joined.iter_rows(named=True):
            cluster_name = row["cluster_name"]
            orig_id = row.get("orig_id", "")
            transcript_id = row.get("transcript_id", "")

            if not orig_id or ":" not in orig_id:
                continue

            # Parse orig_id: ENST...:START-STOP_IDX
            coord_part = orig_id.split(":")[1]  # START-STOP_IDX
            coord_part = re.sub(r"_\d+$", "", coord_part)  # remove _IDX suffix
            parts = coord_part.split("-")
            if len(parts) != 2:
                continue

            try:
                region_start_1 = int(parts[0])  # 1-based
                region_stop_1 = int(parts[1])    # 1-based
            except ValueError:
                continue

            # Get region offset from BED (for 3UTR; for 5UTR offset=0)
            offset = 0
            if transcript_id in offsets and region in offsets[transcript_id]:
                offset = offsets[transcript_id][region]

            # Convert to mRNA-absolute 0-based inclusive
            mRNA_start = offset + region_start_1 - 1
            mRNA_stop = offset + region_stop_1 - 1

            key = (species, cluster_name)
            cluster_loci.setdefault(key, []).append(
                (transcript_id, mRNA_start, mRNA_stop)
            )

    print(f"Built loci index: {len(cluster_loci)} (species, cluster) pairs")
    return cluster_loci
```

- [ ] **Step 2: Quick test**

```bash
cd /Volumes/Masterarbeit
python -c "
from my_module.graphclust_results import _build_cluster_loci
from pathlib import Path
loci = _build_cluster_loci(Path('data'), ['Mouse'], '5UTR')
print(f'Mouse clusters: {len(loci)}')
sample_key = next(iter(loci))
print(f'Sample: {sample_key} → {loci[sample_key][:3]}')
"
```

- [ ] **Step 3: Commit**

```bash
git add my_module/graphclust_results.py
git commit -m "feat: add _build_cluster_loci helper for mRNA coordinate mapping"
```

---

### Task 3: Replace `build_motif_families()` with reciprocal source-locus validation

**Files:**
- Modify: `my_module/graphclust_results.py` (replace existing function)

- [ ] **Step 1: Replace the function**

```python
def build_motif_families(
    all_hits: pl.DataFrame,
    data_root: Path,
    region: str = "5UTR",
    overlap_nt: int = 20,
) -> pl.DataFrame:
    """Group cross-species clusters into motif families via reciprocal source-locus overlap.

    Two clusters (sp_A, cl_A) and (sp_B, cl_B) are connected iff:
    1. cl_A has cmsearch hits in sp_B that overlap a member locus of cl_B
    2. cl_B has cmsearch hits in sp_A that overlap a member locus of cl_A

    Connected components over validated pairwise edges form families.

    Args:
        all_hits: Aggregated cross-species hits (0-based inclusive coords).
        data_root: Root data directory.
        region: Region (5UTR or 3UTR).
        overlap_nt: Minimum overlap in nt for locus match.

    Returns:
        DataFrame with columns: source_species, source_cluster,
        motif_family_id, family_size, family_n_species, family_species_list.
    """
    if all_hits.is_empty():
        return pl.DataFrame()

    # Collect species from hits
    species_in_hits = sorted(
        set(all_hits["source_species"].unique().to_list())
        | set(all_hits["target_species"].unique().to_list())
    )

    # Build cluster loci index (member loci in mRNA coordinates)
    cluster_loci = _build_cluster_loci(data_root, species_in_hits, region)

    # Build cross-species hits index: (source_sp, source_cl, target_sp) → list of (target_name, start, stop)
    # Also build a mapping: target_name → (transcript_id, fragment_mRNA_start)
    # by loading data.names for each target species
    seq_to_transcript: Dict[Tuple[str, str], Tuple[str, int]] = {}
    for sp in species_in_hits:
        gc_dir = _find_graphclust_dir(data_root, sp, region)
        if gc_dir is None:
            continue
        try:
            names = load_data_names(gc_dir)
        except Exception:
            continue

        try:
            from .rfam_overlap import load_region_offsets
        except ImportError:
            from my_module.rfam_overlap import load_region_offsets
        try:
            offsets = load_region_offsets(sp, base_dir=str(data_root))
        except Exception:
            offsets = {}

        for row in names.iter_rows(named=True):
            orig_id = row.get("orig_id", "")
            tid = row.get("transcript_id", "")
            if not orig_id or ":" not in orig_id:
                continue
            coord_part = re.sub(r"_\d+$", "", orig_id.split(":")[1])
            parts = coord_part.split("-")
            if len(parts) != 2:
                continue
            try:
                frag_start_1 = int(parts[0])
            except ValueError:
                continue

            offset = 0
            if tid in offsets and region in offsets[tid]:
                offset = offsets[tid][region]

            # Store: (species, seq_id) → (transcript_id, fragment_mRNA_start_0based)
            seq_to_transcript[(sp, row["seq_id"])] = (tid, offset + frag_start_1 - 1)

    # Index hits: convert fragment-local cmsearch coords to mRNA-absolute
    # hits are already 0-based inclusive (after Task 1 fix)
    hits_by_direction: Dict[Tuple[str, str, str], List[Tuple[str, int, int]]] = {}
    for row in all_hits.iter_rows(named=True):
        src_key = (row["source_species"], row["source_cluster"])
        tgt_sp = row["target_species"]
        tgt_name = row["target_name"]  # e.g. SEQ532

        mapping = seq_to_transcript.get((tgt_sp, tgt_name))
        if mapping is None:
            continue
        tid, frag_mRNA_start = mapping

        # Convert fragment-local hit to mRNA-absolute
        hit_mRNA_start = frag_mRNA_start + row["target_from"]
        hit_mRNA_stop = frag_mRNA_start + row["target_to"]

        direction_key = (row["source_species"], row["source_cluster"], tgt_sp)
        hits_by_direction.setdefault(direction_key, []).append(
            (tid, hit_mRNA_start, hit_mRNA_stop)
        )

    # Collect all source nodes
    all_nodes = sorted(set(
        (row["source_species"], row["source_cluster"])
        for row in all_hits.select("source_species", "source_cluster").unique().iter_rows()
    ))
    node_idx = {n: i for i, n in enumerate(all_nodes)}

    # Union-Find
    parent = list(range(len(all_nodes)))

    def _find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def _union(x: int, y: int) -> None:
        rx, ry = _find(x), _find(y)
        if rx != ry:
            parent[rx] = ry

    def _loci_overlap(
        hits: List[Tuple[str, int, int]],
        loci: List[Tuple[str, int, int]],
    ) -> bool:
        """Check if any hit overlaps any locus by >= overlap_nt."""
        for h_tid, h_start, h_stop in hits:
            for l_tid, l_start, l_stop in loci:
                if h_tid != l_tid:
                    continue
                ov = min(h_stop, l_stop) - max(h_start, l_start) + 1
                if ov >= overlap_nt:
                    return True
        return False

    # Validate pairwise edges
    print("Validating reciprocal source-locus edges...")
    edges_checked = 0
    edges_valid = 0
    for i, (sp_a, cl_a) in enumerate(all_nodes):
        loci_a = cluster_loci.get((sp_a, cl_a), [])
        if not loci_a:
            continue
        for j in range(i + 1, len(all_nodes)):
            sp_b, cl_b = all_nodes[j]
            if sp_a == sp_b:
                continue
            loci_b = cluster_loci.get((sp_b, cl_b), [])
            if not loci_b:
                continue

            # Forward: cl_A hits in sp_B overlap cl_B's loci in sp_B
            hits_a_in_b = hits_by_direction.get((sp_a, cl_a, sp_b), [])
            if not hits_a_in_b:
                continue
            forward = _loci_overlap(hits_a_in_b, loci_b)

            if not forward:
                continue

            # Reverse: cl_B hits in sp_A overlap cl_A's loci in sp_A
            hits_b_in_a = hits_by_direction.get((sp_b, cl_b, sp_a), [])
            if not hits_b_in_a:
                continue
            reverse = _loci_overlap(hits_b_in_a, loci_a)

            edges_checked += 1
            if forward and reverse:
                _union(i, j)
                edges_valid += 1

    # Assign family IDs
    root_members: Dict[int, List[int]] = {}
    for i in range(len(all_nodes)):
        root = _find(i)
        root_members.setdefault(root, []).append(i)

    sorted_families = sorted(root_members.values(), key=len, reverse=True)
    rows: List[Dict] = []
    for fid, members in enumerate(sorted_families):
        species_in_fam = sorted({all_nodes[i][0] for i in members})
        for i in members:
            sp, cl = all_nodes[i]
            rows.append({
                "source_species": sp,
                "source_cluster": cl,
                "motif_family_id": fid,
                "family_size": len(members),
                "family_n_species": len(species_in_fam),
                "family_species_list": species_in_fam,
            })

    result = pl.DataFrame(rows)
    n_fam = result["motif_family_id"].n_unique()
    multi = result.filter(pl.col("family_n_species") > 1)["motif_family_id"].n_unique()
    print(f"Edges: {edges_checked} checked, {edges_valid} valid")
    print(f"Built {n_fam} motif families ({multi} multi-species)")
    return result.sort(["motif_family_id", "source_species", "source_cluster"])
```

**Important:** The function signature now requires `data_root` and `region` parameters. Update the `__init__.py` export and notebook accordingly.

- [ ] **Step 2: Update `__init__.py`**

No change needed — `build_motif_families` is already exported. But the notebook call must be updated to pass `data_root` and `region`.

- [ ] **Step 3: Commit**

```bash
git add my_module/graphclust_results.py
git commit -m "feat: replace build_motif_families with reciprocal source-locus validation"
```

---

### Task 4: Replace `classify_rfam_novel()` with locus-level Rfam overlap

**Files:**
- Modify: `my_module/graphclust_results.py` (replace existing function)

- [ ] **Step 1: Replace the function**

```python
def classify_rfam_novel(
    motif_families: pl.DataFrame,
    all_hits: pl.DataFrame,
    data_root: Path,
    region: str = "5UTR",
    overlap_nt: int = 10,
) -> pl.DataFrame:
    """Classify motif families as Rfam-annotated or novel via locus-level overlap.

    For each family, checks whether any member locus overlaps an Rfam cmscan
    annotation on the same mRNA transcript at the same coordinates.

    Args:
        motif_families: Output of ``build_motif_families()``.
        all_hits: Aggregated cross-species hits (for statistics only).
        data_root: Root data directory.
        region: Region (5UTR or 3UTR).
        overlap_nt: Minimum overlap in nt for Rfam match.

    Returns:
        DataFrame: motif_family_id, family_n_species, family_species_list,
        rfam_family, total_hits, mean_score.
    """
    try:
        from .rfam_overlap import load_region_offsets, load_mRNA_cmscan
    except ImportError:
        from my_module.rfam_overlap import load_region_offsets, load_mRNA_cmscan

    if motif_families.is_empty():
        return pl.DataFrame()

    # Step 1: Per-family hit statistics from all_hits
    hits_stats: Dict[int, Dict] = {}
    if not all_hits.is_empty():
        joined = (
            all_hits
            .join(
                motif_families.select("source_species", "source_cluster", "motif_family_id"),
                on=["source_species", "source_cluster"],
                how="left",
            )
            .drop_nulls("motif_family_id")
        )
        for row in (
            joined.group_by("motif_family_id")
            .agg([pl.len().alias("total_hits"), pl.col("score").mean().alias("mean_score")])
            .iter_rows(named=True)
        ):
            hits_stats[row["motif_family_id"]] = {
                "total_hits": row["total_hits"],
                "mean_score": row["mean_score"],
            }

    # Step 2: Build member loci (reuse _build_cluster_loci)
    species_list = motif_families["source_species"].unique().to_list()
    cluster_loci = _build_cluster_loci(data_root, species_list, region)

    # Step 3: Load Rfam cmscan per species
    rfam_by_species: Dict[str, List[Tuple[str, int, int, str]]] = {}
    for sp in species_list:
        try:
            df = load_mRNA_cmscan(sp, base_dir=str(data_root))
            entries = []
            for _, r in df.iterrows():
                entries.append((
                    str(r["transcript_id"]),
                    int(r["cmscan_start_0"]),
                    int(r["cmscan_end_0"]),
                    str(r["target_name"]),
                ))
            rfam_by_species[sp] = entries
        except Exception:
            continue

    # Step 4: For each family, check member loci vs Rfam
    from collections import Counter

    unique_fam: Dict[int, Dict] = {}
    for row in motif_families.iter_rows(named=True):
        fid = row["motif_family_id"]
        if fid not in unique_fam:
            unique_fam[fid] = {
                "family_n_species": row["family_n_species"],
                "family_species_list": row["family_species_list"],
                "members": [],
            }
        unique_fam[fid]["members"].append((row["source_species"], row["source_cluster"]))

    output_rows: List[Dict] = []
    for fid, finfo in sorted(unique_fam.items()):
        rfam_matches: List[str] = []

        for sp, cl in finfo["members"]:
            loci = cluster_loci.get((sp, cl), [])
            rfam_entries = rfam_by_species.get(sp, [])
            if not loci or not rfam_entries:
                continue

            for l_tid, l_start, l_stop in loci:
                for r_tid, r_start, r_stop, r_name in rfam_entries:
                    if l_tid != r_tid:
                        continue
                    ov = min(l_stop, r_stop) - max(l_start, r_start) + 1
                    if ov >= overlap_nt:
                        rfam_matches.append(r_name)

        rfam_family = "novel"
        if rfam_matches:
            rfam_family = Counter(rfam_matches).most_common(1)[0][0]

        stats = hits_stats.get(fid, {"total_hits": 0, "mean_score": 0.0})
        output_rows.append({
            "motif_family_id": fid,
            "family_n_species": finfo["family_n_species"],
            "family_species_list": finfo["family_species_list"],
            "rfam_family": rfam_family,
            "total_hits": stats["total_hits"],
            "mean_score": stats["mean_score"],
        })

    result = pl.DataFrame(output_rows).sort("family_n_species", descending=True)
    n_novel = (result["rfam_family"] == "novel").sum()
    n_annot = len(result) - n_novel
    print(f"Classified {len(result)} families: {n_annot} Rfam-annotated, {n_novel} novel")
    return result
```

- [ ] **Step 2: Commit**

```bash
git add my_module/graphclust_results.py
git commit -m "feat: replace classify_rfam_novel with locus-level Rfam overlap"
```

---

### Task 5: Update notebook and run

**Files:**
- Modify: `cross_species_5UTR.ipynb`

- [ ] **Step 1: Update the motif families cell**

The `build_motif_families()` call now needs `data_root` and `region`:

```python
families_df = build_motif_families(hits_df, DATA_ROOT, region="5UTR")
```

- [ ] **Step 2: Re-run notebook**

```bash
jupyter nbconvert --to notebook --execute cross_species_5UTR.ipynb \
  --output cross_species_5UTR.ipynb \
  --ExecutePreprocessor.timeout=600 \
  --ExecutePreprocessor.kernel_name=python3
```

- [ ] **Step 3: Commit**

```bash
git add cross_species_5UTR.ipynb my_module/graphclust_results.py
git commit -m "feat: run cross-species analysis with v2 motif families + Rfam classification"
```

---

## Execution Order

1. **Task 1** — coordinate fix (prerequisite, no dependencies)
2. **Task 2** — cluster loci helper (depends on Task 1 for consistent coords)
3. **Task 3** — motif families (depends on Task 2)
4. **Task 4** — Rfam classification (depends on Task 2)
5. **Task 5** — notebook update (depends on Tasks 3+4)

Tasks 3 and 4 are independent of each other but both depend on Task 2.
