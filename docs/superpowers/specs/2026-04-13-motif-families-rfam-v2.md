# Motif Families + Rfam Classification v2

**Date:** 2026-04-13
**Status:** Draft
**Replaces:** Initial build_motif_families() and classify_rfam_novel() implementations

## Problem Statement

Three bugs in the current cross-species analysis functions:

1. **build_motif_families()** merges clusters that hit the same target sequence in a third species, without verifying reciprocal source-locus overlap. This produces one giant 12-species component instead of biologically meaningful families.

2. **classify_rfam_novel()** compares SEQ-fragment IDs with ENST-transcript IDs and uses species-level rather than locus-level matching. Result: 0 novel / 25 annotated — an artefact.

3. **Off-by-one** in overlap calculations: `min(end) - max(start)` without `+1` for 0-based inclusive coordinates.

## Design

### Coordinate Normalisation (prerequisite)

All coordinates are normalised to **0-based inclusive** before any comparison.

| Source | Raw format | Normalisation |
|--------|-----------|---------------|
| `_parse_tblout()` (cmsearch) | 1-based inclusive | `-1` on both start and stop |
| `load_mRNA_cmscan()` (Rfam cmscan) | Already 0-based inclusive | None (existing code handles this) |
| Fragment keys in `data.names` | Verify empirically before implementation | See verification step |
| BED region offsets | 0-based start | Used as additive offset |

**Overlap formula (everywhere):** `overlap = min(endA, endB) - max(startA, startB) + 1`

**Verification step:** Before implementing, compare one known fragment key + BED offset against the FASTA sequence to confirm the coordinate convention.

### build_motif_families() — Reciprocal Source-Locus Validation

**Edge definition:** Two clusters (sp_A, cl_A) and (sp_B, cl_B) are connected iff:

1. cl_A has cmsearch hits in sp_B that overlap a **member locus of cl_B** in sp_B (overlap ≥ threshold)
2. cl_B has cmsearch hits in sp_A that overlap a **member locus of cl_A** in sp_A (overlap ≥ threshold)

Both conditions must hold. Connected components over validated pairwise edges form the families.

**Member loci** are derived from `load_partition()` + `load_data_names()`:
- `load_partition()` gives (fragment_key → cluster assignment)
- `load_data_names()` gives (fragment_key → seq_id, start, stop, transcript_id)
- Together: for each (species, cluster), the set of (transcript_id, start, stop) intervals that belong to that cluster

**Data structures needed:**

```
cluster_loci: Dict[(species, cluster)] → List[(transcript_id, start_0based, stop_0based)]
```

Built from `load_partition()` joined with `load_data_names()` for each species.

**Cross-species hits** are indexed as:

```
hits_index: Dict[(source_species, source_cluster, target_species)] → List[(target_name, hit_start_0based, hit_stop_0based)]
```

**Edge validation algorithm:**

```
for each pair (sp_A, cl_A), (sp_B, cl_B) where sp_A ≠ sp_B:
    # Direction A→B: cl_A's hits in sp_B overlap cl_B's member loci in sp_B
    hits_A_in_B = hits_index[(sp_A, cl_A, sp_B)]
    loci_B = cluster_loci[(sp_B, cl_B)]
    forward_ok = any(overlap(hit, locus) >= threshold for hit in hits_A_in_B for locus in loci_B)

    # Direction B→A: cl_B's hits in sp_A overlap cl_A's member loci in sp_A
    hits_B_in_A = hits_index[(sp_B, cl_B, sp_A)]
    loci_A = cluster_loci[(sp_A, cl_A)]
    reverse_ok = any(overlap(hit, locus) >= threshold for hit in hits_B_in_A for locus in loci_A)

    if forward_ok and reverse_ok:
        union(cl_A, cl_B)
```

**Note:** This is O(n² × m²) in the worst case but n (clusters) ≈ 960 and m (loci per cluster) ≈ 3-7, so it is tractable.

### classify_rfam_novel() — Locus-Level Rfam Overlap

For each motif family, check whether any **member locus** overlaps with an Rfam annotation at the same mRNA position.

**Pipeline (reuses existing rfam_overlap.py infrastructure):**

1. For each family member (species, cluster), get member loci from `cluster_loci`
2. Convert fragment-local coordinates to mRNA-absolute using `load_region_offsets()`
3. Load Rfam cmscan hits for that species via `load_mRNA_cmscan()`
4. Check overlap using the standard formula: `min(endA, endB) - max(startA, startB) + 1 >= overlap_nt`
5. If any member locus in any species overlaps an Rfam hit → annotate family with that Rfam name
6. Otherwise → "novel"

**Important:** The Rfam overlap is against **member loci** (the cluster's own sequences), not against cross-species cmsearch hits. `all_hits` is only used for statistics (total_hits, mean_score).

### Output

**build_motif_families()** returns:
- `source_species`, `source_cluster`, `motif_family_id`, `family_size`, `family_n_species`, `family_species_list`

**classify_rfam_novel()** returns:
- `motif_family_id`, `family_n_species`, `family_species_list`, `rfam_family`, `total_hits`, `mean_score`

## Files Modified

| File | Change |
|------|--------|
| `my_module/graphclust_results.py` | Replace `build_motif_families()` and `classify_rfam_novel()` |
| `my_module/graphclust_results.py` | Fix `_parse_tblout()` to normalise to 0-based inclusive |

## Dependencies

- Existing: `load_partition()`, `load_data_names()`, `load_center_ids()` from graphclust_results.py
- Existing: `load_region_offsets()`, `load_mRNA_cmscan()` from rfam_overlap.py
- Existing: `parse_cmscan_output()` from cmscan_output_parser.py
