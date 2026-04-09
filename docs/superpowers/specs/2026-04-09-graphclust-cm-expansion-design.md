# Design: GraphClust CM Expansion, CMfinder Default, R-scape Integration

**Date:** 2026-04-09  
**Status:** Approved  

---

## Overview

Three GraphClust2-inspired improvements to GraphClust v1:

1. **Iterative CM-Refinement (Stage 8b)** — per-cluster expansion loop after Stage 8: new cmsearch hits are added to the alignment, CMfinder refines the motif, a new CM is built, cmsearch runs again. Repeats until convergence (no new sequences) or max iterations.
2. **CMfinder always active** — `center_model_type` default changed from 1 to 5 in config.
3. **R-scape integration** — R-scape runs on the final `model.stk` after the expansion loop; covariation statistics written to `cluster.stats`.

---

## Architecture

```
Stage 7: alignCenter.pl        (CMfinder now default via center_model_type=5)
    ↓
Stage 8: gc_cmsearch.pl        (unchanged — first cmsearch run)
    ↓  cmsearch.DONE
Stage 8b: gc_cm_expand.pl      (NEW)
    ├─ Loop (max cm_expand_max_iter=3):
    │   1. Extract new hits from CMSEARCH/*.tabresult not yet in alignment
    │   2. If no new sequences → break (convergence)
    │   3. Re-align new seqs into existing alignment with mlocarna
    │   4. CMfinder on expanded set using model.tree.stk as seed
    │   5. cmbuild + cmsearch → new tabresult
    │   DONE-file: expand.iterN.DONE
    └─ After final loop:
        6. R-scape on final model.stk → PDF + .out in MODEL/
        7. Write expand.stats
    expand.DONE
    ↓
Stage 9: gc_results_cluster.pl (reads expand.stats → appends to cluster.stats)
```

---

## New File: `gc_cm_expand.pl`

### Interface

```
perl gc_cm_expand.pl
  --root-dir    <path>     # GraphClust root directory
  --cluster-dir <path>     # CLUSTER/$idx.cluster/
  --db          <path>     # FASTA/data.fasta.scan
  --max-iter    3          # from config: cm_expand_max_iter
  --cpu         N          # thread count
  --verbose                # optional
```

### Internal Logic

**State tracked across iterations:**
- `$current_tabresult` — starts as `CMSEARCH/model.stk.cm.tabresult` (Stage 8 output); updated to `CMSEARCH/expand.iterN/model.stk.cm.tabresult` after each iteration.
- `$known_seq_ids` — set of sequence IDs already in the alignment; grows each iteration.

**Per-iteration:**

1. **Hit extraction** — parse `$current_tabresult`, filter by `cm_min_bitscore` and `cm_bitscore_sig` (same thresholds as `glob_results.pl`). New sequences = hit IDs not in `$known_seq_ids`.
2. **Convergence check** — if no new sequences found, write `expand.DONE` and exit loop.
3. **Re-alignment** — extract new sequences from `data.fasta.scan`, merge with current `MODEL/model.tree.fa` into a combined FASTA, call `GraphClust::mlocarna_center()` (already in `GraphClust.pm`) on the full merged set (existing + new). Full re-alignment, not incremental.
4. **CMfinder** — `cmfinder --g 1.0 -a MODEL/model.tree.stk <expanded.fa> MODEL/model.cmfinder.stk.iter$N`. Identical call to `alignCenter.pl:285`. If CMfinder produces a non-empty result with negative energy, use it as new `model.stk`; otherwise keep existing.
5. **cmbuild + cmsearch** — identical commands to `gc_cmsearch.pl`, written inline. Output: `CMSEARCH/expand.iter$N/model.stk.cm.tabresult`.
6. **Update state** — `$current_tabresult` → new tabresult; `$known_seq_ids` += new IDs; `$total_added` += count new.
7. **DONE-file** — `expand.iter1.DONE`, `expand.iter2.DONE`, ... for restart safety.

**After loop:**

7. **R-scape** — `R-scape --outdir MODEL/ MODEL/model.stk` (identical call to `galaxy_wrappers/CollectResults/gc_res.pl:731`). Produces `model_1.R2R.sto.pdf` and `model.R-scape.out`.
8. **expand.stats** — written to `CLUSTER/$idx.cluster/expand.stats`:
   ```
   RSCAPE_PAIRS 4 RSCAPE_SIGNIF 2 RSCAPE_EVALUE 0.003 EXPAND_ITERS 2 EXPAND_SEQS_ADDED 7
   ```
   - `RSCAPE_PAIRS`: total base pairs tested
   - `RSCAPE_SIGNIF`: pairs with E-value < 0.05
   - `RSCAPE_EVALUE`: minimum E-value found (most significant pair)
   - `EXPAND_ITERS`: number of expansion iterations completed
   - `EXPAND_SEQS_ADDED`: total new sequences added across all iterations

---

## Changes to Existing Files

### `GraphClust_config.pm`

```perl
# Add:
PATH_RSCAPE        => "/Users/yannick/.local/share/mamba/envs/rnaclust/bin/",
cm_expand_max_iter => 3,

# Change:
center_model_type  => 5,    # was: 1
```

### `MASTER_GraphClust.pl`

New Stage 8b block after the Stage 8 block (~line 977), inside the `foreach keys %toDo_models` loop:

```perl
## stage 8b: CM expansion loop (iterative refinement + R-scape)
if (  -e "$CLUSTER_DIR/$clus_idx.cluster/cmsearch.DONE"
  && !-e "$CLUSTER_DIR/$clus_idx.cluster/expand.DONE" ) {

  my $job_name = "Round $CI cluster $clus_idx stage 8b";
  my $job_uuid = "stage8b-$clus_idx";

  my $CMD_expand = [];
  $CMD_expand->[0] = "perl $BIN_DIR/gc_cm_expand.pl";
  $CMD_expand->[1] = "--root-dir $in_ROOTDIR "
                   . "--cluster-dir $curr_cluster_dir "
                   . "--db $FASTA_DIR/data.fasta.scan "
                   . "--max-iter $CONFIG{cm_expand_max_iter} "
                   . "--cpu $stage8_cpu_threads ";
  $CMD_expand->[1] .= "--verbose " if ($in_verbose);

  my $sge_status = job_call(
    $job_name, "$BIN_DIR/gc_cmsearch.sge", $CMD_expand, 1,
    $SGE_ERR_DIR, $in_USE_SGE,
    "$curr_cluster_dir/SGE_log_expand",
    "$EVAL_DIR/times/time.stage.8b.$clus_idx",
    0, $stage8_qsub_opts, $NUM_THREADS, $job_uuid, undef, $stage8_local_slots
  );

  if ( $sge_status->[0] == 1 && $sge_status->[1] == 0 ) {
    system_call("touch $CLUSTER_DIR/$clus_idx.cluster/expand.DONE");
    $trigger_new_partition = 1;
  } elsif ( $sge_status->[1] == 1 ) {
    print "Round $CI cluster $clus_idx stage 8b: SGE job error! Skip...\n";
    $cluster_error++;
  }
}

next if ( !-e "$CLUSTER_DIR/$clus_idx.cluster/expand.DONE" );
```

Also add `expand.DONE` to the Stage 5 skip-check (blacklist logic, ~line 746):
```perl
if ( -e "$CLUSTER_DIR/$clus_idx.cluster/expand.DONE" ) { ... }
```

### `gc_results_cluster.pl`

After the existing RNAz block (~line 225), before `print STAT`:

```perl
## read R-scape and expansion stats if available
my %expand_stats = ( pairs => 0, signif => 0, evalue => '-', iters => 0, added => 0 );
## $clus_idx is the merged cluster index; expand.stats lives under the original cluster(s).
## @orig_clus contains the original cluster indices (could be >1 after merging).
## Use the first original cluster that has an expand.stats file.
my $expand_stats_file = "";
foreach my $orig_cl (@orig_clus) {
  my $f = "$in_root_dir/CLUSTER/$orig_cl.cluster/expand.stats";
  if ( -e $f ) { $expand_stats_file = $f; last; }
}
if ( $expand_stats_file && -e $expand_stats_file ) {
  open( my $EXP, $expand_stats_file );
  my $line = <$EXP>;
  close($EXP);
  chomp $line;
  $expand_stats{pairs}  = $1 if ( $line =~ /RSCAPE_PAIRS (\S+)/ );
  $expand_stats{signif} = $1 if ( $line =~ /RSCAPE_SIGNIF (\S+)/ );
  $expand_stats{evalue} = $1 if ( $line =~ /RSCAPE_EVALUE (\S+)/ );
  $expand_stats{iters}  = $1 if ( $line =~ /EXPAND_ITERS (\S+)/ );
  $expand_stats{added}  = $1 if ( $line =~ /EXPAND_SEQS_ADDED (\S+)/ );
}
```

Append to `print STAT` line:
```perl
print STAT "RSCAPE_PAIRS $expand_stats{pairs} RSCAPE_SIGNIF $expand_stats{signif} "
         . "RSCAPE_EVALUE $expand_stats{evalue} "
         . "EXPAND_ITERS $expand_stats{iters} EXPAND_SEQS_ADDED $expand_stats{added} ";
```

### `alignCenter.pl`

No functional change needed. The `center_model_type=1` fallback on line 64 (when CMfinder binary is missing) remains as a safety net. The new default of 5 in `GraphClust_config.pm` activates CMfinder automatically for all runs where the binary is present.

---

## File Layout After Implementation

```
CLUSTER/$idx.cluster/
├── MODEL/
│   ├── model.stk              (final, after expansion)
│   ├── model.tree.stk         (pre-expansion seed)
│   ├── model.cmfinder.stk     (CMfinder output per iteration)
│   ├── model_1.R2R.sto.pdf    (R-scape covariation plot)
│   └── model.R-scape.out      (R-scape statistics)
├── CMSEARCH/
│   ├── model.stk.cm.tabresult     (Stage 8 initial)
│   ├── expand.iter1/
│   │   └── model.stk.cm.tabresult (Stage 8b iteration 1)
│   └── expand.iter2/
│       └── model.stk.cm.tabresult (Stage 8b iteration 2)
├── cmsearch.DONE
├── expand.iter1.DONE
├── expand.DONE
└── expand.stats
```

---

## Config Keys Added

| Key | Default | Description |
|-----|---------|-------------|
| `PATH_RSCAPE` | `rnaclust/bin/` | Path to R-scape binary |
| `cm_expand_max_iter` | `3` | Max expansion iterations per cluster |
| `center_model_type` | `5` (changed from `1`) | CMfinder active by default |

---

## CMfinder04 Bundling

CMfinder 0.4.1.18 (GPL3) is bundled as `GraphClust/cmfinder_src/` with the following minimal subset:

```
cmfinder_src/
├── cmfinder04/        (C source, 540KB — the actual CMfinder code)
├── lib-infernal-1.1/  (bundled Infernal, 15MB — no dev headers in conda env)
├── lib-vienna-1.4/    (bundled ViennaRNA 1.4, 700KB — system ViennaRNA 2.x has incompatible API)
├── configure          (pre-generated autoconf script)
├── configure.ac / Makefile.am / Makefile.in / aclocal.m4
├── config.guess / config.sub / depcomp / install-sh / missing
└── build.sh           (wrapper script)
```

**Build command (on cluster, Linux x86_64):**
```bash
bash cmfinder_src/build.sh
# or from GraphClust root: make cmfinder
```

Binary output: `cmfinder_src/bin/cmfinder`

After building, set in the run-specific `config` file (not in `GraphClust_config.pm`):
```
PATH_CMFINDER = /path/to/GraphClust/cmfinder_src/bin/
```

**Why not system libraries:** The rnaclust conda env has ViennaRNA 2.7.2 (API incompatible with CMfinder04's v1.4 calls) and Infernal executables only (no headers/static libs). Bundling is necessary.

**Apple Silicon:** CMfinder04 requires SSE/VMX SIMD (bundled Infernal 1.1). ARM64 is not supported. `build.sh` exits with a clear error; `alignCenter.pl` falls back to `center_model_type=1` automatically.

---

## Constraints and Assumptions

- R-scape must be installed in the `rnaclust` conda environment.
- CMfinder binary must be present at `PATH_CMFINDER`; if missing, `alignCenter.pl` falls back to `center_model_type=1` (existing safety net unchanged).
- The expansion loop re-uses `gc_cmsearch.sge` as SGE template (same resource requirements as Stage 8).
- Bitscore filtering in hit extraction uses the same thresholds as `glob_results.pl` (`cm_min_bitscore`, `cm_bitscore_sig`) for consistency.
- Maximum 3 expansion iterations is a hard cap to prevent runaway runtimes; configurable via `cm_expand_max_iter`.
