# GraphClust CM Expansion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add iterative CM-refinement (Stage 8b), R-scape integration, and CMfinder-default activation to GraphClust v1, following GraphClust2's approach.

**Architecture:** New `gc_cm_expand.pl` script handles the expansion loop (hit extraction → mlocarna → CMfinder → cmbuild → cmsearch, repeated until convergence). `MASTER_GraphClust.pl` calls it after Stage 8. `gc_results_cluster.pl` reads expansion stats into `cluster.stats`.

**Tech Stack:** Perl 5, Infernal (cmbuild/cmsearch), LocARNA (mlocarna), CMfinder04, R-scape, ViennaRNA (RNAalifold)

---

## File Map

| File | Action | Responsibility |
|------|--------|----------------|
| `gc_cm_expand.pl` | **Create** | Expansion loop: hit extraction, re-alignment, CMfinder, cmbuild+cmsearch, R-scape, stats |
| `gc_cm_expand.sge` | **Create** | SGE wrapper for gc_cm_expand.pl (same pattern as gc_cmsearch.sge) |
| `motif_plot.py` | **Create** | Visualization of cluster motifs (horizontal barplot, adapted from Galaxy GraphClust2) |
| `MASTER_GraphClust.pl` | **Modify** | Add Stage 8b block, update skip-check at line 746 |
| `gc_results_cluster.pl` | **Modify** | Read `expand.stats`, write R-scape fields to `cluster.stats` |
| `GraphClust_config.pm` | **Already done** | PATH_RSCAPE, cm_expand_max_iter already added |
| `Makefile.am` | **Already done** | cmfinder build target already added |
| `cmfinder_src/` | **Already done** | CMfinder source bundled |

---

### Task 1: Create `gc_cm_expand.sge` (SGE Wrapper)

**Files:**
- Create: `gc_cm_expand.sge`

- [ ] **Step 1: Read existing SGE wrapper for reference**

Read `gc_cmsearch.sge` to understand the wrapper pattern.

- [ ] **Step 2: Create `gc_cm_expand.sge`**

```bash
#!/bin/bash
#$ -S /bin/bash
#$ -e $3
#$ -o $3

## SGE wrapper for gc_cm_expand.pl
## $1 = perl command
## $2 = arguments
## $3 = log dir (used by -e and -o above)

$1 $2
```

This is identical in structure to `gc_cmsearch.sge`. The `$1` is `perl /path/to/gc_cm_expand.pl` and `$2` is the arguments string, both set by MASTER_GraphClust.pl via the `job_call()` function.

- [ ] **Step 3: Make executable**

Run: `chmod +x gc_cm_expand.sge`

- [ ] **Step 4: Commit**

```bash
git add gc_cm_expand.sge
git commit -m "feat: add SGE wrapper for gc_cm_expand.pl (stage 8b)"
```

---

### Task 2: Create `gc_cm_expand.pl` — Scaffolding and Hit Extraction

**Files:**
- Create: `gc_cm_expand.pl`

This task creates the script with argument parsing, config loading, and the hit-extraction logic. The expansion loop and R-scape come in the next tasks.

- [ ] **Step 1: Write the script header, argument parsing, and config loading**

```perl
#!/usr/bin/perl -w
use warnings;
use strict;

## GraphClust Stage 8b: Iterative CM expansion and R-scape
## After Stage 8 (cmsearch), this script:
## 1. Extracts new cmsearch hits not in the original model alignment
## 2. Re-aligns the expanded set with mlocarna
## 3. Refines with CMfinder (if available)
## 4. Rebuilds the CM and re-scans
## 5. Repeats until convergence or max iterations
## 6. Runs R-scape on the final model alignment

use List::Util qw/ min max /;
use FindBin;
use lib "$FindBin::Bin";
use GraphClust;

use vars qw(%CONFIG);
*CONFIG = \%GraphClust::CONFIG;

use Getopt::Long;
use Cwd qw(abs_path getcwd);
use File::Path qw(make_path);

my $binDir = "$FindBin::Bin";

my $in_root_dir;
my $in_cluster_dir;
my $db;
my $max_iter = 3;
my $verbose;
my $cpu_threads = 1;

usage()
  unless GetOptions(
  "root-dir=s"    => \$in_root_dir,
  "cluster-dir=s" => \$in_cluster_dir,
  "db=s"          => \$db,
  "max-iter=i"    => \$max_iter,
  "cpu=i"         => \$cpu_threads,
  "verbose"       => \$verbose,
  );

sub usage { die "Usage: gc_cm_expand.pl --root-dir DIR --cluster-dir DIR --db FILE [--max-iter N] [--cpu N] [--verbose]\n"; }

die "No --root-dir specified\n" unless $in_root_dir;
die "No --cluster-dir specified\n" unless $in_cluster_dir;
die "No --db specified\n" unless $db && -e $db;

################################################################################
## Load options from root-dir/config file
%GraphClust::CONFIG = readConfigFile( "$in_root_dir/config", 1 );
SUBSECTION("options loaded for gc_cm_expand.pl");
printConfig( \%CONFIG ) if ($verbose);

my $infernal_path  = $CONFIG{PATH_INFERNAL};
my $cmfinder_path  = $CONFIG{PATH_CMFINDER};
my $rscape_path    = $CONFIG{PATH_RSCAPE};
my $vrna_path      = $CONFIG{PATH_VRNA};
my $tmp_path       = $CONFIG{PATH_TMP};
my $cm_min_bitscore = $CONFIG{cm_min_bitscore};
my $cm_bitscore_sig = $CONFIG{cm_bitscore_sig};
my $cm_max_eval     = $CONFIG{cm_max_eval};
my $cm_top_only     = $CONFIG{cm_top_only};

my $model_dir    = "$in_cluster_dir/MODEL";
my $cmsearch_dir = "$in_cluster_dir/CMSEARCH";

## Infernal 1.1 compatibility (same logic as gc_cmsearch.pl)
my $cmsearch_help = `$infernal_path/cmsearch -h 2>&1`;
my $cmsearch_filter_opt = "--nohmm";
$cmsearch_filter_opt = "--fil-no-hmm" if ( $cmsearch_help =~ /--fil-no-hmm/ );
my $cmsearch_table_opt = "--tblout";
$cmsearch_table_opt = "--tabfile" if ( $cmsearch_help !~ /--tblout/ && $cmsearch_help =~ /--tabfile/ );
my $cmsearch_noalign_opt = "--noali";
$cmsearch_noalign_opt = "--noalign" if ( $cmsearch_help !~ /--noali/ && $cmsearch_help =~ /--noalign/ );
my $cmsearch_cpu_opt = "";
$cmsearch_cpu_opt = " --cpu $cpu_threads " if ( $cmsearch_help =~ /--cpu/ );

## Check CMfinder availability
my $use_cmfinder = 0;
if ( $cmfinder_path && $cmfinder_path ne "false" ) {
  if ( -e "${cmfinder_path}cmfinder" || -e "$cmfinder_path/cmfinder" ) {
    $use_cmfinder = 1;
    print "CMfinder available at $cmfinder_path\n";
  }
}
print "CMfinder: " . ( $use_cmfinder ? "enabled" : "disabled (binary not found)" ) . "\n";

################################################################################
## Helper: extract sequence IDs from a FASTA file

sub get_fasta_ids {
  my $fa_file = shift;
  my @fa = GraphClust::read_fasta_file($fa_file);
  return { map { $_ => 1 } @{ $fa[1] } };
}

################################################################################
## Helper: extract new hits from a tabresult file, excluding known IDs

sub extract_new_hits {
  my ( $tabresult_file, $known_ids_href ) = @_;

  my $cm_hits = GraphClust::read_CM_tabfile_ext(
    $tabresult_file, $cm_min_bitscore, $cm_max_eval, $cm_bitscore_sig, "expand"
  );

  return [] if ( !$cm_hits || !@{$cm_hits} );

  my @new_hits = ();
  my %seen_seqids = ();
  foreach my $hit ( @{$cm_hits} ) {
    next if exists $known_ids_href->{ $hit->{SEQID} };
    next if exists $seen_seqids{ $hit->{SEQID} };
    $seen_seqids{ $hit->{SEQID} } = 1;
    push( @new_hits, $hit );
  }

  return \@new_hits;
}

################################################################################
## Main logic starts below (continued in Task 3)

print "gc_cm_expand.pl: cluster-dir=$in_cluster_dir max-iter=$max_iter\n";
```

- [ ] **Step 2: Verify script loads without errors**

Run: `perl -c gc_cm_expand.pl`
Expected: `gc_cm_expand.pl syntax OK`

- [ ] **Step 3: Commit**

```bash
git add gc_cm_expand.pl
git commit -m "feat: gc_cm_expand.pl scaffolding — arg parsing, config, hit extraction"
```

---

### Task 3: `gc_cm_expand.pl` — Expansion Loop

**Files:**
- Modify: `gc_cm_expand.pl` (append after "Main logic starts below")

- [ ] **Step 1: Implement the expansion loop**

Append the following to `gc_cm_expand.pl`, replacing the placeholder comment:

```perl
################################################################################
## Expansion loop

## State: track which sequences are already in the alignment
my $model_fa_file = "$model_dir/model.tree.fa";
die "No model FASTA: $model_fa_file\n" unless -e $model_fa_file;

my $known_ids = get_fasta_ids($model_fa_file);
print "Initial model sequences: " . scalar( keys %{$known_ids} ) . "\n";

## Find the initial tabresult from Stage 8
my @initial_tabs = glob("$cmsearch_dir/*.tabresult");
die "No tabresult found in $cmsearch_dir\n" unless @initial_tabs;
my $current_tabresult = $initial_tabs[0];
print "Initial tabresult: $current_tabresult\n";

## Read the scan FASTA (all input sequences) for extracting new seqs
my @fa_scan = GraphClust::read_fasta_file("$db");

my $total_added = 0;
my $iterations_done = 0;

for ( my $iter = 1; $iter <= $max_iter; $iter++ ) {

  ## Skip if this iteration was already completed (restart safety)
  if ( -e "$in_cluster_dir/expand.iter$iter.DONE" ) {
    print "Iteration $iter already done, skipping.\n";

    ## Update state for next iteration
    my $iter_dir = "$cmsearch_dir/expand.iter$iter";
    my @iter_tabs = glob("$iter_dir/*.tabresult");
    $current_tabresult = $iter_tabs[0] if @iter_tabs;

    ## Re-read the model FASTA that was expanded in this iteration
    $known_ids = get_fasta_ids($model_fa_file) if -e $model_fa_file;
    $iterations_done = $iter;
    next;
  }

  SUBSECTION("Expansion iteration $iter / $max_iter");

  ## Step 1: Extract new hits
  my $new_hits = extract_new_hits( $current_tabresult, $known_ids );
  my $num_new = scalar( @{$new_hits} );
  print "Iteration $iter: found $num_new new sequence(s) from cmsearch hits.\n";

  ## Step 2: Convergence check
  if ( $num_new == 0 ) {
    print "No new sequences — convergence reached at iteration $iter.\n";
    $iterations_done = $iter;
    last;
  }

  ## Step 3: Extract new sequences and merge with current model FASTA
  my $expand_fa_file = "$model_dir/expand.iter$iter.fa";
  my @new_seq_ids = map { $_->{SEQID} } @{$new_hits};

  open( my $OUT_FA, ">$expand_fa_file" ) or die "Cannot write $expand_fa_file\n";

  ## Write existing model sequences
  my @model_fa = GraphClust::read_fasta_file($model_fa_file);
  foreach my $id ( @{ $model_fa[1] } ) {
    print $OUT_FA ">$id $model_fa[2]->{$id}\n$model_fa[0]->{$id}\n";
  }

  ## Write new sequences from scan FASTA
  foreach my $seq_id (@new_seq_ids) {
    if ( exists $fa_scan[0]->{$seq_id} ) {
      print $OUT_FA ">$seq_id $fa_scan[2]->{$seq_id}\n$fa_scan[0]->{$seq_id}\n";
    } else {
      print "WARNING: sequence $seq_id not found in scan FASTA, skipping.\n";
    }
  }
  close($OUT_FA);

  ## Step 4: Re-align with mlocarna
  my $expand_aln_dir = "$model_dir/expand.iter$iter";
  if ( !-e "$expand_aln_dir/results/result.aln" ) {
    print "Running mlocarna on expanded set (" . ( scalar(@{ $model_fa[1] }) + $num_new ) . " sequences)...\n";
    GraphClust::mlocarna_center( $expand_fa_file, $expand_aln_dir, "$in_cluster_dir/dp", 1 );
  }

  my $expanded_aln = "$expand_aln_dir/results/result.aln";
  if ( !-e $expanded_aln ) {
    print "WARNING: mlocarna failed for iteration $iter, stopping expansion.\n";
    $iterations_done = $iter;
    last;
  }

  ## Run alifold on new alignment
  GraphClust::aln2alifold( $expanded_aln, $tmp_path, $vrna_path );

  ## Step 5: CMfinder refinement (if available)
  my $final_aln = $expanded_aln;

  if ($use_cmfinder) {
    print "Running CMfinder on expanded set...\n";
    my $cmf_stk = "$model_dir/model.cmfinder.iter$iter.stk";
    my $cmf_cm  = "$model_dir/model.cmfinder.iter$iter.cm";

    my $cmf_cmd = "${cmfinder_path}cmfinder --g 1.0 -a $model_dir/model.tree.stk $expand_fa_file $cmf_cm > $cmf_stk";
    system_call( $cmf_cmd, 1 );
    system("rm -f $cmf_cm");

    if ( -e $cmf_stk && !-z $cmf_stk ) {
      open( my $STK_IN, $cmf_stk );
      my $model_stk = GraphClust::read_stockholm($STK_IN);
      close($STK_IN);
      open( my $ALN_OUT, ">$model_dir/model.cmfinder.iter$iter.aln" );
      GraphClust::write_clustal( $ALN_OUT, $model_stk );
      close($ALN_OUT);
      GraphClust::aln2alifold( "$model_dir/model.cmfinder.iter$iter.aln", $tmp_path, $vrna_path );

      ## Check if CMfinder model has structure (negative energy)
      open( my $ALI, "$model_dir/model.cmfinder.iter$iter.aln.alifold" );
      my @cmf_energy = <$ALI>;
      chomp(@cmf_energy);
      close($ALI);

      my $cmf_en = 0;
      if ( $cmf_energy[2] =~ /\(\s*(\S+)\s+=\s+(\S+)\s+\+\s+(\S+)\)/ ) {
        $cmf_en = $1;
      }

      if ( $cmf_en < 0 ) {
        print "CMfinder model energy: $cmf_en — using CMfinder alignment.\n";
        $final_aln = "$model_dir/model.cmfinder.iter$iter.aln";
      } else {
        print "CMfinder model unstable (energy=$cmf_en) — keeping mlocarna alignment.\n";
      }
    } else {
      print "CMfinder produced no output — keeping mlocarna alignment.\n";
    }
  }

  ## Step 6: Create new model.stk
  my $stk_cmd = "$binDir/mloc2stockholm.pl -file $final_aln";
  $stk_cmd .= " -split_input yes -con_struct $final_aln.alifold";
  system_call($stk_cmd);
  system("cp $final_aln.sth $model_dir/model.stk");

  ## Update model.tree.fa with expanded set
  system("cp $expand_fa_file $model_fa_file");

  ## Step 7: cmbuild + cmsearch
  my $iter_dir = "$cmsearch_dir/expand.iter$iter";
  make_path($iter_dir);

  my $cm_file = "$iter_dir/model.stk.cm";
  system_call("$infernal_path/cmbuild -F $cm_file $model_dir/model.stk");

  my $tab_file = "$iter_dir/model.stk.cm.tabresult";
  my $tmp_tab  = "$iter_dir/model.stk.cm.tab_tmp";
  my $cms_cmd = "$infernal_path/cmsearch -g $cmsearch_filter_opt $cmsearch_noalign_opt $cmsearch_cpu_opt ";
  $cms_cmd .= "$cmsearch_table_opt $tmp_tab ";
  $cms_cmd .= "--toponly " if ($cm_top_only);
  my $min_bs = min( 10, $cm_min_bitscore );
  $cms_cmd .= "-T $min_bs " if ( $cm_bitscore_sig == 1 );
  $cms_cmd .= "$cm_file $db";
  system("rm -f $tmp_tab $tab_file");
  system_call($cms_cmd);
  system("mv $tmp_tab $tab_file");

  ## Step 8: Update state
  $current_tabresult = $tab_file;
  foreach my $id (@new_seq_ids) { $known_ids->{$id} = 1; }
  $total_added += $num_new;
  $iterations_done = $iter;

  ## DONE file for restart safety
  system("touch $in_cluster_dir/expand.iter$iter.DONE");
  print "Iteration $iter complete: added $num_new sequences.\n";
}

print "Expansion loop finished after $iterations_done iteration(s), $total_added total sequences added.\n";
```

- [ ] **Step 2: Verify syntax**

Run: `perl -c gc_cm_expand.pl`
Expected: `gc_cm_expand.pl syntax OK`

- [ ] **Step 3: Commit**

```bash
git add gc_cm_expand.pl
git commit -m "feat: gc_cm_expand.pl expansion loop — mlocarna + CMfinder + cmbuild/cmsearch"
```

---

### Task 4: `gc_cm_expand.pl` — R-scape and Stats Output

**Files:**
- Modify: `gc_cm_expand.pl` (append after expansion loop)

- [ ] **Step 1: Append R-scape and stats-writing code**

```perl
################################################################################
## R-scape on final model alignment

my $rscape_pairs  = 0;
my $rscape_signif = 0;
my $rscape_evalue = "-";

if ( $rscape_path && $rscape_path ne "false" && -e "$model_dir/model.stk" ) {
  my $rscape_bin = "${rscape_path}R-scape";
  $rscape_bin = "$rscape_path/R-scape" if ( !-e $rscape_bin && -e "$rscape_path/R-scape" );

  if ( -e $rscape_bin ) {
    print "Running R-scape on $model_dir/model.stk ...\n";
    my $rscape_cmd = "$rscape_bin --outdir $model_dir/ $model_dir/model.stk";
    my $rscape_out = readpipe("$rscape_cmd 2>&1");
    print $rscape_out if $verbose;

    ## Parse R-scape output for covariation statistics
    ## R-scape prints lines like:
    ##   # 12 pairs analyzed
    ##   # 3 pairs with E-value < 0.05
    ##   # min E-value 0.00034
    if ( $rscape_out =~ /(\d+)\s+pairs\s+analyzed/ ) {
      $rscape_pairs = $1;
    }
    if ( $rscape_out =~ /(\d+)\s+pairs\s+with\s+E-value/ ) {
      $rscape_signif = $1;
    }
    ## Find the minimum E-value from significant pairs
    while ( $rscape_out =~ /^\s*\*\s+\S+\s+\S+\s+(\S+)\s+/mg ) {
      my $ev = $1;
      if ( $rscape_evalue eq "-" || $ev < $rscape_evalue ) {
        $rscape_evalue = $ev;
      }
    }

    print "R-scape: $rscape_pairs pairs tested, $rscape_signif significant, min E-value=$rscape_evalue\n";
  } else {
    print "R-scape binary not found at $rscape_bin — skipping.\n";
  }
} else {
  print "R-scape: skipped (not configured or no model.stk).\n";
}

################################################################################
## Write expand.stats

my $stats_file = "$in_cluster_dir/expand.stats";
open( my $STATS, ">$stats_file" ) or die "Cannot write $stats_file\n";
print $STATS "RSCAPE_PAIRS $rscape_pairs RSCAPE_SIGNIF $rscape_signif RSCAPE_EVALUE $rscape_evalue "
           . "EXPAND_ITERS $iterations_done EXPAND_SEQS_ADDED $total_added\n";
close($STATS);
print "Stats written to $stats_file\n";

################################################################################
## Done
print "gc_cm_expand.pl finished successfully.\n";
```

- [ ] **Step 2: Verify syntax**

Run: `perl -c gc_cm_expand.pl`
Expected: `gc_cm_expand.pl syntax OK`

- [ ] **Step 3: Commit**

```bash
git add gc_cm_expand.pl
git commit -m "feat: gc_cm_expand.pl R-scape integration and expand.stats output"
```

---

### Task 5: Modify `MASTER_GraphClust.pl` — Stage 8b Block

**Files:**
- Modify: `MASTER_GraphClust.pl:746-749` (skip-check)
- Modify: `MASTER_GraphClust.pl:977` (insert Stage 8b block after Stage 8)

- [ ] **Step 1: Update the skip-check at line 746**

Replace:
```perl
      if ( -e "$CLUSTER_DIR/$clus_idx.cluster/cmsearch.DONE" ) {
        print "Round $CI cluster $clus_idx stages 6-8 already finished!\n";
        delete $toDo_models{$clus_idx};
        next;
      }
```

With:
```perl
      if ( -e "$CLUSTER_DIR/$clus_idx.cluster/expand.DONE" ) {
        print "Round $CI cluster $clus_idx stages 6-8b already finished!\n";
        delete $toDo_models{$clus_idx};
        next;
      }
```

This changes the "fully done" check from `cmsearch.DONE` to `expand.DONE`, so that
a cluster is only considered complete once expansion (or its skip) has finished.

- [ ] **Step 2: Insert Stage 8b block after Stage 8 (after `}  ## fi stage 8`)**

After line ~977 (`}    ## fi stage 8`), before `}    ## foreach keys %toDo_models`, insert:

```perl
      ####################################################################################################
      ## stage 8b: iterative CM expansion + R-scape

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
          $job_name, "$BIN_DIR/gc_cm_expand.sge", $CMD_expand, 1,
          $SGE_ERR_DIR, $in_USE_SGE,
          "$curr_cluster_dir/SGE_log_expand",
          "$EVAL_DIR/times/time.stage.8b.$clus_idx",
          0, $stage8_qsub_opts, $NUM_THREADS, $job_uuid, undef, $stage8_local_slots
        );

        if ( $sge_status->[0] == 1 && $sge_status->[1] == 0 ) {
          system_call("touch $CLUSTER_DIR/$clus_idx.cluster/expand.DONE");
          $trigger_new_partition = 1;

        } elsif ( $sge_status->[1] == 1 ) {
          print "Round $CI cluster $clus_idx stage 8b: SGE job generated some error! Skip...\n";
          ## Still mark as done so pipeline continues (expansion is enhancement, not critical)
          system_call("touch $CLUSTER_DIR/$clus_idx.cluster/expand.DONE");
          $cluster_error++;
        }
      }

      next if ( !-e "$CLUSTER_DIR/$clus_idx.cluster/expand.DONE" );
```

- [ ] **Step 3: Check that `$stage8_cpu_threads` and `$stage8_qsub_opts` are in scope**

These variables are defined inside the Stage 8 block (around line 949-959). They need to be accessible in the Stage 8b block too. Verify they are declared in the same scope (inside the `foreach` over `%toDo_models`). If they're scoped inside the `if (!-e cmsearch.DONE)` block, move them before that block.

Check lines ~949-959 in MASTER_GraphClust.pl. The variables `$stage8_cpu_threads`, `$stage8_qsub_opts`, `$stage8_local_slots` are defined inside the Stage 8 `if` block. Move the declarations BEFORE the Stage 8 `if` block:

```perl
      ## CPU/SGE settings for stage 8 and 8b
      my $stage8_cpu_threads = 1;
      my $stage8_qsub_opts   = "";
      my $stage8_local_slots = 1;

      if ( $in_USE_SGE && $SGE_PE_THREADS > 1 ) {
        $stage8_cpu_threads = $SGE_PE_THREADS;
        $stage8_qsub_opts .= " -pe \"$in_SGE_PE_name\" 1-$SGE_PE_THREADS ";
      } elsif ( $NUM_THREADS > 1 ) {
        $stage8_cpu_threads = $NUM_THREADS;
        $stage8_local_slots = $NUM_THREADS;
      }
```

This replaces the existing declarations that are currently inside the Stage 8 `if` block.

- [ ] **Step 4: Verify syntax**

Run: `perl -c MASTER_GraphClust.pl`
Expected: `MASTER_GraphClust.pl syntax OK`

- [ ] **Step 5: Commit**

```bash
git add MASTER_GraphClust.pl
git commit -m "feat: add Stage 8b (gc_cm_expand) to MASTER pipeline"
```

---

### Task 6: Modify `gc_results_cluster.pl` — Read expand.stats into cluster.stats

**Files:**
- Modify: `gc_results_cluster.pl:220-231` (after RNAz block, before `print STAT`)

- [ ] **Step 1: Add expand.stats reading after the RNAz block**

After the block at ~line 225 (`print "RNAz: " . join( ":", @{$top5_rnaz} ) . "\n";`) and BEFORE `open( STAT, ">$clus_dir/cluster.stats" );` at ~line 228, insert:

```perl
  ## Read R-scape and CM expansion stats (from Stage 8b)
  my %expand_stats = ( pairs => 0, signif => 0, evalue => '-', iters => 0, added => 0 );
  foreach my $orig_cl (@orig_clus) {
    my $f = "$in_root_dir/CLUSTER/$orig_cl.cluster/expand.stats";
    if ( -e $f ) {
      open( my $EXP, $f );
      my $line = <$EXP>;
      close($EXP);
      chomp $line;
      $expand_stats{pairs}  = $1 if ( $line =~ /RSCAPE_PAIRS (\S+)/ );
      $expand_stats{signif} = $1 if ( $line =~ /RSCAPE_SIGNIF (\S+)/ );
      $expand_stats{evalue} = $1 if ( $line =~ /RSCAPE_EVALUE (\S+)/ );
      $expand_stats{iters}  = $1 if ( $line =~ /EXPAND_ITERS (\S+)/ );
      $expand_stats{added}  = $1 if ( $line =~ /EXPAND_SEQS_ADDED (\S+)/ );
      last;
    }
  }
```

- [ ] **Step 2: Append R-scape and expansion fields to `print STAT`**

Find the `print STAT` line at ~line 230-231:
```perl
  print STAT "CLUSTER $clus_idx SEQS " . scalar(@clus_keys) . " ";
  print STAT "IDS_UNIQUE " . scalar(@ids_unique) . " MODELS " . scalar(@orig_clus) . " MPI_TOP5 $top5_rnaz->[0] SCI_TOP5 $top5_rnaz->[9] ";
```

Add AFTER the second `print STAT` (and before the `if ($evaluate)` block):

```perl
  print STAT "RSCAPE_PAIRS $expand_stats{pairs} RSCAPE_SIGNIF $expand_stats{signif} "
           . "RSCAPE_EVALUE $expand_stats{evalue} "
           . "EXPAND_ITERS $expand_stats{iters} EXPAND_SEQS_ADDED $expand_stats{added} ";
```

- [ ] **Step 3: Verify syntax**

Run: `perl -c gc_results_cluster.pl`
Expected: `gc_results_cluster.pl syntax OK`

- [ ] **Step 4: Commit**

```bash
git add gc_results_cluster.pl
git commit -m "feat: include R-scape and expansion stats in cluster.stats output"
```

---

### Task 7: Add `gc_cm_expand.pl` and `gc_cm_expand.sge` to `Makefile.am` bin_SCRIPTS

**Files:**
- Modify: `Makefile.am`

- [ ] **Step 1: Add the new scripts to `bin_SCRIPTS` and `EXTRA_DIST`**

In `bin_SCRIPTS`, after `gc_cmsearch.pl gc_cmsearch.sge`:
```
gc_cm_expand.pl gc_cm_expand.sge \
```

In `EXTRA_DIST`, after `gc_cmsearch.pl gc_cmsearch.sge`:
```
gc_cm_expand.pl gc_cm_expand.sge \
```

- [ ] **Step 2: Commit**

```bash
git add Makefile.am
git commit -m "build: add gc_cm_expand.pl and .sge to Makefile.am"
```

---

### Task 8: Integration Test with Example Data

**Files:**
- No new files; uses existing `examples/` data

- [ ] **Step 1: Verify all scripts pass syntax check**

```bash
cd /Volumes/Masterarbeit/GraphClust
perl -c gc_cm_expand.pl
perl -c MASTER_GraphClust.pl
perl -c gc_results_cluster.pl
```

Expected: All report `syntax OK`.

- [ ] **Step 2: Create a minimal test fixture**

Create a temporary directory structure mimicking a post-Stage-8 cluster:

```bash
TEST_ROOT=/tmp/gc_expand_test
mkdir -p $TEST_ROOT/FASTA $TEST_ROOT/CLUSTER/1.1.cluster/MODEL $TEST_ROOT/CLUSTER/1.1.cluster/CMSEARCH

# Copy example data if available, or create minimal test files
# This needs actual FASTA and tabresult data from an existing GraphClust run.
# On the cluster, use an existing completed run directory.
```

- [ ] **Step 3: Run gc_cm_expand.pl in dry mode**

On the cluster (after building CMfinder), run on a real completed GraphClust run:

```bash
# Find a completed cluster directory
CLUSTER_DIR=/path/to/real/GraphClust_run/CLUSTER/1.1.cluster

perl gc_cm_expand.pl \
  --root-dir /path/to/real/GraphClust_run \
  --cluster-dir $CLUSTER_DIR \
  --db /path/to/real/GraphClust_run/FASTA/data.fasta.scan \
  --max-iter 1 \
  --cpu 2 \
  --verbose
```

Expected: Script runs, extracts hits, attempts mlocarna, writes `expand.stats`.

- [ ] **Step 4: Verify expand.stats was written**

```bash
cat $CLUSTER_DIR/expand.stats
```

Expected: Line with `RSCAPE_PAIRS ... EXPAND_ITERS ... EXPAND_SEQS_ADDED ...`

- [ ] **Step 5: Final commit**

```bash
git add -A
git commit -m "test: verify gc_cm_expand.pl integration with real data"
```

---

### Task 9: Create `motif_plot.py` — Cluster Motif Visualization in Results Phase

**Files:**
- Create: `motif_plot.py`
- Modify: `gc_results_cluster.pl` (add call to motif_plot.py at end of per-cluster processing)

Adapted from GraphClust2's Galaxy `MotifFinderPlot.py`. Creates a horizontal bar chart showing motif positions per sequence, coloured by cluster. Called by `gc_results_cluster.pl` at the end of Stage 9 for each result cluster.

- [ ] **Step 1: Create `motif_plot.py`**

```python
#!/usr/bin/env python
"""
Cluster motif visualisation for GraphClust results.
Reads RESULTS/*/cluster.all files and produces a motif overview plot (PDF).

Adapted from GraphClust2 Galaxy MotifFinderPlot.py (bgruening/galaxytools).

Usage:
  python motif_plot.py <results_dir> [output_file]

  results_dir  — path to RESULTS/ directory (contains subdirs like 1/, 2/, ...)
  output_file  — optional, default: <results_dir>/motif_plot.pdf
"""

import glob
import itertools
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import pyplot as plt  # noqa: E402

# Use a colourblind-friendly categorical palette (Okabe-Ito inspired)
PALETTE = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#999999",
]


def parse_clusters(results_dir):
    """Read all cluster.all files from results_dir/*/cluster.all."""
    cluster_files = sorted(glob.glob(os.path.join(results_dir, "*", "cluster.all")))
    if not cluster_files:
        print("WARNING: no cluster.all files found in %s" % results_dir)
        return {}, {}, {}, {}

    palette = itertools.cycle(PALETTE)
    ranges = defaultdict(list)
    colors = defaultdict(list)
    orig_names = {}
    cluster_nums = {}

    for cluster_file in cluster_files:
        cluster_color = next(palette)
        with open(cluster_file) as fh:
            for line in fh:
                parts = line.strip().split()
                if len(parts) < 11 or parts[1] != "RESULT":
                    continue
                seq, start, end, strand = parts[0].split("#")
                ranges[seq].append((int(start), max(1, int(end) - int(start) + 1)))
                colors[seq].append(cluster_color)
                cluster_label = "cluster-%s" % parts[2]
                cluster_nums[cluster_label] = cluster_color
                if "ORIGHEAD" in parts:
                    idx = parts.index("ORIGHEAD")
                    if idx + 1 < len(parts):
                        orig_names[seq] = parts[idx + 1]
                elif seq not in orig_names:
                    orig_names[seq] = seq

    return ranges, colors, orig_names, cluster_nums


def plot_bar(ranges, colors, orig_names, cluster_nums, output_file):
    """Create horizontal bar chart of motif positions."""
    if not ranges:
        print("No data to plot.")
        return

    sorted_keys = sorted(ranges.keys())
    fig, ax = plt.subplots(figsize=(10, max(3, 0.4 * len(sorted_keys))))

    for i, k in enumerate(sorted_keys):
        ax.broken_barh(ranges[k], (i - 0.25, 0.5), facecolors=colors[k])

    ax.set_xlim(0)
    ax.set_xlabel("Position in sequence")
    ax.set_yticks(range(len(sorted_keys)))
    short_labels = []
    for k in sorted_keys:
        name = orig_names.get(k, k)
        if len(name) > 30:
            name = name[:27] + "..."
        short_labels.append("%s — %s" % (k, name))
    ax.set_yticklabels(short_labels, fontsize=7)
    ax.grid(True, alpha=0.3)

    patches = [
        mpatches.Patch(color=cluster_nums[lab], label=lab)
        for lab in sorted(cluster_nums)
    ]
    ax.legend(handles=patches, loc="upper right", fontsize=7, framealpha=0.8)

    plt.tight_layout()
    plt.savefig(output_file, bbox_inches="tight")
    plt.close()
    print("Motif plot saved to %s" % output_file)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python motif_plot.py <results_dir> [output_file]")
        sys.exit(1)

    results_dir = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else os.path.join(results_dir, "motif_plot.pdf")

    r, c, n, cn = parse_clusters(results_dir)
    plot_bar(r, c, n, cn, output_file)
```

- [ ] **Step 2: Make executable**

Run: `chmod +x motif_plot.py`

- [ ] **Step 3: Add call in `gc_results_cluster.pl`**

At the very end of `gc_results_cluster.pl` (after the `foreach` loop over `@res_todo` closes, before the script ends), add:

```perl
## Generate motif overview plot (requires Python + matplotlib)
my $motif_plot_script = "$bin_dir/motif_plot.py";
if ( -e $motif_plot_script ) {
  my $plot_out = "$in_root_dir/RESULTS/motif_plot.pdf";
  system_call("python $motif_plot_script $in_root_dir/RESULTS $plot_out", 1);
  print "Motif plot: $plot_out\n" if ( -e $plot_out );
}
```

- [ ] **Step 4: Add to `Makefile.am`**

In `bin_SCRIPTS`, add:
```
motif_plot.py \
```

In `EXTRA_DIST`, add:
```
motif_plot.py \
```

- [ ] **Step 5: Verify**

Run: `python motif_plot.py --help 2>&1 || python3 motif_plot.py --help 2>&1`
Expected: Usage line printed.

- [ ] **Step 6: Commit**

```bash
git add motif_plot.py gc_results_cluster.pl Makefile.am
git commit -m "feat: add motif_plot.py for cluster visualisation in results phase"
```
