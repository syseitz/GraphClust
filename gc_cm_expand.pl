#!/usr/bin/perl -w
use warnings;
use strict;

## GraphClust post-stage-8 script for iterative cluster expansion via cmsearch

use List::Util qw/ min max /;
use FindBin;
use lib "$FindBin::Bin";
use GraphClust;

use vars qw(%CONFIG);
*CONFIG = \%GraphClust::CONFIG;

use Getopt::Long;
use Cwd qw/ abs_path /;
use File::Path qw/ make_path /;

my $binDir = "$FindBin::Bin";

################################################################################
## Argument parsing
################################################################################

my $in_root_dir;
my $in_cluster_dir;
my $db;
my $max_iter    = 3;
my $cpu_threads = 1;
my $verbose;

usage()
  unless GetOptions(
  "root-dir=s"    => \$in_root_dir,
  "cluster-dir=s" => \$in_cluster_dir,
  "db=s"          => \$db,
  "max-iter=i"    => \$max_iter,
  "cpu=i"         => \$cpu_threads,
  "verbose"       => \$verbose,
  "task-id=i"     => \(my $task_id),  # passed by SGE wrapper, not used directly
  );

################################################################################
## Load config from root-dir
################################################################################

%GraphClust::CONFIG = readConfigFile( "$in_root_dir/config", 1 );
SUBSECTION("options loaded for gc_cm_expand.pl");
printConfig( \%CONFIG ) if ($verbose);

################################################################################
## Extract config variables
################################################################################

my $infernal_path  = $CONFIG{PATH_INFERNAL};
my $cmfinder_path  = $CONFIG{PATH_CMFINDER};
my $rscape_path    = $CONFIG{PATH_RSCAPE};
my $vrna_path      = $CONFIG{PATH_VRNA};
my $tmp_path       = $CONFIG{PATH_TMP};
my $cm_min_bitscore = $CONFIG{cm_min_bitscore};
my $cm_bitscore_sig = $CONFIG{cm_bitscore_sig};
my $cm_max_eval     = $CONFIG{cm_max_eval};
my $cm_top_only     = $CONFIG{cm_top_only};

$cpu_threads = 1 if ( !$cpu_threads || $cpu_threads < 1 );

################################################################################
## Infernal 1.0.2 / 1.1 compatibility detection
################################################################################

my $cmsearch_help = `cmsearch -h 2>&1`;
if ( $CONFIG{PATH_INFERNAL} && $CONFIG{PATH_INFERNAL} ne "false" ) {
  $cmsearch_help = `$CONFIG{PATH_INFERNAL}/cmsearch -h 2>&1`;
}
my $cmcalibrate_help = `cmcalibrate -h 2>&1`;
if ( $CONFIG{PATH_INFERNAL} && $CONFIG{PATH_INFERNAL} ne "false" ) {
  $cmcalibrate_help = `$CONFIG{PATH_INFERNAL}/cmcalibrate -h 2>&1`;
}

my $cmsearch_filter_opt = "--fil-no-hmm";
$cmsearch_filter_opt = "--nohmm"
  if ( $cmsearch_help !~ /--fil-no-hmm/ && $cmsearch_help =~ /--nohmm/ );

my $cmsearch_table_opt = "--tabfile";
$cmsearch_table_opt = "--tblout" if ( $cmsearch_help =~ /--tblout/ );

my $cmsearch_noalign_opt = "--noalign";
$cmsearch_noalign_opt = "--noali"
  if ( $cmsearch_help !~ /--noalign/ && $cmsearch_help =~ /--noali/ );

my $cmsearch_cpu_opt    = "";
my $cmcalibrate_cpu_opt = "";
$cmsearch_cpu_opt    = " --cpu $cpu_threads " if ( $cmsearch_help    =~ /--cpu/ );
$cmcalibrate_cpu_opt = " --cpu $cpu_threads " if ( $cmcalibrate_help =~ /--cpu/ );

my $OPTS_cmsearch   = " -g $cmsearch_filter_opt $cmsearch_noalign_opt $cmsearch_cpu_opt";
my $OPTS_cmbuild    = "";
my $OPTS_cmcalibrate = " -L 0.01 $cmcalibrate_cpu_opt";

################################################################################
## CMfinder availability check
################################################################################

my $use_cmfinder = 0;
if ( $cmfinder_path && $cmfinder_path ne "false" ) {
  if ( -e "${cmfinder_path}cmfinder" || -e "$cmfinder_path/cmfinder" ) {
    $use_cmfinder = 1;
  }
}

################################################################################
## Validate required arguments
################################################################################

die "No --root-dir specified or directory does not exist.\n\n"
  if ( !$in_root_dir || !-d $in_root_dir );
die "No --cluster-dir specified or directory does not exist.\n\n"
  if ( !$in_cluster_dir || !-d $in_cluster_dir );
die "No --db specified or file does not exist.\n\n"
  if ( !$db || !-e $db );

################################################################################
## Print startup info
################################################################################

print "gc_cm_expand.pl -- iterative cluster expansion\n";
print "  cluster-dir : $in_cluster_dir\n";
print "  max-iter    : $max_iter\n";
print "  db          : $db\n";
print "  use_cmfinder: $use_cmfinder\n";

################################################################################
## Helper: collect sequence IDs present in a FASTA file
################################################################################

sub get_fasta_ids {
  my ($fa_file) = @_;

  my @result   = GraphClust::read_fasta_file($fa_file);
  my $ids_aref = $result[1];    ## element [1] is \@order (array of IDs)

  my %id_set = ();
  foreach my $id ( @{$ids_aref} ) {
    $id_set{$id} = 1;
  }
  return \%id_set;
}

################################################################################
## Helper: parse a cmsearch tabresult and return hits not yet in the cluster
################################################################################

sub extract_new_hits {
  my ( $tabresult_file, $known_ids_href ) = @_;

  my $all_hits = GraphClust::read_CM_tabfile_ext(
    $tabresult_file,
    $cm_min_bitscore,
    $cm_max_eval,
    $cm_bitscore_sig,
    "expand"
  );

  return [] if ( !$all_hits || !@{$all_hits} );

  ## Filter out hits already present in the cluster; deduplicate by SEQID
  my %seen_seqid = ();
  my @new_hits   = ();
  foreach my $hit ( @{$all_hits} ) {
    my $seqid = $hit->{SEQID};
    next if ( exists $known_ids_href->{$seqid} );
    next if ( exists $seen_seqid{$seqid} );
    $seen_seqid{$seqid} = 1;
    push( @new_hits, $hit );
  }

  return \@new_hits;
}

################################################################################
## Derived paths
################################################################################

my $model_dir    = "$in_cluster_dir/MODEL";
my $cmsearch_dir = "$in_cluster_dir/CMSEARCH";

################################################################################
## State initialisation
################################################################################

SUBSECTION("gc_cm_expand.pl -- initialising expansion state");

my $known_ids       = get_fasta_ids("$model_dir/model.tree.fa");
my ($current_tabresult) = glob("$cmsearch_dir/*.tabresult");

die "No initial tabresult found in $cmsearch_dir\n" if ( !$current_tabresult || !-e $current_tabresult );

my @fa_scan      = GraphClust::read_fasta_file($db);
## fa_scan[0] = \%seqs  (seqid -> sequence)
## fa_scan[2] = \%headers (seqid -> header text after the id)

my $total_added    = 0;
my $iterations_done = 0;

print "  initial tabresult : $current_tabresult\n";
print "  sequences in model: " . scalar( keys %{$known_ids} ) . "\n";
print "  sequences in db   : " . scalar( keys %{ $fa_scan[0] } ) . "\n\n";

################################################################################
## Expansion loop
################################################################################

for my $iter ( 1 .. $max_iter ) {

  SUBSECTION("Expansion iteration $iter / $max_iter");

  ##############################################################################
  ## Restart safety: skip if this iteration was already completed
  ##############################################################################

  my $done_flag = "$in_cluster_dir/expand.iter$iter.DONE";

  if ( -e $done_flag ) {
    print "  expand.iter$iter.DONE exists — skipping, updating state.\n";

    ## update tabresult pointer to what this iteration produced
    my ($iter_tab) = glob("$cmsearch_dir/expand.iter$iter/*.tabresult");
    $current_tabresult = $iter_tab if ( $iter_tab && -e $iter_tab );

    ## re-read known IDs (the model.tree.fa was updated at the end of this iter)
    $known_ids = get_fasta_ids("$model_dir/model.tree.fa");

    $iterations_done = $iter;
    next;
  }

  ##############################################################################
  ## Step 1 — extract new hits not yet in the cluster
  ##############################################################################

  my $new_hits = extract_new_hits( $current_tabresult, $known_ids );
  my $num_new  = scalar( @{$new_hits} );

  print "  new hits found: $num_new\n";

  ##############################################################################
  ## Step 2 — convergence check
  ##############################################################################

  if ( $num_new == 0 ) {
    print "  No new hits — expansion converged after $iter iteration(s).\n";
    $iterations_done = $iter;
    last;
  }

  ##############################################################################
  ## Step 3 — build expanded FASTA combining known model seqs and new hits
  ##############################################################################

  my $expand_fa_file = "$model_dir/expand.iter$iter.fa";

  ## read current model sequences
  my @model_fa = GraphClust::read_fasta_file("$model_dir/model.tree.fa");

  open( my $FA_OUT, ">$expand_fa_file" )
    or die "Cannot write $expand_fa_file: $!\n";

  ## write sequences already in the model
  foreach my $id ( @{ $model_fa[1] } ) {
    my $head = $model_fa[2]->{$id};
    $head = "" if ( !defined $head );
    print $FA_OUT ">$id $head\n$model_fa[0]->{$id}\n";
  }

  ## write the newly discovered sequences from the scan database
  foreach my $hit ( @{$new_hits} ) {
    my $seqid = $hit->{SEQID};
    if ( !exists $fa_scan[0]->{$seqid} ) {
      warn "  WARNING: SEQID $seqid not found in scan database, skipping.\n";
      next;
    }
    my $head = $fa_scan[2]->{$seqid};
    $head = "" if ( !defined $head );
    print $FA_OUT ">$seqid $head\n$fa_scan[0]->{$seqid}\n";
  }

  close($FA_OUT);

  print "  Wrote expanded FASTA: $expand_fa_file\n";

  ##############################################################################
  ## Step 4 — multiple alignment with mlocarna (locarna-P, use_locP = 1)
  ##############################################################################

  my $mloc_dir = "$model_dir/expand.iter$iter";
  my $dp_dir   = "$in_cluster_dir/dp";

  GraphClust::mlocarna_center( $expand_fa_file, $mloc_dir, $dp_dir, 1 );

  my $result_aln = "$mloc_dir/results/result.aln";

  unless ( -e $result_aln ) {
    warn "  WARNING: mlocarna did not produce result.aln for iteration $iter — stopping expansion.\n";
    last;
  }

  GraphClust::aln2alifold( $result_aln, $tmp_path, $vrna_path );

  ## final_aln starts as the mlocarna result; may be overridden by CMfinder below
  my $final_aln = $result_aln;

  ##############################################################################
  ## Step 5 — optional CMfinder refinement
  ##############################################################################

  if ($use_cmfinder) {

    my $cmf_stk = "$model_dir/expand.iter$iter.cmfinder.stk";
    my $cmf_cm  = "$model_dir/expand.iter$iter.cmfinder.cm";

    system_call(
      "${cmfinder_path}cmfinder --g 1.0 -a $model_dir/model.tree.stk $expand_fa_file $cmf_cm > $cmf_stk",
      1
    );

    system("rm -f $cmf_cm");

    if ( -e $cmf_stk && !-z $cmf_stk ) {

      open( my $STK_IN, $cmf_stk ) or die "Cannot open $cmf_stk: $!\n";
      my $model_stk = GraphClust::read_stockholm($STK_IN);
      close($STK_IN);

      my $cmf_aln = "$model_dir/expand.iter$iter.cmfinder.aln";
      open( my $ALN_OUT, ">$cmf_aln" ) or die "Cannot write $cmf_aln: $!\n";
      GraphClust::write_clustal( $ALN_OUT, $model_stk );
      close($ALN_OUT);

      GraphClust::aln2alifold( $cmf_aln, $tmp_path, $vrna_path );

      open( my $ALI, "$cmf_aln.alifold" ) or die "Cannot open $cmf_aln.alifold: $!\n";
      my @cmf_energy_lines = <$ALI>;
      chomp(@cmf_energy_lines);
      close($ALI);

      my $cmf_energy = 0;
      if ( $cmf_energy_lines[2] =~ /\(\s*(\S+)\s+=\s+(\S+)\s+\+\s+(\S+)\)/ ) {
        $cmf_energy = $1;
      }

      if ( $cmf_energy < 0 ) {
        print "  CMfinder model energy: $cmf_energy — using CMfinder alignment.\n";
        $final_aln = $cmf_aln;
      } else {
        print "  CMfinder energy $cmf_energy >= 0 — keeping mlocarna alignment.\n";
      }
    }
  }

  ##############################################################################
  ## Step 6 — build new model.stk and update model.tree.fa
  ##############################################################################

  system_call(
    "perl $binDir/mloc2stockholm.pl -file $final_aln -split_input yes -con_struct $final_aln.alifold",
    $verbose
  );

  ## mloc2stockholm produces $final_aln with .aln replaced by .sth
  ( my $sth_file = $final_aln ) =~ s/\.aln$/.sth/;

  system("cp $sth_file $model_dir/model.stk");

  ## update model.tree.fa with the expanded sequence set for the next iteration
  system("cp $expand_fa_file $model_dir/model.tree.fa");

  ##############################################################################
  ## Step 7 — cmbuild + cmsearch
  ##############################################################################

  my $iter_search_dir = "$cmsearch_dir/expand.iter$iter";
  make_path($iter_search_dir);

  my $cm_file  = "$iter_search_dir/model.stk.cm";
  my $tmp_tab  = "$iter_search_dir/model.stk.cm.tab_tmp";
  my $tab_file = "$iter_search_dir/model.stk.cm.tabresult";

  system_call( "$infernal_path/cmbuild -F $cm_file $model_dir/model.stk", $verbose );

  my $cmd_cms = "$infernal_path/cmsearch -g $cmsearch_filter_opt $cmsearch_noalign_opt $cmsearch_cpu_opt";
  $cmd_cms .= " $cmsearch_table_opt $tmp_tab";
  $cmd_cms .= " --toponly"        if ($cm_top_only);
  my $min_bitscore = List::Util::min( 10, $cm_min_bitscore );
  $cmd_cms .= " -T $min_bitscore" if ( $cm_bitscore_sig == 1 );
  $cmd_cms .= " $cm_file $db";

  system("rm -f $tmp_tab");
  system_call( $cmd_cms, $verbose );
  system("mv $tmp_tab $tab_file");

  ##############################################################################
  ## Step 8 — update state and mark iteration complete
  ##############################################################################

  $current_tabresult = $tab_file;

  ## add all newly incorporated SEQIDs to the known-ids set
  foreach my $hit ( @{$new_hits} ) {
    $known_ids->{ $hit->{SEQID} } = 1;
  }

  $total_added    += $num_new;
  $iterations_done = $iter;

  system("touch $done_flag");

  print "  Iteration $iter done — added $num_new sequences (total so far: $total_added).\n";
}

################################################################################

print "Expansion loop finished after $iterations_done iteration(s), $total_added total sequences added.\n";

################################################################################
## R-scape covariation analysis
################################################################################

my $rscape_pairs  = 0;
my $rscape_signif = 0;
my $rscape_evalue = "-";

if ( defined $rscape_path && $rscape_path ne "false" && -f "$model_dir/model.stk" ) {

  ## locate the R-scape binary — try bare concatenation first, then with a slash
  my $rscape_bin = "";
  if ( -x "${rscape_path}R-scape" ) {
    $rscape_bin = "${rscape_path}R-scape";
  } elsif ( -x "$rscape_path/R-scape" ) {
    $rscape_bin = "$rscape_path/R-scape";
  }

  if ( $rscape_bin ne "" ) {
    print "Running R-scape: $rscape_bin --outdir $model_dir/ $model_dir/model.stk\n" if $verbose;

    my $rscape_out = readpipe("$rscape_bin --outdir $model_dir/ $model_dir/model.stk 2>&1");

    if ($verbose) {
      print "R-scape output:\n$rscape_out\n";
    }

    ## parse total pairs analysed
    if ( $rscape_out =~ /(\d+)\s+pairs\s+analyzed/ ) {
      $rscape_pairs = $1;
    }

    ## parse number of significant pairs
    if ( $rscape_out =~ /(\d+)\s+pairs\s+with\s+E-value/ ) {
      $rscape_signif = $1;
    }

    ## find minimum E-value from lines for significant pairs (marked with *)
    my $min_eval = undef;
    while ( $rscape_out =~ /^\s*\*\s+\S+\s+\S+\s+(\S+)\s+/mg ) {
      my $ev = $1;
      if ( !defined $min_eval || $ev < $min_eval ) {
        $min_eval = $ev;
      }
    }
    $rscape_evalue = $min_eval if defined $min_eval;

    print "R-scape summary: $rscape_pairs pairs analysed, $rscape_signif significant, "
        . "best E-value = $rscape_evalue\n";
  } else {
    print "R-scape binary not found at $rscape_path — skipping covariation analysis.\n" if $verbose;
  }
} else {
  print "R-scape skipped (path not set or model.stk absent).\n" if $verbose;
}

################################################################################
## Write expand.stats
################################################################################

open( my $stats_fh, ">", "$in_cluster_dir/expand.stats" )
    or die "Cannot write expand.stats: $!";
print $stats_fh "RSCAPE_PAIRS $rscape_pairs "
              . "RSCAPE_SIGNIF $rscape_signif "
              . "RSCAPE_EVALUE $rscape_evalue "
              . "EXPAND_ITERS $iterations_done "
              . "EXPAND_SEQS_ADDED $total_added\n";
close($stats_fh);

print "gc_cm_expand.pl finished successfully.\n";

