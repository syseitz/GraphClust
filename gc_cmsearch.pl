#!/usr/bin/perl -w
use warnings;
use strict;

## GraphClust stage 8 script for cmsearch

use List::Util qw/ min max /;
use FindBin;
use lib "$FindBin::Bin";
use GraphClust;

use vars qw(%CONFIG);
*CONFIG = \%GraphClust::CONFIG;

use Getopt::Long;

my $binDir = "$FindBin::Bin";

my $stk_file;
my $db;
my $tgtdir;
my $evaluate;
my $use_cmcalibrate;
my $in_root_dir;
my $verbose;
my $cpu_threads = 1;


## infernal 1.1
# my $OPTS_cmsearch = " -g --rfam --noali --cpu 1";
## infernal 1.1: --refine
# my $OPTS_cmbuild = "";
# my $OPTS_cmcalibrate = " -L 0.01 --cpu 1 ";

## infernal 1.0.2 / 1.1 compatibility
my $cmsearch_help = `cmsearch -h 2>&1`;
if ( $CONFIG{PATH_INFERNAL} && $CONFIG{PATH_INFERNAL} ne "false" ) {
  $cmsearch_help = `$CONFIG{PATH_INFERNAL}/cmsearch -h 2>&1`;
}
my $cmcalibrate_help = `cmcalibrate -h 2>&1`;
if ( $CONFIG{PATH_INFERNAL} && $CONFIG{PATH_INFERNAL} ne "false" ) {
  $cmcalibrate_help = `$CONFIG{PATH_INFERNAL}/cmcalibrate -h 2>&1`;
}
my $cmsearch_filter_opt = "--fil-no-hmm";
$cmsearch_filter_opt = "--nohmm" if ( $cmsearch_help !~ /--fil-no-hmm/ && $cmsearch_help =~ /--nohmm/ );
my $cmsearch_table_opt = "--tabfile";
$cmsearch_table_opt = "--tblout" if ( $cmsearch_help =~ /--tblout/ );
my $cmsearch_noalign_opt = "--noalign";
$cmsearch_noalign_opt = "--noali" if ( $cmsearch_help !~ /--noalign/ && $cmsearch_help =~ /--noali/ );
my $cmsearch_cpu_opt = "";
my $cmcalibrate_cpu_opt = "";

my $OPTS_cmsearch = " -g $cmsearch_filter_opt $cmsearch_noalign_opt ";
my $OPTS_cmbuild = "";
my $OPTS_cmcalibrate = " -L 0.01 ";

usage()
  unless GetOptions(
  "root=s"           => \$in_root_dir,
  "tgtdir=s"         => \$tgtdir,
  "stk=s"            => \$stk_file,
  "db=s"             => \$db,
  "cpu=i"            => \$cpu_threads,
  "calibrate"        => \$use_cmcalibrate,
  "verbose"          => \$verbose
  );

################################################################################
## Load options from root-dir/config file
%GraphClust::CONFIG = readConfigFile( "$in_root_dir/config", 1 );
SUBSECTION("options loaded for gc_cmsearch.pl");
printConfig( \%CONFIG ) if ($verbose);

my $cm_top_only = 1;
$cm_top_only = $CONFIG{cm_top_only};
my $infernal_path = $CONFIG{PATH_INFERNAL};
my $cm_bitscore_sig = $CONFIG{cm_bitscore_sig};
my $cm_min_bitscore = $CONFIG{cm_min_bitscore};

$cpu_threads = 1 if ( !$cpu_threads || $cpu_threads < 1 );
$cmsearch_cpu_opt = " --cpu $cpu_threads " if ( $cmsearch_help =~ /--cpu/ );
$cmcalibrate_cpu_opt = " --cpu $cpu_threads " if ( $cmcalibrate_help =~ /--cpu/ );
$OPTS_cmsearch .= $cmsearch_cpu_opt;
$OPTS_cmcalibrate .= $cmcalibrate_cpu_opt;

## check top_only option consistency
## if we have revcompl seqs than we should also scan on them
## remember: scan is done on data.fasta.scan which are the original input seqs
## and they dont contain the added revcompl seqs in general
if ( $CONFIG{input_add_revcompl} && $cm_top_only ) {
  print "\nATTENTION! Option 'cm_top_only' changed to 0 (false)!\n";
  print " This is necessary because option 'input_add_revcompl' is 1 (true)!\n\n";
  $cm_top_only = 0;
}

################################################################################

mkdir($tgtdir);

die "No input file specified....\n\n" if ( !$stk_file || !-e $stk_file);
die "No input database specified....\n" if ( !$db || !-e $db );

################################
### build the covariance model |
################################

my @file_line = split( /\//, $stk_file );
my $target_file = $file_line[ scalar(@file_line) - 1 ];

my $cmd_cmbuild = "$infernal_path/cmbuild $OPTS_cmbuild -F $tgtdir/$target_file.cm $stk_file";

system($cmd_cmbuild);

####################################
### calibrate the covariance model |
####################################

if ($use_cmcalibrate ) {

  my $cmd_cal = "$infernal_path/cmcalibrate $OPTS_cmcalibrate $tgtdir/$target_file.cm";
  system($cmd_cal);
}
#####################################################################################################
### search the database with the calibrated covariance model and report the findings into a tabfile |
#####################################################################################################
### Each non-# prefixed line of this file corresponds to a hit, and each such line has 9 fields:    |
### <model name> the name of the CM used for the search, 					    |
### <target name> the name of the target sequence the hit was found in, 			    |
### <target coord - start> the start position of the hit in the target sequence, 		    |
### <target coord - stop> the end position of hit in the target sequence, 			    |
### <query coord - start> the start position of the hit in the query model, 			    |
### <query coord - stop> the end position of hit in the query sequence, 			    |
### <bit sc> the bit score of the hit, 								    |
### <E-value> the E-value of the hit (if available, ”-” if not),				    |
### <GC> the percentage of G and C residues in the hit within the target sequence.		    |
#####################################################################################################
## add '-o $e/$target_file.cm.alignments' if also alignment file should be created

my $tmp_tab_file = "$tgtdir/$target_file.cm.tab_tmp";
my $cmd_cms = "$infernal_path/cmsearch $OPTS_cmsearch $cmsearch_table_opt $tmp_tab_file ";
$cmd_cms .= "--toponly " if ($cm_top_only);
my $min_bitscore = min(10,$cm_min_bitscore);
$cmd_cms .= "-T $min_bitscore " if ($cm_bitscore_sig == 1);
$cmd_cms .= " $tgtdir/$target_file.cm $db ";

system("rm -f $tmp_tab_file");
system("rm -f $tgtdir/$target_file.cm.tabresult");

## run cmsearch
system( $cmd_cms );

## postprocessing
system("mv $tmp_tab_file $tgtdir/$target_file.cm.tabresult ");
