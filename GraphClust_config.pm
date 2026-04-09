package GraphClust_config;

use strict;
use warnings;

require Exporter;

################################################################################
## default config

our %CONFIG = (

  PATH_LOCARNA       => "/Users/yannick/.local/share/mamba/envs/rnaclust/bin/",
  PATH_RNASHAPES     => "/Users/yannick/.local/share/mamba/envs/rnaclust/bin/",
  PATH_VRNA          => "/Users/yannick/.local/share/mamba/envs/rnaclust/bin/",
  PATH_INFERNAL      => "/Users/yannick/.local/share/mamba/envs/rnaclust/bin/",
  PATH_R             => "/Users/yannick/.local/share/mamba/envs/rnaclust/bin/",
  PATH_RNAZ          => "/Users/yannick/.local/share/mamba/envs/rnaclust/bin/",
  PATH_RSCAPE        => "/Users/yannick/.local/share/mamba/envs/rnaclust/bin/",
  PATH_BLASTCLUST    => "/Users/yannick/.local/share/mamba/envs/rnaclust/bin/",
  PATH_MMSEQS2       => "false",
  PATH_OCTAVE        => "false",
  cm_expand_max_iter => 3,
  ## CMfinder04: bundled source in cmfinder_src/, build with: bash cmfinder_src/build.sh
  PATH_CMFINDER      => "/Volumes/Masterarbeit/GraphClust/cmfinder_src/bin/",

  VERSION_INFO       =>  "GraphClust 0.7.6",

  ## paths not automatically configured
  ## please adapt to your own system
  PATH_TMP           => "/var/tmp/",
  ## path to qsub for SGE/OGE job submission (sun/oracle grid engine)
  PATH_SGE           => "/opt/sge-6.0/bin/lx24-amd64/",
  ## hostname which is allowed to submit SGE/OGE jobs, ssh is used to login
  SGE_HOSTNAME       => "biui",
);
