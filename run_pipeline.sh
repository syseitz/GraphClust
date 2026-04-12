#!/bin/bash
# Run GraphClust pipeline with rnaclust conda environment.
# Works on both macOS and Linux — auto-detects conda location.
# Usage: ./run_pipeline.sh [pipeline args...]

# Find conda
CONDA_SH=""
for candidate in \
  "$HOME/miniforge3/etc/profile.d/conda.sh" \
  "$HOME/miniconda3/etc/profile.d/conda.sh" \
  "$HOME/.local/share/mamba/etc/profile.d/conda.sh" \
  "$CONDA_PREFIX/../etc/profile.d/conda.sh" \
  "/opt/conda/etc/profile.d/conda.sh"; do
  if [ -f "$candidate" ]; then
    CONDA_SH="$candidate"
    break
  fi
done

if [ -z "$CONDA_SH" ]; then
  # Try to find it via which conda
  CONDA_BIN=$(which conda 2>/dev/null)
  if [ -n "$CONDA_BIN" ]; then
    CONDA_SH="$(dirname "$(dirname "$CONDA_BIN")")/etc/profile.d/conda.sh"
  fi
fi

if [ -z "$CONDA_SH" ] || [ ! -f "$CONDA_SH" ]; then
  echo "ERROR: Cannot find conda. Please ensure conda/miniforge/miniconda is installed."
  exit 1
fi

source "$CONDA_SH"
# Try activating by name first, then search common mamba/conda env paths
conda activate rnaclust 2>/dev/null || {
  for envdir in \
    "$HOME/.local/share/mamba/envs/rnaclust" \
    "$HOME/miniforge3/envs/rnaclust" \
    "$HOME/miniconda3/envs/rnaclust" \
    "$HOME/mambaforge/envs/rnaclust"; do
    if [ -d "$envdir" ]; then
      conda activate "$envdir" && break
    fi
  done
  if [ -z "$CONDA_PREFIX" ] || [[ "$CONDA_PREFIX" != *rnaclust* ]]; then
    echo "ERROR: Cannot activate rnaclust environment"
    exit 1
  fi
}

# Ensure conda env binaries come before Homebrew/system perl
export PATH="$CONDA_PREFIX/bin:$PATH"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Use patched mlocarna + MatchProbs.pm from locarna_patches/ if available.
# These contain fixes for write_reliability_bars and multi_exists (shared hash bug).
PATCHES_DIR="$SCRIPT_DIR/locarna_patches"
if [ -d "$PATCHES_DIR/bin" ]; then
  export PATH="$PATCHES_DIR/bin:$PATH"
fi
if [ -d "$PATCHES_DIR/lib/perl" ]; then
  export PERL5LIB="$PATCHES_DIR/lib/perl${PERL5LIB:+:$PERL5LIB}"
fi

# Set Perl lib path for remaining MLocarna modules from conda env
LOCARNA_BASE=$(dirname "$(which mlocarna 2>/dev/null)")
if [ -n "$LOCARNA_BASE" ] && [ -d "$LOCARNA_BASE/../lib/perl" ]; then
  export PERL5LIB="${PERL5LIB:+$PERL5LIB:}${LOCARNA_BASE}/../lib/perl"
fi

perl "$SCRIPT_DIR/MASTER_GraphClust.pl" "$@"
