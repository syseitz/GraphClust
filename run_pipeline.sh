#!/bin/bash
# Run GraphClust pipeline with rnaclust conda environment.
# Usage: ./run_pipeline.sh [pipeline args...]

# Clean stale mamba lockfiles
rm -f "$HOME/.cache/mamba/proc"/*.json 2>/dev/null
rm -f "$HOME/.cache/mamba/proc/proc.lock" 2>/dev/null

# Activate rnaclust environment using full path
eval "$(/Users/yannick/miniforge3/bin/conda shell.bash hook)"
conda activate /Users/yannick/.local/share/mamba/envs/rnaclust

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
perl "$SCRIPT_DIR/MASTER_GraphClust.pl" "$@"
