#!/bin/bash
# Reset GraphClust pipeline to a specific round and stage.
# Usage: ./reset_pipeline.sh --root-dir <path> --round <N> --stage <N>
#
# Examples:
#   ./reset_pipeline.sh --root-dir /path/to/GraphClust_80 --round 1 --stage 8
#     → Keeps Stage 0-7 of Round 1, deletes Stage 8+ and all later rounds
#
#   ./reset_pipeline.sh --root-dir /path/to/GraphClust_80 --round 2 --stage 5
#     → Keeps Round 1 complete, deletes Stage 5+ of Round 2 and Round 3
#
#   ./reset_pipeline.sh --root-dir /path/to/GraphClust_80 --round 1 --stage 2
#     → Keeps Stage 0-1, deletes GSPAN, SVECTOR, CLUSTER, everything
#
# Stage reference:
#   0  = Initialisation (FASTA dir)
#   1  = Preprocessing (data.fasta)
#   2  = Graph Generation (GSPAN)
#   3  = NSPDK Vectors (SVECTOR)
#   5  = NSPDK Clustering (fast_cluster)
#   6  = LocARNA Align (pp.DONE, cluster.DONE)
#   7  = Build Model (model.DONE)
#   8  = CM Search (cmsearch.DONE, expand.DONE)
#   9  = Collect Results (round.DONE)
#  10  = Final Results (final_partition.soft)

set -euo pipefail

# Parse arguments
ROOT=""
ROUND=""
STAGE=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --root-dir) ROOT="$2"; shift 2 ;;
    --round)    ROUND="$2"; shift 2 ;;
    --stage)    STAGE="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

if [[ -z "$ROOT" || -z "$ROUND" || -z "$STAGE" ]]; then
  echo "Usage: $0 --root-dir <path> --round <N> --stage <N>"
  exit 1
fi

if [[ ! -d "$ROOT" ]]; then
  echo "ERROR: Root directory does not exist: $ROOT"
  exit 1
fi

# Safety: check no pipeline is running
if pgrep -f "MASTER_GraphClust.*$(basename "$ROOT")" > /dev/null 2>&1; then
  echo "ERROR: Pipeline is still running! Stop it first."
  exit 1
fi

echo "========================================"
echo " GraphClust Pipeline Reset"
echo " Root:  $ROOT"
echo " Reset: Round $ROUND, Stage $STAGE"
echo "========================================"
echo ""

# Read config for num_clusters
NUM_CLUSTERS=$(grep "^GLOBAL_num_clusters" "$ROOT/config" | awk '{print $2}')
MAX_ROUNDS=$(grep "^GLOBAL_iterations" "$ROOT/config" | awk '{print $2}')
echo "Config: $NUM_CLUSTERS clusters, $MAX_ROUNDS rounds"
echo ""

########################################
# Delete all rounds > ROUND
########################################

for r in $(seq $((ROUND + 1)) $MAX_ROUNDS); do
  if ls "$ROOT/CLUSTER/$r."*.cluster 1>/dev/null 2>&1; then
    echo "Deleting Round $r cluster directories..."
    find "$ROOT/CLUSTER" -maxdepth 1 -name "$r.*.cluster" -type d -exec rm -rf {} +
  fi
  rm -f "$ROOT/$r.round.DONE"
  rm -f "$ROOT/SVECTOR/data.svector.fast_cluster.$r.DONE"
  rm -f "$ROOT/SVECTOR/data.svector.blacklist.$r"
  rm -f "$ROOT/SVECTOR/data.svector.$r.fast_cluster"
  rm -f "$ROOT/SVECTOR/data.svector.$r.fast_cluster_sim"
done

########################################
# Delete final results (Stage 10)
########################################

if [[ $STAGE -le 10 ]]; then
  echo "Deleting final results..."
  rm -rf "$ROOT/RESULTS/partitions"
fi

########################################
# Delete round.DONE for current round (Stage 9)
########################################

if [[ $STAGE -le 9 ]]; then
  echo "Deleting Round $ROUND round.DONE marker..."
  rm -f "$ROOT/$ROUND.round.DONE"
  rm -rf "$ROOT/EVAL/partitions"
fi

########################################
# Delete Stage 8 (cmsearch + expand) for current round
########################################

if [[ $STAGE -le 8 ]]; then
  echo "Deleting Stage 8 (CM Search) for Round $ROUND..."
  for d in "$ROOT"/CLUSTER/$ROUND.*.cluster; do
    [[ -d "$d" ]] || continue
    rm -f "$d/cmsearch.DONE" "$d/expand.DONE"
    rm -rf "$d/CMSEARCH" "$d/SGE_log_expand"
  done
fi

########################################
# Delete Stage 7 (model) for current round
########################################

if [[ $STAGE -le 7 ]]; then
  echo "Deleting Stage 7 (Build Model) for Round $ROUND..."
  for d in "$ROOT"/CLUSTER/$ROUND.*.cluster; do
    [[ -d "$d" ]] || continue
    rm -f "$d/model.DONE" "$d/model_build.DONE"
    rm -rf "$d/MODEL" "$d/SGE_log_MSA" "$d/SUBTREES" "$d/TREE"
    rm -rf "$d/maligs.loc_new" "$d/stage7.cmd"
  done
fi

########################################
# Delete Stage 6 (LocARNA align) for current round
########################################

if [[ $STAGE -le 6 ]]; then
  echo "Deleting Stage 6 (LocARNA Align) for Round $ROUND..."
  for d in "$ROOT"/CLUSTER/$ROUND.*.cluster; do
    [[ -d "$d" ]] || continue
    rm -f "$d/pp.DONE" "$d/cluster.DONE"
    rm -rf "$d/pp" "$d/dp" "$d/bestSubtrees"
  done
fi

########################################
# Delete Stage 5 (NSPDK Clustering) for current round
########################################

if [[ $STAGE -le 5 ]]; then
  echo "Deleting Stage 5 (NSPDK Clustering) for Round $ROUND..."
  # Delete cluster directories entirely
  find "$ROOT/CLUSTER" -maxdepth 1 -name "$ROUND.*.cluster" -type d -exec rm -rf {} +
  rm -f "$ROOT/SVECTOR/data.svector.fast_cluster.$ROUND.DONE"
  rm -f "$ROOT/SVECTOR/data.svector.$ROUND.fast_cluster"
  rm -f "$ROOT/SVECTOR/data.svector.$ROUND.fast_cluster_sim"
fi

########################################
# Delete Stage 3 (NSPDK Vectors)
########################################

if [[ $STAGE -le 3 && $ROUND -eq 1 ]]; then
  echo "Deleting Stage 3 (NSPDK Vectors)..."
  rm -f "$ROOT/SVECTOR/svector.groups.DONE"
  find "$ROOT/SVECTOR" -name "*.svector.*" -not -name "*.blacklist.*" -delete 2>/dev/null
  rm -rf "$ROOT/SVECTOR/sge_log"
fi

########################################
# Delete Stage 2 (Graph Generation)
########################################

if [[ $STAGE -le 2 && $ROUND -eq 1 ]]; then
  echo "Deleting Stage 2 (Graph Generation)..."
  rm -f "$ROOT/GSPAN/gspan.DONE"
  find "$ROOT/GSPAN" -name "*.gspan.bz2" -delete 2>/dev/null
  rm -rf "$ROOT/GSPAN/sge_log"
fi

########################################
# Delete Stage 1 (Preprocessing)
########################################

if [[ $STAGE -le 1 && $ROUND -eq 1 ]]; then
  echo "Deleting Stage 1 (Preprocessing)..."
  rm -rf "$ROOT/FASTA"
fi

########################################
# Delete Stage 0 (Init)
########################################

if [[ $STAGE -le 0 && $ROUND -eq 1 ]]; then
  echo "Deleting Stage 0 (Init) — full reset..."
  rm -rf "$ROOT/FASTA" "$ROOT/GSPAN" "$ROOT/SVECTOR" "$ROOT/CLUSTER" "$ROOT/EVAL" "$ROOT/RESULTS" "$ROOT/SGE_LOG"
fi

########################################
# Clean up JOB.ERROR files everywhere
########################################

find "$ROOT" -name "JOB.ERROR" -delete 2>/dev/null

########################################
# Summary
########################################

echo ""
echo "========================================"
echo " Reset complete!"
echo " Pipeline is ready to start at:"
echo "   Round $ROUND, Stage $STAGE"
echo ""
echo " Start with:"
echo "   nohup bash $(dirname "$0")/run_pipeline.sh --root-dir $ROOT --config $ROOT/config --threads 10 > $ROOT/pipeline.log 2>&1 &"
echo "========================================"
